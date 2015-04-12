
# Module for working with curve files
# Curve files are csv deliniated with the following columns:
#    entry - the index in the file, starts at 0
#    lane  - id of the lane, starts at 0. Not necessarily ordered left-to-right
#    s     - distance along curve [m]
#    x     - x-pos in UTM [m]
#    y     - y-pos in UTM [m]
#    t     - heading [rad] [0-2pi]
#    k     - curvature [1/m]
#    k_d   - derivative of curvature [1/(m^2)]

module Curves

# ----------------------------------------------------------------

using DataFrames

# ----------------------------------------------------------------

export SIND, XIND, YIND, TIND, KIND, KDIND, CURVE_PT_VEC_LEN
export Curve, CurveSet
export fit_curve
export curve_at, closest_point_ind_to_curve, closest_point_ind_to_curve_guess
export closest_point_extind_to_curve, closest_point_extind_to_curve_guess
export pt_to_frenet_xy, pt_to_frenet_xyy, pt_to_roadway_xyy
export traj_to_frenet, extind_offset_from_end
export project_to_curves, nonexistant_curves, lane_rhs_distances
export lanewidths, interp_offset
export closest_point_to_curveset
export curveset_ordering
export arcdist_to_extind
export proj

# ----------------------------------------------------------------

const SIND  = 1
const XIND  = 2
const YIND  = 3
const TIND  = 4
const KIND  = 5
const KDIND = 6
const CURVE_PT_VEC_LEN = 6

immutable Curve	
	lane :: Int             # id of the lane, starts at 0. Not necessarily ordered left-to-right
	s    :: Vector{Float64} # distance along curve [m]
	x    :: Vector{Float64} # x-pos in UTM [m]
	y    :: Vector{Float64} # y-pos in UTM [m]
	t    :: Vector{Float64} # heading [rad] [0-2pi]
	k    :: Vector{Float64} # curvature [1/m]
	k_d  :: Vector{Float64} # derivative of curvature [1/(m^2)]
	
	function Curve()
		new(-1,Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])
	end
	function Curve(
		lane :: Int,
		s    :: Vector{Float64},
		x    :: Vector{Float64},
		y    :: Vector{Float64},
		t    :: Vector{Float64},
		k    :: Vector{Float64},
		k_d  :: Vector{Float64}
		)
		N = length(s)
		@assert(N == length(x) == length(y) == length(t) == length(k) == length(k_d))
		new(lane, s, x, y, t, k, k_d)
	end
end

typealias CurveSet Vector{Curve}

const THRESHOLD_TO_NOT_BOTHER_INTERPOLATING = 1e-4

# ----------------------------------------------------------------

include("splines.jl")

function load( filename::String; load_kdtrees::Bool = true )
	# LOAD: loads a curve file, returning an array of curves index by their curve id's + 1

	@assert( splitext(basename(filename))[2] == ".csv" )

	df = readtable(filename)

	lane_arr = convert(Array{Int},     df[:lane])
	s_arr    = convert(Array{Float64}, df[:s])
	x_arr    = convert(Array{Float64}, df[:x])
	y_arr    = convert(Array{Float64}, df[:y])
	t_arr    = convert(Array{Float64}, df[:t])
	k_arr    = convert(Array{Float64}, df[:k])
	kd_arr   = convert(Array{Float64}, df[:k_d])

	lanes = unique(lane_arr)

	# ensure lanes 0 -> max all exist
	@assert( isempty(setdiff(lanes, [0:maximum(lanes)])) )

	retval = Array(Curve, maximum(lanes+1))
	for laneid = 0 : maximum(lanes)
		lane_inds = find(x->x==laneid, lane_arr)
		curve = Curve(laneid, s_arr[lane_inds], x_arr[lane_inds], y_arr[lane_inds],
				              t_arr[lane_inds], k_arr[lane_inds], kd_arr[lane_inds])
		retval[laneid+1] = curve
	end
	
	return retval
end

function fit_curve(
	pts            :: AbstractMatrix{Float64}, # 2×n
	lane_id        :: Integer,
	desired_distance_between_samples :: Real;
	max_iterations :: Int = 50,
	epsilon        :: Float64 = 1e-4,
	n_intervals_in_arclen :: Int = 100
	)
	
	@assert(size(pts, 1) == 2)

	spline_coeffs = fit_cubic_spline(pts)

	L = calc_curve_length(spline_coeffs[1], spline_coeffs[2], n_intervals_per_segment=n_intervals_in_arclen)
	n = int(round(L/desired_distance_between_samples))+1

	s_arr = linspace(0.0,L,n)
	t_arr = calc_curve_param_given_arclen(spline_coeffs[1], spline_coeffs[2], s_arr, 
			curve_length=L, max_iterations=max_iterations, epsilon=epsilon, n_intervals_in_arclen=n_intervals_in_arclen)
	
	x_arr = sample_spline(spline_coeffs[1], t_arr)
	y_arr = sample_spline(spline_coeffs[2], t_arr)
	θ_arr = sample_spline_theta(spline_coeffs[1], spline_coeffs[2], t_arr)
	κ_arr = sample_spline_curvature(spline_coeffs[1], spline_coeffs[2], t_arr)
	κd_arr = sample_spline_derivative_of_curvature(spline_coeffs[1], spline_coeffs[2], t_arr)

	@assert(!any(s->isnan(s), s_arr))
	@assert(!any(s->isnan(s), x_arr))
	@assert(!any(s->isnan(s), y_arr))
	@assert(!any(s->isnan(s), θ_arr))

	Curve(lane_id, s_arr, x_arr, y_arr, θ_arr, κ_arr, κd_arr)
end


curve_at( curve::Curve, ind::Int) = [curve.s[ind], curve.x[ind], curve.y[ind], curve.t[ind], curve.k[ind], curve.k_d[ind]]
function curve_at( curve::Curve, extind::Float64 )

	@assert( extind >= 1.0 && extind <= length(curve) )

	ind_lo = ifloor(extind)
	t      = extind - ind_lo
	p_lo   = curve_at(curve, ind_lo)

	if t < THRESHOLD_TO_NOT_BOTHER_INTERPOLATING
		return p_lo
	elseif 1.0 - t < THRESHOLD_TO_NOT_BOTHER_INTERPOLATING
		return curve_at(curve, ind_lo+1)
	end

	p_hi = curve_at(curve, ind_lo+1)
	return lerp( p_lo, p_hi, t )
end

function closest_point_ind_to_curve( curve::Curve, x::Float64, y::Float64 )
	# find the index of the point closest to the curve
	# picks the lower one in a tie

	return indmin([norm([x - curve.x[i], y - curve.y[i]]) for i = 1 : length(curve)])
end
function closest_point_ind_to_curve_guess( curve::Curve, x::Real, y::Real, ind_guess::Int )
	# like closest_point_ind_to_curve, but starts from ind_guess and picks the local min

	distat = (i) -> norm([x - curve.x[i], y - curve.y[i]])
	start_dist = distat(ind_guess)

	up_ind_prev = ind_guess
	up_ind = min(ind_guess + 1, length(curve))
	dist_up_prev = start_dist
	dist_up = distat(up_ind)
	while dist_up < dist_up_prev
		up_ind_prev = up_ind
		up_ind = min(up_ind + 1, length(curve))
		dist_up_prev = dist_up
		dist_up = distat(up_ind)
	end

	lo_ind_prev = ind_guess
	lo_ind = max(ind_guess - 1, 1)
	dist_lo_prev = start_dist
	dist_lo = distat(lo_ind)
	while dist_lo < dist_lo_prev
		lo_ind_prev = lo_ind
		lo_ind = max(lo_ind - 1, 1)
		dist_lo_prev = dist_lo
		dist_lo = distat(lo_ind)
	end

	i = indmin([dist_lo_prev, start_dist, dist_up_prev])
	return [lo_ind_prev, ind_guess, up_ind_prev][i]
end
function closest_point_ind_to_curve( curve::Curve, s::Float64 )
	# find the index of the s-value closest to the curve. 
	# picks the lower one in a tie

	return _binary_search( curve.s, s )
end
function closest_point_ind_to_curve_guess( curve::Curve, s::Real, ind_guess::Int )
	# like closest_point_ind_to_curve, but starts from ind_guess and picks the local min

	distat = (i) -> abs(curve.s[i] - s)
	start_dist = distat(ind_guess)

	up_ind_prev = ind_guess
	up_ind = min(ind_guess + 1, length(curve))
	dist_up_prev = start_dist
	dist_up = distat(up_ind)
	while dist_up < dist_up_prev
		up_ind_prev = up_ind
		up_ind = min(up_ind + 1, length(curve))
		dist_up_prev = dist_up
		dist_up = distat(up_ind)
	end

	lo_ind_prev = ind_guess
	lo_ind = max(ind_guess - 1, 1)
	dist_lo_prev = start_dist
	dist_lo = distat(lo_ind)
	while dist_lo < dist_lo_prev
		lo_ind_prev = lo_ind
		lo_ind = max(lo_ind - 1, 1)
		dist_lo_prev = dist_lo
		dist_lo = distat(lo_ind)
	end

	i = indmin([dist_lo_prev, start_dist, dist_up_prev])
	return [lo_ind_prev, ind_guess, up_ind_prev][i]
end

function closest_point_extind_to_curve( curve::Curve, x::Real, y::Real )

	ind = closest_point_ind_to_curve(curve, x, y)

	pt = [x, y]

	# interpolate farther
	if ind > 1 && ind < length(curve)
		t_lo = proj( [curve.x[ind-1], curve.y[ind-1]], [curve.x[ind], curve.y[ind]], pt )
		t_hi = proj( [curve.x[ind], curve.y[ind]], [curve.x[ind+1], curve.y[ind+1]], pt )

		p_lo = lerp( [curve.x[ind-1], curve.y[ind-1]], [curve.x[ind], curve.y[ind]], t_lo )
		p_hi = lerp( [curve.x[ind], curve.y[ind]], [curve.x[ind+1], curve.y[ind+1]], t_hi )

		d_lo = norm( p_lo - pt )
		d_hi = norm( p_hi - pt )

		return d_lo < d_hi ? ind-1+t_lo : ind+t_hi
	elseif ind == 1
		return ind + proj( [curve.x[ind], curve.y[ind]], [curve.x[ind+1], curve.y[ind+1]], pt )
	else
		return ind -1 + proj( [curve.x[ind-1], curve.y[ind-1]], [curve.x[ind], curve.y[ind]], pt )
	end
end
function closest_point_extind_to_curve_guess( curve::Curve, x::Real, y::Real, extind_guess::Real )

	ind = closest_point_ind_to_curve_guess(curve, x, y, int(extind_guess))

	pt = [x, y]

	# interpolate farther
	if ind > 1 && ind < length(curve)
		t_lo = proj( [curve.x[ind-1], curve.y[ind-1]], [curve.x[ind], curve.y[ind]], pt )
		t_hi = proj( [curve.x[ind], curve.y[ind]], [curve.x[ind+1], curve.y[ind+1]], pt )

		p_lo = lerp( [curve.x[ind-1], curve.y[ind-1]], [curve.x[ind], curve.y[ind]], t_lo )
		p_hi = lerp( [curve.x[ind], curve.y[ind]], [curve.x[ind+1], curve.y[ind+1]], t_hi )

		d_lo = norm( p_lo - pt )
		d_hi = norm( p_hi - pt )

		return d_lo < d_hi ? ind-1+t_lo : ind+t_hi
	elseif ind == 1
		return ind + proj( [curve.x[ind], curve.y[ind]], [curve.x[ind+1], curve.y[ind+1]], pt )
	else
		return ind -1 + proj( [curve.x[ind-1], curve.y[ind-1]], [curve.x[ind], curve.y[ind]], pt )
	end
end
function closest_point_extind_to_curve( curve::Curve, s::Real )

	ind = closest_point_ind_to_curve( curve, s )

	pt = [s]

	# interpolate farther
	if ind > 1 && ind < length(curve)
		t_lo = proj( [curve.s[ind-1]], [curve.s[ind]], pt)
		t_hi = proj( [curve.s[ind]], [curve.s[ind+1]], pt)

		p_lo = lerp( [curve.s[ind-1]], [curve.s[ind]], t_lo)
		p_hi = lerp( [curve.s[ind]], [curve.s[ind+1]], t_hi)

		d_lo = norm( p_lo - pt )
		d_hi = norm( p_hi - pt )

		return d_lo < d_hi ? ind-1+t_lo : ind+t_hi
	elseif ind == 1
		return ind + proj( [curve.s[ind]], [curve.s[ind+1]], pt )
	else
		return ind-1 + proj( [curve.s[ind-1]], [curve.s[ind]], pt )
	end
end
function closest_point_extind_to_curve_guess( curve::Curve, s::Real, extind_guess::Real )

	ind = closest_point_ind_to_curve_guess( curve, s, int(extind_guess) )

	pt = [s]

	# interpolate farther
	if ind > 1 && ind < length(curve)
		t_lo = proj( [curve.s[ind-1]], [curve.s[ind]], pt)
		t_hi = proj( [curve.s[ind]], [curve.s[ind+1]], pt)

		p_lo = lerp( [curve.s[ind-1]], [curve.s[ind]], t_lo)
		p_hi = lerp( [curve.s[ind]], [curve.s[ind+1]], t_hi)

		d_lo = norm( p_lo - pt )
		d_hi = norm( p_hi - pt )

		return d_lo < d_hi ? ind-1+t_lo : ind+t_hi
	elseif ind == 1
		return ind + proj( [curve.s[ind]], [curve.s[ind+1]], pt )
	else
		return ind-1 + proj( [curve.s[ind-1]], [curve.s[ind]], pt )
	end
end

function pt_to_frenet_xy( curvept::Vector{Float64}, x::Real, y::Real)

	@assert(length(curvept) == CURVE_PT_VEC_LEN)

	s = curvept[SIND]
	d = norm(curvept[[XIND,YIND]] - [x, y])
	t = curvept[TIND]

	dyaw = mod2pi( atan2( y-curvept[YIND], x-curvept[XIND] ) - t )

	mod2pi2 = x->begin
		val = mod2pi(x)
		if val > pi
			val -= 2pi
		end
		return val
	end

	on_left_side = abs(mod2pi2(dyaw - pi/2)) < abs(mod2pi2(dyaw - 3pi/2))
	d *= on_left_side ? 1.0 : -1.0 # left side is positive, right side is negative

	return (s,d)
end
function pt_to_frenet_xy( curve::Curve, x::Real, y::Real, guess::Nothing=nothing )

	extind = closest_point_extind_to_curve(curve, x, y)

	if extind == 1.0 || extind == length(curve)
		warn("Curve.pt_to_frenet: results are unreliable at curve start and end")
	end

	ptoncurve = curve_at( curve, extind )
	s,d = pt_to_frenet_xy(ptoncurve, x, y)
	return (s,d,extind)
end
function pt_to_frenet_xy( curve::Curve, x::Real, y::Real, guess::Real )

	extind = closest_point_extind_to_curve_guess(curve, x, y, guess)

	if extind == 1.0 || extind == length(curve)
		warn("Curve.pt_to_frenet: results are unreliable at curve start and end")
	end

	ptoncurve = curve_at( curve, extind )
	s,d = pt_to_frenet_xy(ptoncurve, x, y)
	return (s,d,extind)
end
function pt_to_frenet_xyy( curvept::Vector{Float64}, x::Real, y::Real, yaw::Real )

	@assert(length(curvept) == CURVE_PT_VEC_LEN)

	s = curvept[SIND]
	d = norm(curvept[[XIND,YIND]] - [x, y])
	t = curvept[TIND]

	dyaw = mod2pi( atan2( y-curvept[YIND], x-curvept[XIND] ) - t )

	mod2pi2 = x->begin
		val = mod2pi(x)
		if val > pi
			val -= 2pi
		end
		return val
	end

	on_left_side = abs(mod2pi2(dyaw - pi/2)) < abs(mod2pi2(dyaw - 3pi/2))
	d *= on_left_side ? 1.0 : -1.0 # left side is positive, right side is negative

	yaw = mod2pi(yaw-t)
	if yaw > pi
		yaw -= 2pi
	end

	return (s,d,yaw)
end
function pt_to_frenet_xyy( curve::Curve, x::Real, y::Real, yaw::Real, guess::Nothing=nothing )

	extind = closest_point_extind_to_curve(curve, x, y)

	if extind == 1.0 || extind == length(curve)
		warn("Curve.pt_to_frenet: results are unreliable at curve start and end")
	end

	ptoncurve = curve_at( curve, extind )
	s,d,yaw = pt_to_frenet_xyy(ptoncurve, x, y, yaw)
	return (s,d,yaw,extind)
end
function pt_to_frenet_xyy( curve::Curve, x::Real, y::Real, yaw::Real, guess::Real )

	extind = closest_point_extind_to_curve_guess(curve, x, y, guess)

	if extind == 1.0 || extind == length(curve)
		warn("Curve.pt_to_frenet: results are unreliable at curve start and end")
	end

	ptoncurve = curve_at( curve, extind )
	s,d,yaw = pt_to_frenet_xyy(ptoncurve, x, y, yaw)
	return (s,d,yaw,extind)
end

function extind_offset_from_end( curve::Curve, offset::Float64 )
	# given a curve, find the extind values that are an offset distance from start and end
	# lo -> extind such that curve_at( curve, lo ).s = offset
	# hi -> extind such that curve_at( curve, hi ).s = max(curve.s) - offset

	lo = closest_point_extind_to_curve_guess( curve, curve.s[1]   + offset, 1.0          )::Float64
	hi = closest_point_extind_to_curve_guess( curve, curve.s[end] - offset, length(curve))::Float64
	(lo, hi)
end

function project_to_curves(
	curveset :: CurveSet,
	x        :: Real,
	y        :: Real,
	yaw      :: Real
	)
	
	# project the point to each curve, returning a list of projections

	pts = Array(Tuple, length(curveset)) # (x,y,yaw,extind)
	for (i,curve) in enumerate(curveset)
		pts[i] = pt_to_frenet_xyy(curve, x, y, yaw)
		if pts[i][end] == 1.0 || pts[i][end] == length(curve)
			warn("Curve.pt_to_frenet: results are unreliable at curve start and end")
		end
	end
	pts
end
function project_to_curves(
	curveset :: CurveSet,
	x        :: Real,
	y        :: Real,
	yaw      :: Real,
	guesses  :: Vector
	)
	
	# project the point to each curve, returning a list of projections

	pts = Array(Tuple, length(curveset)) # (x,y,yaw,extind)
	for (i,curve) in enumerate(curveset)
		pts[i] = pt_to_frenet_xyy(curve, x, y, yaw, guesses[i])
		if pts[i][end] == 1.0 || pts[i][end] == length(curve)
			warn("Curve.pt_to_frenet: results are unreliable at curve start and end")
		end
	end
	pts
end
function lanewidths( d_lane::Vector{Float64} )
	if isempty(d_lane)
		return Float64[]
	elseif length(d_lane) == 1
		return [0.0]
	elseif length(d_lane) == 2
		return [d_lane[2] - d_lane[1]]
	else
		return d_lane[2:end] - [0, d_lane[2:end-1]]
	end
end

function pt_to_roadway_xyy(
	curveset :: CurveSet, # assumes it is ordered
	x        :: Real, 
	y        :: Real, 
	θ        :: Real
	)

	# project the point to each curve
	# if the point is between two curves, weight the actual distance
	# returns  ((x,y,yaw), guess_ret)

	pts = [pt_to_frenet_xyy(curve, x, y, θ, nothing) for curve in curveset]
	d_lane = lane_rhs_distances(curveset, pts)
	temppt = interp_offset(pts, d_lane)
	((temppt[1],temppt[2],temppt[3]), map(pt->pt[4], pts))
end
function pt_to_roadway_xyy(
	curveset :: CurveSet, # assumes it is ordered
	x        :: Real, 
	y        :: Real, 
	θ        :: Real,
	guesses  :: Vector{Float64}
	)

	# project the point to each curve
	# if the point is between two curves, weight the actual distance
	# returns  ((x,y,yaw), guess_ret)

	pts = [pt_to_frenet_xyy(curve, x, y, θ, guesses) for curve in curveset]
	d_lane = lane_rhs_distances(curveset, pts)
	temppt = interp_offset(pts, d_lane)
	((temppt[1],temppt[2],temppt[3]), map(pt->pt[4], pts))
end

function arcdist_to_extind(curve::Curve, s::Real; clamp::Bool = true)
	# find the closest extind to the given arc-distance along the curve
	# will clamp to a valid extind if clamp is set [default]
	# otherwise returns -1.0 if the arc-distance is out of range

	# NOTE(tim): this assumes that the curve is ordered

	arr = curve.s
    lo, hi = 1, length(arr)

    if s < arr[lo]
        return clamp ? 1.0 : -1.0
    elseif s > arr[hi]
        return clamp ? float64(hi) : -1.0
    end

    while lo < hi-1
        mid::Int = div(lo + hi, 2)
        v = arr[mid]
        if isapprox(v, s)
            return float64(mid)
        elseif v < s
            lo = mid
        elseif s < v
            hi = mid
        end
    end
    slo, shi = arr[lo], arr[hi]
    lo + (s - slo) / (shi - slo)
end

# ----------------------------------------------------------------

function closest_point_to_curveset( curveset::CurveSet, x::Real, y::Real )

	best_dist     = Inf
	best_curveind = 0
	best_pt       = []
	best_extind   = 0.0
	for i = 1 : length(curveset)

		extind = closest_point_extind_to_curve( curveset[i], x, y )
		pt = curve_at(curveset[i], extind)[[XIND,YIND]]
		dist = norm([x,y]-pt)
		if dist < best_dist
			best_dist, best_curveind, best_pt, best_extind = (dist, i, pt, extind)
		end
	end

	return (best_pt, best_curveind, best_extind)
end

function curveset_ordering( curveset::CurveSet )

	# compute the right -> left ordering of the curves

	ncurves = length(curveset)

	mid_ind = int(length(curveset[1])/2) # pick the middle of the first curve

	pts = Array(Vector,ncurves)
	pts[1] = curve_at(curveset[1], mid_ind)[[XIND,YIND]]
	for i = 2 : ncurves
		extind = closest_point_extind_to_curve(curveset[i], pts[1][1], pts[1][2])
		pts[i] = curve_at(curveset[i], extind)[[XIND,YIND]] - pts[1]
	end

	# we now have a set of points, project them onto their cross-section axis
	along = curve_at(curveset[1], mid_ind + 1.0)[[XIND,YIND]] - pts[1] # point one tick along the curve
	pts[1] = [0.0,0.0]
	along = along ./ norm(along) # unit vector
	R = [ 0.0 -1.0;
	      1.0  0.0] # 90-degree rotation to left
	axis = R*along # axis pointing towards positive s-direction

	projs = zeros(ncurves)
	for i = 1 : ncurves
		projs[i] = dot(axis, pts[i])
	end

	# println("axis: ",  axis)
	# println("along: ", along)
	# res = sortperm(projs)
	# for i = 1 : ncurves
	# 	@printf("%d: (%10.3f %10.3f) @ %10.3f\n", i, pts[i][1], pts[i][2], projs[i])
	# end
	# println("retval: ", res)

	return sortperm(projs) # returns lane ordering, lowest is rightmost lane
end

# ----------------------------------------------------------------

Base.length( curve::Curve ) = length(curve.x)
function lerp( P0::Array, P1::Array, t::Real )

	return P0 + (P1 - P0)*t
end
function proj( P0::Array, P1::Array, Q::Array )

	b = P1 - P0
	a = Q  - P0

	c = b * dot(a,b) / dot(b,b)

	t = 0.0
	if b[1] != 0.0
		t = c[1] / b[1]
	else
		t = c[2] / b[2]
	end

	return max(0.0, min(t, 1.0))
end

function _binary_search( arr::Vector, f::Function, a::Int=1, b::Int=length(arr) )

	# INPUT: Function f, endpoint values a, b, tolerance TOL, maximum iterations NMAX
	# CONDITIONS: a < b, either f(a) < 0 and f(b) > 0 or f(a) > 0 and f(b) < 0
	# OUTPUT: value which differs from a root of f(x)=0 by less than TOL

	n = 1
	while true

		if b == a
			return a
		elseif b == a + 1
			return f(b) < -f(a) ? b : a
		end

		c = div(a+b,2)
		fc = f(arr[c])
		if isapprox(fc, 0.0) 
			return c
		elseif fc < 0.0
			a = c
		else
			b = c
		end
	end
end
function _binary_search( arr::Vector{Float64}, target::Float64, a::Int=1, b::Int=length(arr) )

	# INPUT: Function f, endpoint values a, b, tolerance TOL, maximum iterations NMAX
	# CONDITIONS: a < b, either f(a) < 0 and f(b) > 0 or f(a) > 0 and f(b) < 0
	# OUTPUT: value which differs from a root of f(x)=0 by less than TOL

	n = 1
	while true
		if b == a
			return a
		elseif b == a + 1
			return arr[b]-target < target-arr[a] ? b : a
		end

		c = div( a+b,2 )
		fc = arr[c] - target

		if isapprox(fc, 0.0) 
			return c
		elseif fc < 0.0
			a = c
		else
			b = c
		end
	end
end

# ----------------------------------------------------------------

end # end module