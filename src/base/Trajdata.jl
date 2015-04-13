# Traj Data

# code for the use of trajectory data files
# the DataFrame version contains raw info obtained from csv files
# the PrimaryDataset type contains processed data

module Trajdata

using DataFrames
using DataArrays

import Base: start, next, done, get

const NUMBER_REGEX = r"(-)?(0|[1-9]([\d]*))(\.[\d]*)?((e|E)(\+|-)?[\d]*)?"

export PrimaryDataset
export CARIND_EGO, CARID_EGO
export CONTROL_STATUS_START, CONTROL_STATUS_ACTIVE, CONTROL_STATUS_ACTIVE_REINT
export CONTROL_STATUS_CONTROL_ACTIVE, CONTROL_STATUS_AUTO, CONTROL_STATUS_HANDOVER_TO_DRIVER
export CONTROL_STATUS_HANDOVER_TO_SYSTEM, CONTROL_STATUS_PREPARE_ACTUATORS, CONTROL_STATUS_READY
export CONTROL_STATUS_WAITING, CONTROL_STATUS_FAILURE, CONTROL_STATUS_INIT, CONTROL_STATUS_PASSIVE
export validfind_inbounds, frameind_inbounds, nframeinds, nvalidfinds
export sec_per_frame, closest_frameind, closest_validfind
export carid2ind, carind2id, validfind2frameind, frameind2validfind, jumpframe
export getc, gete, gete_validfind, setc!, sete!, idinframe, indinframe, getc_symb, getc_has, setc_has!
export carind_exists, carid_exists, frame2frameind, get_max_num_cars, get_max_carind, get_maxcarind
export get_carids, get_valid_frameinds, load_trajdata, export_trajdata
export remove_car_from_frame!, remove_cars_from_frame!, remove_car!, find_slot_for_car!
export ncars_in_frame, add_car_slot!, frames_contain_carid, extract_continuous_tracks, spanning_validfinds		
export isnum

# ------------------------

#=
CONTROL STATUS VARIABLE FOR BOSCH VEHICLE IS CENSORED
=#

const CARIND_EGO                        = -1
const CARID_EGO                         = -1

type PrimaryDataset

	df_ego_primary     :: DataFrame              # [frameind,  :feature] -> value
	df_other_primary   :: DataFrame              # [validfind, :feature] -> value, missing features are NA
	dict_trajmat       :: Dict{Uint32,DataFrame} # [carid] -> trajmat
	dict_other_idmap   :: Dict{Uint32,Uint16}    # [carid] -> matind
	mat_other_indmap   :: Matrix{Int16}          # [validfind, matind] -> carind, -1 if not present
	ego_car_on_freeway :: BitArray{1}            # [frameind] -> bool
	validfind2frameind :: Vector{Int32}          # [validfind] -> frameind
	frameind2validfind :: Vector{Int32}          # [frameind] -> validfind, 0 if not applicable
	maxcarind          :: Int16                  # max ind of car in frame

	PrimaryDataset(
		df_ego_primary     :: DataFrame, 
		df_other_primary   :: DataFrame,
		dict_trajmat       :: Dict{Uint32,DataFrame},
		dict_other_idmap   :: Dict{Uint32,Uint16},
		mat_other_indmap   :: Matrix{Int16},
		ego_car_on_freeway :: BitArray{1}
		) = begin

		n_frames = size(df_ego_primary, 1)
		n_valid  = size(df_other_primary, 1)

		validfind2frameind = zeros(Int32, n_valid)
		frameind2validfind = zeros(Int32, n_frames)

		validfind = 0
		for frameind = 1 : n_frames
			if ego_car_on_freeway[frameind]
				validfind += 1
				validfind2frameind[validfind] = frameind
				frameind2validfind[frameind] = validfind
			end
		end

		maxcarind = -1
		while haskey(df_other_primary, symbol(@sprintf("id_%d", maxcarind+1)))
			maxcarind += 1
		end

		new(df_ego_primary, 
			df_other_primary,
			dict_trajmat,
			dict_other_idmap,
			mat_other_indmap, 
			ego_car_on_freeway, 
			validfind2frameind,
			frameind2validfind,
			maxcarind)
	end
end

export IterOtherCarindsInFrame, IterAllCarsInFrame

type IterOtherCarindsInFrame
	# An iterator over other cars in a frame
	# iterators from carind 0 to maxcarind
	# does not return ego carind

	# the state is the carind

	pdset     :: PrimaryDataset
	validfind :: Int
end
type IterAllCarindsInFrame
	# An iterator over all cars in a frame
	# iterators from carind -1 to maxcarind
	# does not return ego carind

	# the state is the carind

	pdset     :: PrimaryDataset
	validfind :: Int
end

start(I::IterOtherCarindsInFrame) = 0
function done(I::IterOtherCarindsInFrame, carind::Int)
	# NOTE(tim): the state is the carind
	!validfind_inbounds(I.pdset, I.validfind) || !indinframe(I.pdset, carind, I.validfind)
end
next(I::IterOtherCarindsInFrame, carind::Int) = (carind, carind+1)

start(I::IterAllCarindsInFrame) = -1
function done(I::IterAllCarindsInFrame, carind::Int)
	# NOTE(tim): the state is the carind
	if !validfind_inbounds(I.pdset, I.validfind)
		return true
	end
	if carind == CARIND_EGO
		return false
	end
	!indinframe(I.pdset, carind, I.validfind)
end
next(I::IterAllCarindsInFrame, carind::Int) = (carind, carind+1)

nframeinds( trajdata::DataFrame ) = size(trajdata, 1)
nframeinds(  pdset::PrimaryDataset ) = length(pdset.frameind2validfind)
nvalidfinds( pdset::PrimaryDataset ) = length(pdset.validfind2frameind)

function carid2ind( trajdata::DataFrame, carid::Integer, frameind::Integer)
	@assert(haskey(trajdata, symbol(@sprintf("has_%d", carid))))
	int(trajdata[symbol(@sprintf("has_%d", carid))][frameind])
end
function carid2ind( pdset::PrimaryDataset, carid::Integer, validfind::Integer )
	# returns -1 if not in the frame
	@assert(haskey(pdset.dict_other_idmap, carid))
	matind = pdset.dict_other_idmap[carid]
	pdset.mat_other_indmap[validfind, matind]
end
function carind2id( trajdata::DataFrame, carind::Integer, frameind::Integer)
	if carind == CARIND_EGO
		return CARID_EGO
	end
	@assert(haskey(trajdata, symbol(@sprintf("id_%d", carind))))
	trajdata[frameind, symbol(@sprintf("id_%d", carind))]
end
function carind2id( pdset::PrimaryDataset, carind::Integer, validfind::Integer )
	if carind == CARIND_EGO
		return CARID_EGO
	end
	@assert(haskey(pdset.df_other_primary, symbol(@sprintf("id_%d", carind))))
	pdset.df_other_primary[validfind, symbol(@sprintf("id_%d", carind))]
end

validfind_inbounds( pdset::PrimaryDataset, validfind::Integer ) = 1 <= validfind <= length(pdset.validfind2frameind)
frameind_inbounds( pdset::PrimaryDataset, frameind::Integer ) = 1 <= frameind <= length(pdset.frameind2validfind)
function frameind2validfind( pdset::PrimaryDataset, frameind::Integer )
	# returns 0 if it does not exist
	if !frameind_inbounds(pdset, frameind)
		return 0
	end
	pdset.frameind2validfind[frameind]
end
function validfind2frameind( pdset::PrimaryDataset, validfind::Integer )
	# returns 0 if it does not exist
	if !validfind_inbounds(pdset, validfind)
		return 0
	end
	pdset.validfind2frameind[validfind]
end
function jumpframe( pdset::PrimaryDataset, validfind::Integer, delta::Integer )
	# get the validfind that has a frameind of + delta
	# returns 0 if it does not exist
	frameind = validfind2frameind(pdset, validfind)
	jump_frameind = frameind + delta
	frameind2validfind(pdset, jump_frameind)
end

function sec_per_frame( pdset::PrimaryDataset )
	# estimate the difference between frames assuming that is it more or less consistent
	# across the entire sequence
	nframes = nframeinds(pdset)
	t0 = gete(pdset, :time, 1)::Float64
	tf = gete(pdset, :time, nframes)::Float64

	(tf - t0) / (nframes-1)
end
function closest_frameind( pdset::PrimaryDataset, time::Float64 )
	# returns the frameind that is closest to the given time
	# uses the fact that Δt is more or less consistent

	N = nframeinds(pdset)

	Δt = sec_per_frame(pdset)
	t0 = gete(pdset, :time, 1)::Float64
	find = clamp(int((time - t0) / Δt + 1), 1, N)
	t  = gete(pdset, :time, find)::Float64

	tp = find < N ? gete(pdset, :time, find + 1)::Float64 : t
	if abs(tp - time) < abs(t-time)
		t, find = tp, find+1
		tp = find < N ? gete(pdset, :time, find + 1)::Float64 : t
		while abs(tp - time) < abs(t-time)
			t, find = tp, find+1
			tp = find < N ? gete(pdset, :time, find + 1)::Float64 : t
		end
		return find
	end

	tn = find > 1 ? gete(pdset, :time, find - 1)::Float64 : t
	if abs(tn - time) < abs(t-time)
		t, find = tn, find-1
		tn = find > 1 ? gete(pdset, :time, find - 1)::Float64 : t
		while abs(tn - time) < abs(t-time)
			t, find = tn, find-1
			tn = find > 1 ? gete(pdset, :time, find - 1)::Float64 : t
		end
		return find
	end

	find
end
function closest_validfind( pdset::PrimaryDataset, time::Float64 )
	# returns the validfind that is closest	to the given time
	find = closest_frameind(pdset, time)
	vind = frameind2validfind(pdset, find)
	if vind != 0
		return vind
	end

	if find < N
		find_forward = find + 1
		while frameind2validfind(pdset, find_forward) == 0
			find_forward += 1
			if find_forward ≥ N
				break
			end
		end
		fp = frameind2validfind(pdset, find_forward) == 0 ? 0 : find_forward
	else
		fp = 0
	end

	isvalid_n = false
	if find > 1
		find_back = find - 1
		while frameind2validfind(pdset, find_back) == 0
			find_back -= 1
			if find_back ≤ 1
				break
			end
		end
		fn = frameind2validfind(pdset, find_back) == 0 ? 0 : find_back
	else
		fn = 0
	end

	
	@assert(fn != fp)

	if fn == 0
		return fp
	elseif fp == 0
		return fn
	end

	tp = gete(pdset, :time, fp)::Float64
	tn = gete(pdset, :time, fn)::Float64

	dp = abs(tp - time)
	dn = abs(tn - time)

	return dp ≤ dn ? fp : fn
end


getc_symb( str::String, carind::Integer ) = symbol(@sprintf("%s_%d", str, carind))
function getc( trajdata::DataFrame, str::String, carind::Integer, frameind::Integer )
	@assert(haskey(trajdata, symbol(@sprintf("id_%d", carind))))
	trajdata[frameind, symbol(@sprintf("%s_%d", str, carind))]
end
function getc( pdset::PrimaryDataset, str::String, carind::Integer, validfind::Integer )
	@assert(haskey(pdset.df_other_primary, symbol(@sprintf("id_%d", carind))))
	pdset.df_other_primary[validfind, symbol(@sprintf("%s_%d", str, carind))]
end
function gete( pdset::PrimaryDataset, sym::Symbol, frameind::Integer )
	pdset.df_ego_primary[frameind, sym]
end
function gete_validfind( pdset::PrimaryDataset, sym::Symbol, validfind::Integer )
	pdset.df_ego_primary[validfind2frameind(pdset, validfind), sym]
end
function get(pdset::PrimaryDataset, str::String, carind::Integer, validfind::Integer)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		return gete(pdset, symbol(str), frameind)
	end
	getc(pdset, str, carind, validfind)
end
function get(pdset::PrimaryDataset, sym::Symbol, str::String, carind::Integer, frameind::Integer, validfind::Integer)
	# a more efficient version
	carind == CARIND_EGO ? gete(pdset, sym, frameind) : getc(pdset, str, carind, validfind)
end
function setc!( trajdata::DataFrame, str::String, carind::Integer, frameind::Integer, val )
	@assert(haskey(trajdata, symbol(@sprintf("id_%d", carind))))
	trajdata[symbol(@sprintf("%s_%d", str, carind))][frameind] = val
end
function sete!( pdset::PrimaryDataset, sym::Symbol, frameind::Integer, val )
	@assert(1 <= frameind && frameind <= length(pdset.frameind2validfind))
	pdset.df_ego_primary[frameind, sym] = val
end
function setc!( pdset::PrimaryDataset, str::String, carind::Integer, validfind::Integer, val )
	@assert(haskey(pdset.df_other_primary, symbol(@sprintf("id_%d", carind))))
	pdset.df_other_primary[validfind, symbol(@sprintf("%s_%d", str, carind))] = val
end
function idinframe( pdset::PrimaryDataset, carid::Integer, validfind::Integer )
	if !haskey(pdset.dict_other_idmap, carid)
		return false
	end
	if validfind < 1 || validfind > size(pdset.mat_other_indmap, 1)
		return false
	end
	matind = pdset.dict_other_idmap[carid]
	pdset.mat_other_indmap[validfind, matind] != -1
end
function indinframe( pdset::PrimaryDataset, carind::Integer, validfind::Integer )
	@assert(carind >= 0 && carind <= pdset.maxcarind)
	id_ = pdset.df_other_primary[validfind, symbol(@sprintf("id_%d", carind))]
	!isa(id_, NAtype)
end
function getc_has( trajdata::DataFrame, carid::Integer, frameind::Integer )
	@assert(haskey(trajdata, symbol(@sprintf("has_%d", carid))))
	return int(trajdata[symbol(@sprintf("has_%d", carid))][frameind])
end
function getc_has( trajdata::DataFrame, carid::Integer )
	@assert(haskey(trajdata, symbol(@sprintf("has_%d", carid))))
	return int(trajdata[symbol(@sprintf("has_%d", carid))])
end
function setc_has!( trajdata::DataFrame, carid::Integer, frameind::Integer, carind::Integer )
	@assert(haskey(trajdata, symbol(@sprintf("has_%d", carid))))
	@assert(carind >= -1)
	trajdata[symbol(@sprintf("has_%d", carid))][frameind] = carind
end
function get_max_num_cars( trajdata::DataFrame )
	i = -1
	while haskey(trajdata, symbol(@sprintf("id_%d", i+1)))
		i += 1
	end
	return i
end
function get_max_carind( trajdata::DataFrame )
	i = -1
	while haskey(trajdata, symbol(@sprintf("id_%d", i+1)))
		i += 1
	end
	return i
end
get_maxcarind( pdset::PrimaryDataset ) = pdset.maxcarind
function get_carids( trajdata::DataFrame )
	# obtain a set of car ids ::Set{Int}
	carids = Set{Int}()
	for name in names(trajdata)
		str = string(name)
		if ismatch(r"^has_(\d)+$", str)
			m = match(r"(\d)+", str)
			id = int(m.match)
			push!(carids, id)
		end
	end
	carids
end
get_carids( pdset::PrimaryDataset ) = Set(keys(pdset.dict_other_idmap))
function carind_exists( trajdata::DataFrame, carind::Integer, frameind::Integer )
	@assert( carind >= 0 )
	@assert( haskey(trajdata, symbol(@sprintf("id_%d", carind))) )
	!isa(trajdata[frameind, symbol(@sprintf("id_%d", carind))], NAtype)
end
function carid_exists( trajdata::DataFrame, carid::Integer, frameind::Integer )

	return getc_has( trajdata, carid, frameind ) != -1
end

function get_valid_frameinds( pdset::PrimaryDataset ) 

	all_inds = [1:size(pdset.df_ego_primary,1)]
	all_inds[pdset.ego_car_on_freeway]
end
frame2frameind( trajdata::DataFrame, frame::Integer ) = findfirst(x->x == frame, array(trajdata[:frame], -999))

function remove_car_from_frame!( pdset::PrimaryDataset, carind::Integer, validfind::Integer )
	@assert( carind >= 0 && carind <= pdset.maxcarind )
	@assert(validfind >= 0 && validfind <= length(pdset.validfind2frameind))

	# removes the car with the specified carind from the frame, shifting the other cars down

	# pull the carid from this frame
	id = getc( pdset, "id", carind, validfind )
	mat_ind = pdset.dict_other_idmap[id]
	pdset.mat_other_indmap[validfind, mat_ind] = -1

	loop_cind = carind
	while loop_cind < pdset.maxcarind

		if !indinframe(pdset, loop_cind, validfind)
			break
		end

		setc!( pdset, "id",        loop_cind, validfind, getc( pdset, "id",        loop_cind+1, validfind ))
		setc!( pdset, "posFx",     loop_cind, validfind, getc( pdset, "posFx",     loop_cind+1, validfind ))
		setc!( pdset, "posFy",     loop_cind, validfind, getc( pdset, "posFy",     loop_cind+1, validfind ))
		setc!( pdset, "posFyaw",   loop_cind, validfind, getc( pdset, "posFyaw",   loop_cind+1, validfind ))
		setc!( pdset, "velFx",     loop_cind, validfind, getc( pdset, "velFx",     loop_cind+1, validfind ))
		setc!( pdset, "velFy",     loop_cind, validfind, getc( pdset, "velFy",     loop_cind+1, validfind ))
		setc!( pdset, "lane",      loop_cind, validfind, getc( pdset, "lane",      loop_cind+1, validfind ))
		setc!( pdset, "nlr",       loop_cind, validfind, getc( pdset, "nlr",       loop_cind+1, validfind ))
		setc!( pdset, "nll",       loop_cind, validfind, getc( pdset, "nll",       loop_cind+1, validfind ))
		setc!( pdset, "curvature", loop_cind, validfind, getc( pdset, "curvature", loop_cind+1, validfind ))
		setc!( pdset, "d_cl",      loop_cind, validfind, getc( pdset, "d_cl",      loop_cind+1, validfind ))
		setc!( pdset, "id",        loop_cind, validfind, getc( pdset, "id",        loop_cind+1, validfind ))
		setc!( pdset, "t_inview",  loop_cind, validfind, getc( pdset, "t_inview",  loop_cind+1, validfind ))

		# account for downshift in mat_other_indmap
		id = getc( pdset, "id", loop_cind, validfind )
		if !isa(id, NAtype)
			mat_ind = pdset.dict_other_idmap[id]
			pdset.mat_other_indmap[validfind, mat_ind] = loop_cind
		end

		loop_cind += 1
	end

	if loop_cind == pdset.maxcarind
		# delete the last entry
		setc!( pdset, "id",        loop_cind, validfind, NA )
		setc!( pdset, "posFx",     loop_cind, validfind, NA )
		setc!( pdset, "posFy",     loop_cind, validfind, NA )
		setc!( pdset, "posFyaw",   loop_cind, validfind, NA )
		setc!( pdset, "velFx",     loop_cind, validfind, NA )
		setc!( pdset, "velFy",     loop_cind, validfind, NA )
		setc!( pdset, "lane",      loop_cind, validfind, NA )
		setc!( pdset, "nlr",       loop_cind, validfind, NA )
		setc!( pdset, "nll",       loop_cind, validfind, NA )
		setc!( pdset, "curvature", loop_cind, validfind, NA )
		setc!( pdset, "d_cl",      loop_cind, validfind, NA )
		setc!( pdset, "id",        loop_cind, validfind, NA )
		setc!( pdset, "t_inview",  loop_cind, validfind, NA )
	end
end
function remove_car_from_frame!( trajdata::DataFrame, carid::Integer, frameind::Integer )

	# returns true if successfully removed

	# remove the car from the specified frame
	@assert( haskey(trajdata, symbol(@sprintf("has_%d", carid))) )
	@assert( frameind > 0 && frameind <= size(trajdata,1))

	# ind = trajdata[symbol(@sprintf("has_%d", carid))][frameind]
	carind = carid2ind(trajdata, carid, frameind)
	if carind == -1
		return false
	end

	while haskey(trajdata, symbol(@sprintf("id_%d", carind+1)))
		# move it down
		setc!( trajdata, "id",            carind, frameind, getc( trajdata, "id",            carind+1, frameind) )
		setc!( trajdata, "posEx",         carind, frameind, getc( trajdata, "posEx",         carind+1, frameind) )
		setc!( trajdata, "posEy",         carind, frameind, getc( trajdata, "posEy",         carind+1, frameind) )
		setc!( trajdata, "velEx",         carind, frameind, getc( trajdata, "velEx",         carind+1, frameind) )
		setc!( trajdata, "velEy",         carind, frameind, getc( trajdata, "velEy",         carind+1, frameind) )
		setc!( trajdata, "posGx",         carind, frameind, getc( trajdata, "posGx",         carind+1, frameind) )
		setc!( trajdata, "posGy",         carind, frameind, getc( trajdata, "posGy",         carind+1, frameind) )
		setc!( trajdata, "velGx",         carind, frameind, getc( trajdata, "velGx",         carind+1, frameind) )
		setc!( trajdata, "velGy",         carind, frameind, getc( trajdata, "velGy",         carind+1, frameind) )
		setc!( trajdata, "yawG",          carind, frameind, getc( trajdata, "yawG",          carind+1, frameind) )
		setc!( trajdata, "lane",          carind, frameind, getc( trajdata, "lane",          carind+1, frameind) )
		setc!( trajdata, "lane_offset",   carind, frameind, getc( trajdata, "lane_offset",   carind+1, frameind) )
		setc!( trajdata, "lane_tangent",  carind, frameind, getc( trajdata, "lane_tangent",  carind+1, frameind) )
		setc!( trajdata, "angle_to_lane", carind, frameind, getc( trajdata, "angle_to_lane", carind+1, frameind) )
		setc!( trajdata, "velLs",         carind, frameind, getc( trajdata, "velLs",         carind+1, frameind) )
		setc!( trajdata, "velLd",         carind, frameind, getc( trajdata, "velLd",         carind+1, frameind) )
		setc!( trajdata, "posRLs",        carind, frameind, getc( trajdata, "posRLs",        carind+1, frameind) )
		setc!( trajdata, "posRLd",        carind, frameind, getc( trajdata, "posRLd",        carind+1, frameind) )

		# account from the downshift in has_ field
		if carind_exists( trajdata, carind+1, frameind )
			carid2 = carind2id( trajdata, carind+1, frameind )
			setc_has!( trajdata, carid2, frameind, carind )
		end

		carind += 1
	end

	# delete the last entry
	setc!( trajdata, "id",            carind, frameind, NA )
	setc!( trajdata, "posEx",         carind, frameind, NA )
	setc!( trajdata, "posEy",         carind, frameind, NA )
	setc!( trajdata, "velEx",         carind, frameind, NA )
	setc!( trajdata, "velEy",         carind, frameind, NA )
	setc!( trajdata, "posGx",         carind, frameind, NA )
	setc!( trajdata, "posGy",         carind, frameind, NA )
	setc!( trajdata, "velGx",         carind, frameind, NA )
	setc!( trajdata, "velGy",         carind, frameind, NA )
	setc!( trajdata, "yawG",          carind, frameind, NA )
	setc!( trajdata, "lane",          carind, frameind, NA )
	setc!( trajdata, "lane_offset",   carind, frameind, NA )
	setc!( trajdata, "lane_tangent",  carind, frameind, NA )
	setc!( trajdata, "angle_to_lane", carind, frameind, NA )
	setc!( trajdata, "velLs",         carind, frameind, NA )
	setc!( trajdata, "velLd",         carind, frameind, NA )
	setc!( trajdata, "posRLs",        carind, frameind, NA )
	setc!( trajdata, "posRLd",        carind, frameind, NA )

	# pull the car from has_ field
	setc_has!( trajdata, carid, frameind, -1 )

	return true
end
function remove_cars_from_frame!( trajdata::DataFrame, frameind::Integer )

	maxncars = get_max_num_cars(trajdata)

	for carind = 0 : maxncars

		# account from the removal in has_ field
		carid = carind2id(trajdata, carind, frameind)
		if carind_exists( trajdata, carind, frameind )

			setc_has!( trajdata, carid, frameind, -1)

			setc!( trajdata, "id",            carind, frameind, NA )
			setc!( trajdata, "posEx",         carind, frameind, NA )
			setc!( trajdata, "posEy",         carind, frameind, NA )
			setc!( trajdata, "velEx",         carind, frameind, NA )
			setc!( trajdata, "velEy",         carind, frameind, NA )
			setc!( trajdata, "posGx",         carind, frameind, NA )
			setc!( trajdata, "posGy",         carind, frameind, NA )
			setc!( trajdata, "velGx",         carind, frameind, NA )
			setc!( trajdata, "velGy",         carind, frameind, NA )
			setc!( trajdata, "yawG",          carind, frameind, NA )
			setc!( trajdata, "lane",          carind, frameind, NA )
			setc!( trajdata, "lane_offset",   carind, frameind, NA )
			setc!( trajdata, "lane_tangent",  carind, frameind, NA )
			setc!( trajdata, "angle_to_lane", carind, frameind, NA )
			setc!( trajdata, "velLs",         carind, frameind, NA )
			setc!( trajdata, "velLd",         carind, frameind, NA )
			setc!( trajdata, "posRLs",        carind, frameind, NA )
			setc!( trajdata, "posRLd",        carind, frameind, NA )
		end
	end
end
function remove_cars_from_frame!( pdset::PrimaryDataset, validfind::Integer )

	for carind = 0 : pdset.maxcarind

		id = getc( pdset, "id", loop_cind, validfind )
		if !isa(id, NAtype)

			mat_ind = pdset.dict_other_idmap[id]
			pdset.mat_other_indmap[validfind, mat_ind] = -1

			setc!( pdset, "id",        carind, validfind, NA )
			setc!( pdset, "posFx",     carind, validfind, NA )
			setc!( pdset, "posFy",     carind, validfind, NA )
			setc!( pdset, "posFyaw",   carind, validfind, NA )
			setc!( pdset, "velFx",     carind, validfind, NA )
			setc!( pdset, "velFy",     carind, validfind, NA )
			setc!( pdset, "lane",      carind, validfind, NA )
			setc!( pdset, "nlr",       carind, validfind, NA )
			setc!( pdset, "nll",       carind, validfind, NA )
			setc!( pdset, "curvature", carind, validfind, NA )
			setc!( pdset, "d_cl",      carind, validfind, NA )
			setc!( pdset, "id",        carind, validfind, NA )
			setc!( pdset, "t_inview",  carind, validfind, NA )
		end
	end
end
function remove_car!( pdset::PrimaryDataset, carid::Integer )

	@assert(haskey(pdset.dict_other_idmap, carid))
	matind = pdset.dict_other_idmap[carid]

	# removes all instances of the given car from the pdset
	for validfind = 1 : length(pdset.validfind2frameind)
		carind = pdset.mat_other_indmap[validfind, matind]
		if carind != -1
			remove_car_from_frame!( pdset, carind, validfind )
		end
	end

	# remove id from mat_other_indmap 
	pdset.mat_other_indmap = pdset.mat_other_indmap[:,[1:matind-1,matind+1:end]]

	# shift everyone else down
	keyset = collect(keys(pdset.dict_other_idmap))
	for key in keyset
		matind2 = pdset.dict_other_idmap[key]
		if matind2 > matind
			pdset.dict_other_idmap[key] = matind2-1
		end
	end

	# remove id from dict_other_idmap
	delete!(pdset.dict_other_idmap, carid)

	pdset
end

function find_slot_for_car!( trajdata::DataFrame, frameind::Integer, maxncars = -2 )

	# returns the index within the frame to which the car can be added
	# adds a new slot if necessary

	if maxncars < -1
		maxncars = get_max_num_cars(trajdata)
	end

	ncars = ncars_in_frame( trajdata, frameind, maxncars)
	if ncars > maxncars
		# we have no extra slots, make one
		return add_car_slot!( trajdata )
	else
		# we have enough slots; use the first one
		return ncars
	end
end
function find_slot_for_car!( pdset::DataFrame, validfind::Integer )

	# returns the index within the frame to which the car can be added
	# adds a new slot if necessary

	ncars = ncars_in_frame( pdset, validfind )
	if ncars > pdset.maxcarind
		# we have no extra slots, make one
		return add_car_slot!( pdset )
	else
		# we have enough slots; use the first one
		return ncars # this is correct since ncars-1 is the last ind and we want +1 from that
	end
end


function ncars_in_frame( trajdata::DataFrame, frameind::Integer, maxncars = -2 )

	if maxncars < -1
		maxncars = get_max_num_cars(trajdata)
	end

	count = 0
	for carind = 0 : maxncars
		if carind_exists(trajdata, carind, frameind)
			count += 1
		end
	end
	count
end
function ncars_in_frame( pdset::PrimaryDataset, validfind::Integer )

	# returns the number of cars in a given frame, not including the ego car
	count = 0
	for carind = 0 : pdset.maxcarind
		if indinframe(pdset, carind, validfind)
			count += 1
		end
	end
	count
end

function add_car_slot!( trajdata::DataFrame, maxncars = -2 )
	# increases the number of observed car column sets by 1

	if maxncars < -1
		maxncars = get_max_num_cars(trajdata)
	end

	carind = maxncars+1

	na_arr = Array(Any, size(trajdata,1))
	fill!(na_arr, NA)

	trajdata[symbol(@sprintf("id_%d",            carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("posEx_%d",         carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("posEy_%d",         carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("velEx_%d",         carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("velEy_%d",         carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("posGx_%d",         carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("posGy_%d",         carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("velGx_%d",         carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("velGy_%d",         carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("yawG_%d",          carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("lane_%d",          carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("lane_offset_%d",   carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("lane_tangent_%d",  carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("angle_to_lane_%d", carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("velLs_%d",         carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("velLd_%d",         carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("posRLs_%d",        carind))] = copy(na_arr)
	trajdata[symbol(@sprintf("posRLd_%d",        carind))] = copy(na_arr)

	return maxncars+1
end
function add_car_slot!( pdset::PrimaryDataset )
	# increases the number of observed car column sets by 1
	# returns the index of the newly added slot

	carind = pdset.maxcarind + 1 # new index to add

	N = size(pdset.df_other_primary,1)

	float64_arr = DataArray(Float64, N)
	int8_arr = DataArray(Int8, N)
	uint32_arr = DataArray()

	trajdata[symbol(@sprintf("posFx_%d",     carind))] = DataArray(Float64, N)
	trajdata[symbol(@sprintf("posFy_%d",     carind))] = DataArray(Float64, N)
	trajdata[symbol(@sprintf("posFyaw_%d",   carind))] = DataArray(Float64, N)
	trajdata[symbol(@sprintf("velFx_%d",     carind))] = DataArray(Float64, N)
	trajdata[symbol(@sprintf("velFy_%d",     carind))] = DataArray(Float64, N)
	trajdata[symbol(@sprintf("lane_%d",      carind))] = DataArray(Int8,    N)
	trajdata[symbol(@sprintf("nlr_%d",       carind))] = DataArray(Int8,    N)
	trajdata[symbol(@sprintf("nll_%d",       carind))] = DataArray(Int8,    N)
	trajdata[symbol(@sprintf("curvature_%d", carind))] = DataArray(Float64, N)
	trajdata[symbol(@sprintf("d_cl_%d",      carind))] = DataArray(Float64, N)
	trajdata[symbol(@sprintf("id_%d",        carind))] = DataArray(Uint32,  N)
	trajdata[symbol(@sprintf("t_inview_%d",  carind))] = DataArray(Float64, N)

	return carind
end

function frames_contain_carid( trajdata::DataFrame, carid::Integer, startframeind::Int, horizon::Int; frameskip::Int=1 )
		
	if startframeind+horizon > size(trajdata,1)
		return false
	end

	symb = symbol(@sprintf("has_%d", carid))
	if haskey(trajdata, symb)

		for j = 0 : frameskip : horizon
			if trajdata[symb][startframeind+j] == -1
				return false
			end
		end
	end
	return true
end

function extract_continuous_tracks( trajdata::DataFrame )

	segdict = Dict{Integer,Array{Integer,2}}() # carid -> segment array
	for key in names(trajdata)

		if ismatch(r"^has_(\d)+$", string(key))

			carid = int(match(r"(\d)+", string(key)).match)

			arr = array(trajdata[symbol(@sprintf("has_%d", carid))],0)

			for frameind = 1 : length(arr)
				if arr[frameind] != -1
					if !haskey(segdict, carid)
						# start a new segment
						segdict[carid] = [frameind frameind]'
					elseif segdict[carid][2,end] == frameind-1
						# extend the segment further by one
						segdict[carid][2,end] = frameind
					else
						# start a new segment
						segdict[carid] = hcat(segdict[carid], [frameind frameind]')
					end
				end
			end
		end
	end

	return segdict
end

function spanning_validfinds( pdset::PrimaryDataset, carid::Integer )
	# obtain the min and max validfind where this car can be found
	mat_ind = pdset.dict_other_idmap[carid]
	ind_lo = findfirst(x->x!=-1, pdset.mat_other_indmap[:,mat_ind])
	ind_hi = findfirst(x->x!=-1, reverse(pdset.mat_other_indmap[:,mat_ind]))
	(ind_lo, ind_hi)
end

# Helper Functions
# -------------------------
function quat2euler{T <: Real}( quat::Array{T,1} )
	# convert a quaternion to roll-pitch-yaw
	
	den = norm(quat)
	w = quat[1]/den
	x = quat[2]/den
	y = quat[3]/den
	z = quat[4]/den

	roll  = atan2(y*z+w*x, 0.5-x^2-y^2)
	pitch = asin(-2*(x*z + w*y))
	yaw   = atan2(x*y+w*z, 0.5-y^2-z^2)

	return (roll, pitch, yaw)
end
function euler2quat( roll::Real, pitch::Real, yaw::Real )

	cr = cos(roll/2.0)
	sr = sin(roll/2.0)
	cp = cos(pitch/2.0)
	sp = sin(pitch/2.0)
	cy = cos(yaw/2.0)
	sy = sin(yaw/2.0)

	w2 = cr*cp*cy + sr*sp*sy
	x2 = sr*cp*cy - cr*sp*sy
	y2 = cr*sp*cy + sr*cp*sy
	z2 = cr*cp*sy - sr*sp*cy

	return [w2, x2, y2, z2]
end
function global2ego( egocarGx::Real, egocarGy::Real, egocarYaw::Real, posGx::Real, posGy::Real )

	pt = [posGx, posGy]
	eo = [egocarGx, egocarGy]

	# translate to be relative to ego car
	pt -= eo

	# rotate to be in ego car frame
	R = [ cos(egocarYaw) sin(egocarYaw);
	     -sin(egocarYaw) cos(egocarYaw)]
	pt = R*pt

	return (pt[1], pt[2]) # posEx, posEy
end
function global2ego( egocarGx::Real, egocarGy::Real, egocarYaw::Real, posGx::Real, posGy::Real, yawG::Real )

	posEx, posEy = global2ego( egocarGx, egocarGy, egocarYaw, posGx, posGy)
	yawE = mod2pi(yawG - egocarYaw)

	return (posEx, posEy, yawE)
end
function ego2global( egocarGx::Real, egocarGy::Real, egocarYaw::Real, posEx::Real, posEy::Real )

	pt = [posEx, posEy]
	eo = [egocarGx, egocarGy]

	# rotate
	R = [ cos(egocarYaw) -sin(egocarYaw);
	      sin(egocarYaw)  cos(egocarYaw)]
	pt = R*pt

	# translate
	pt += eo

	return (pt[1], pt[2]) # posGx, posGy
end
function ego2global( egocarGx::Real, egocarGy::Real, egocarYaw::Real, posEx::Real, posEy::Real, yawE::Real )

	posGx, posGy = ego2global( egocarGx, egocarGy, egocarYaw, posEx, posEy)
	yawG = mod2pi(yawG + egocarYaw)

	return (posGx, posGy, yawG)
end
isnum(s::String) = ismatch(NUMBER_REGEX, s) || s == "Inf" || s == "inf" || s == "-Inf" || s == "-inf" || s == "+Inf" || s == "+inf"

end # end module
