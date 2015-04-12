##############################################################################
##
##  Features
##
##  Allows for the extraction of feature values for cars in scenes
##
##############################################################################

module Features

import Base: symbol, get
import Base.Meta: quot

using Trajdata
using DataArrays
using LaTeXStrings
using StreetMap
using Curves

# ----------------------------------
# exports

export AbstractFeature, InherentFeature, MarkovFeature, IteratedFeature
export RichVehicle
export FRAME_PER_SEC, SEC_PER_FRAME
export YAW, POSFX, POSFY, SPEED, VELFX, VELFY, DELTA_SPEED_LIMIT, TURNRATE, TURNRATE_GLOBAL, ACC, ACCFX, ACCFY, ISEGO
export CL, D_CL, D_ML, D_MR, SCENEVELFX, D_ONRAMP, D_OFFRAMP, TIMETOCROSSING_RIGHT, TIMETOCROSSING_LEFT
export D_MERGE, D_SPLIT, ID
export A_REQ_STAYINLANE, N_LANE_R, N_LANE_L, HAS_LANE_R, HAS_LANE_L, LANECURVATURE
export INDFRONT, HAS_FRONT, D_X_FRONT, D_Y_FRONT, V_X_FRONT, V_Y_FRONT, YAW_FRONT, TURNRATE_FRONT, A_REQ_FRONT, TTC_X_FRONT, TIMEGAP_X_FRONT
export INDREAR,  HAS_REAR,  D_X_REAR,  D_Y_REAR,  V_X_REAR,  V_Y_REAR,  YAW_REAR,  TURNRATE_REAR,  A_REQ_REAR,  TTC_X_REAR,  TIMEGAP_X_REAR
export INDLEFT,             D_X_LEFT,  D_Y_LEFT,  V_X_LEFT,  V_Y_LEFT,  YAW_LEFT,  TURNRATE_LEFT,  A_REQ_LEFT,  TTC_X_LEFT,  TIMEGAP_X_LEFT
export INDRIGHT,            D_X_RIGHT, D_Y_RIGHT, V_X_RIGHT, V_Y_RIGHT, YAW_RIGHT, TURNRATE_RIGHT, A_REQ_RIGHT, TTC_X_RIGHT, TIMEGAP_X_RIGHT
export GAINING_ON_FRONT

export FUTURETURNRATE_250MS, FUTUREACCELERATION_250MS, FUTUREDESIREDANGLE_250MS, FUTUREDESIREDSPEED_250MS, FUTUREACCELCONTROL_250MS
export FUTURETURNRATE_500MS, FUTUREACCELERATION_500MS, FUTUREDESIREDANGLE_500MS, FUTUREDESIREDSPEED_500MS, FUTUREACCELCONTROL_500MS
export FUTUREDELTAY2S, FUTUREDELTAY1S, FUTUREDELTAY_250MS
export TIMETOLANECROSSING, TIMESINCELANECROSSING, LANECHANGEDIR, LANECHANGE2S, LANECHANGE1S, LANECHANGE500MS

export TIME_CONSECUTIVE_BRAKE, TIME_CONSECUTIVE_ACCEL
export      PASTACC250MS,      PASTACC500MS,      PASTACC750MS,      PASTACC1S
export PASTTURNRATE250MS, PASTTURNRATE500MS, PASTTURNRATE750MS, PASTTURNRATE1S
export    PASTVELFY250MS,    PASTVELFY500MS,    PASTVELFY750MS,    PASTVELFY1S
export     PASTD_CL250MS,     PASTD_CL500MS,     PASTD_CL750MS,     PASTD_CL1S

export     MAXACCFX250MS,     MAXACCFX500MS,     MAXACCFX750MS,     MAXACCFX1S,     MAXACCFX1500MS,     MAXACCFX2S,     MAXACCFX2500MS,     MAXACCFX3S,     MAXACCFX4S
export     MAXACCFY250MS,     MAXACCFY500MS,     MAXACCFY750MS,     MAXACCFY1S,     MAXACCFY1500MS,     MAXACCFY2S,     MAXACCFY2500MS,     MAXACCFY3S,     MAXACCFY4S
export  MAXTURNRATE250MS,  MAXTURNRATE500MS,  MAXTURNRATE750MS,  MAXTURNRATE1S,  MAXTURNRATE1500MS,  MAXTURNRATE2S,  MAXTURNRATE2500MS,  MAXTURNRATE3S,  MAXTURNRATE4S
export    MEANACCFX250MS,    MEANACCFX500MS,    MEANACCFX750MS,    MEANACCFX1S,    MEANACCFX1500MS,    MEANACCFX2S,    MEANACCFX2500MS,    MEANACCFX3S,    MEANACCFX4S
export    MEANACCFY250MS,    MEANACCFY500MS,    MEANACCFY750MS,    MEANACCFY1S,    MEANACCFY1500MS,    MEANACCFY2S,    MEANACCFY2500MS,    MEANACCFY3S,    MEANACCFY4S
export MEANTURNRATE250MS, MEANTURNRATE500MS, MEANTURNRATE750MS, MEANTURNRATE1S, MEANTURNRATE1500MS, MEANTURNRATE2S, MEANTURNRATE2500MS, MEANTURNRATE3S, MEANTURNRATE4S
export     STDACCFX250MS,     STDACCFX500MS,     STDACCFX750MS,     STDACCFX1S,     STDACCFX1500MS,     STDACCFX2S,     STDACCFX2500MS,     STDACCFX3S,     STDACCFX4S
export     STDACCFY250MS,     STDACCFY500MS,     STDACCFY750MS,     STDACCFY1S,     STDACCFY1500MS,     STDACCFY2S,     STDACCFY2500MS,     STDACCFY3S,     STDACCFY4S
export  STDTURNRATE250MS,  STDTURNRATE500MS,  STDTURNRATE750MS,  STDTURNRATE1S,  STDTURNRATE1500MS,  STDTURNRATE2S,  STDTURNRATE2500MS,  STDTURNRATE3S,  STDTURNRATE4S

export SUBSET_EMERGENCY, SUBSET_FREE_FLOW, SUBSET_CAR_FOLLOWING, SUBSET_LANE_CROSSING, SUBSET_SUSTAINED_CROSSING
export SUBSET_AT_SIXTYFIVE, SUBSET_AUTO

export description, units, isint, isbool, lowerbound, upperbound, symbol, lsymbol, couldna, get, symbol2feature
export allfeatures
export NA_ALIAS

# ----------------------------------
# constants

const FRAME_PER_SEC  = 20                     # source frames per second
const SEC_PER_FRAME  = 1.0/FRAME_PER_SEC
const NA_ALIAS = Inf                          # value assigned to floating point values when not available
const THRESHOLD_DY_CONSIDERED_IN_FRONT =  2.0 # number of meters of horizontal separation which allows two cars to be considered in front of / behind each other
const THRESHOLD_A_REQ                  = 10.0 # maximum required absolute value of required acceleration (m/s^2)
const THRESHOLD_TIME_TO_CROSSING       = 10.0 # maximum time to crossing
const THRESHOLD_TIME_TO_COLLISION      = 10.0 # maximum time to collision
const THRESHOLD_TIMEGAP                = 10.0 # maximum timegap
const THRESHOLD_TIMETOLANECROSSING     = 10.0 # [s]
const THRESHOLD_TIMESINCELANECROSSING  = 10.0 # [s]
const THRESHOLD_TIMECONSECUTIVEBRAKE   = 10.0 # [s]
const THRESHOLD_TIMECONSECUTIVEACCEL   = 10.0 # [s]
const SPEED_LIMIT                      = 29.06 # [m/s]
const KP_DESIRED_ANGLE                 = 1.0 # unitary feedback is the general case
const KP_DESIRED_SPEED                 = 0.2 # 

# ----------------------------------
# types

abstract AbstractFeature
abstract InherentFeature      <: AbstractFeature # are always directly accessible from the scene data (posFx, vel) (not stored in meta)
abstract MarkovFeature        <: AbstractFeature # only depend on the current state, but must be extracted (d_x_front, etc.)
abstract IteratedFeature      <: AbstractFeature # require aggregated knowledge over time. They call an update function
abstract UnextractableFeature <: AbstractFeature # features that do not support get()

type BoundingBox
	# a rectangle with length and width
	l :: Float64
	w :: Float64
end
type PointSE2
	# a point in the 2D special euclidean group
	# represents a position and rotation
	x   :: Float64
	y   :: Float64
	yaw :: Float64 # [rad]
end

const CAR_WIDTH_STANDARD  = 2.0 # [m]
const CAR_LENGTH_STANDARD = 4.6 # [m]
const LANE_WIDTH_STANDARD = 3.5 # [m]

abstract AbstractVehicle
type RichVehicle <: AbstractVehicle
	pos  :: PointSE2    # location in the frenet frame (s,d) [m, rad]
	spd  :: Real        # speed [m/s]
	box  :: BoundingBox # the collision bounding box [m]
	RichVehicle(x::Real, y::Real, θ::Real, v::Real) = new(PointSE2(x,y,θ), v, BoundingBox(CAR_WIDTH_STANDARD,CAR_LENGTH_STANDARD), Dict{Symbol,Float64}())
end

type TrafficScene
	road      :: Road            # road definition
	cars      :: Vector{AbstractVehicle}
	dist_fore :: Float64 # distance ahead of the ego car in which other agents are considered [m]
	dist_rear :: Float64 # distance behind the ego car in which other agents are considered [m]
end

# ----------------------------------
# globals
sym2ftr = Dict{Symbol,AbstractFeature}()

meta = Dict{Int,Dict{Symbol,Float64}}() # carid -> meta map
meta_validfind = 0                      # the validfind corresponding to the current meta datascene

# ----------------------------------
# abstract feature functions

symbol2feature(sym::Symbol)    = sym2ftr[sym]
description(::AbstractFeature) = "None"
units(      ::AbstractFeature) = "None"
isint(      ::AbstractFeature) = error("Not implemented!") # true if the feature takes on only integer values
isbool(     ::AbstractFeature) = error("Not implemented!") # true if the feature takes on only boolean values (0,1) (will also have isint = true)
upperbound( ::AbstractFeature) = error("Not implemented!") # max value of the feature
lowerbound( ::AbstractFeature) = error("Not implemented!") # min value of the feature
couldna(    ::AbstractFeature) = true                      # whether the given variable could be NA (naalias is Inf)
symbol(     ::AbstractFeature) = error("Not implemented!")
lsymbol(    ::AbstractFeature) = error("Not implemented!")
get(        ::AbstractFeature, ::TrafficScene, carind::Int) = error("Not implemented!")
get(        ::AbstractFeature, ::PrimaryDataset, carind::Int, validfind::Int) = error("Not implemented!")

allfeatures() = collect(values(sym2ftr)) # array of all features

 get(F::InherentFeature, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int) =  get(F, pdset, carind, validfind)
_get(F::AbstractFeature, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int) = _get(F, pdset, carind, validfind)

function get(F::MarkovFeature, scene::TrafficScene, carind::Int)
	car = scene.cars[carind]
	sym = symbol(F)
	if haskey(car.meta, F)
		return car.meta[F]
	end
	value = _get(F, scene, carind)
	car.meta[sym] = value
	@assert(!isa(value, NAtype))
	value
end
function get(F::MarkovFeature, pdset::PrimaryDataset, carind::Int, validfind::Int)
	# a wrapper for get that automatically does the meta[] stuff for you
	sym = symbol(F)
	if checkmeta(carind, validfind, sym)
		return meta[carind][sym]
	end
	value = _get(F, pdset, carind, validfind)
	meta[carind][sym] = value
	@assert(!isa(value, NAtype))
	value
end
function get(F::UnextractableFeature, pdset::PrimaryDataset, carind::Int, validfind::Int)
	# a wrapper for get that automatically does the meta[] stuff for you
	sym = symbol(F)
	if checkmeta(carind, validfind, sym)
		return meta[carind][sym]
	end
	value = _get(F, pdset, carind, validfind)
	meta[carind][sym] = value
	@assert(!isa(value, NAtype))
	value
end
function get(F::MarkovFeature, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	# a wrapper for get that automatically does the meta[] stuff for you
	sym = symbol(F)
	if checkmeta(carind, validfind, sym)
		return meta[carind][sym]
	end
	value = _get(F, pdset, sn, carind, validfind)
	meta[carind][sym] = value
	@assert(!isa(value, NAtype))
	value
end
function get(F::UnextractableFeature, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	# a wrapper for get that automatically does the meta[] stuff for you
	sym = symbol(F)
	if checkmeta(carind, validfind, sym)
		return meta[carind][sym]
	end
	value = _get(F, pdset, sn, carind, validfind)
	meta[carind][sym] = value
	@assert(!isa(value, NAtype))
	value
end


# ----------------------------------

function create_feature_basics{F<:AbstractFeature}( 
	name         :: String,
	feature_type :: Type{F},
	unit         :: String,
	isint        :: Bool,
	isbool       :: Bool,
	ub           :: Float64,
	lb           :: Float64,
	could_be_na  :: Bool,
	sym          :: Symbol,
	lstr         :: LaTeXString,
	desc         :: String
	)

	for feature in values(sym2ftr)
		@assert(desc != description(feature), "desc: $name -> $feature")
		@assert(sym  != symbol(feature), "symb: $name -> $feature")
		@assert(lstr != lsymbol(feature), "lstr: $name -> $feature")
	end
	@assert(ub >= lb)

	feature_name = symbol("Feature_" * name)
	const_name   = symbol(uppercase(name))
	sym_feature  = quot(sym)

	@eval begin
		immutable $feature_name <: $feature_type end
		const       $const_name  = ($feature_name)()
		description( ::$feature_name)  = $desc
		units(       ::$feature_name)  = $unit
		isint(       ::$feature_name)  = $isint
		isbool(      ::$feature_name)  = $isbool
		upperbound(  ::$feature_name)  = $ub
		lowerbound(  ::$feature_name)  = $lb
		couldna(     ::$feature_name)  = $could_be_na
		symbol(      ::$feature_name)  = $sym_feature
		lsymbol(     ::$feature_name)  = $lstr
		sym2ftr[symbol(  $const_name)] = $const_name
	end
end

# ----------------------------------
# inherent features

create_feature_basics("Yaw", InherentFeature, "rad", false, false, Inf, -Inf, false, :yaw, L"\psi", "angle relative to the closest lane")
get(         ::Feature_Yaw, scene::TrafficScene, carind::Int) = scene.cars[carind].pos.yaw
function get(::Feature_Yaw, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		@assert(frameind != 0)
		return gete(pdset, :posFyaw, frameind)
	end
	getc(pdset, "posFyaw", carind, validfind)
end

create_feature_basics("PosFx", InherentFeature, "m", false, false, Inf, -Inf, false, :posFx, L"p^F_x", "x position in the frenet frame")
get(         ::Feature_PosFx, scene::TrafficScene, carind::Int) = scene.cars[carind].pos.x
function get(::Feature_PosFx, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		return gete(pdset, :posFx, frameind)
	end
	getc(pdset, "posFx", carind, validfind)
end

create_feature_basics("PosFy", InherentFeature, "m", false, false, Inf, -Inf, false, :posFy, L"p^F_y", "y position in the frenet frame")
get(         ::Feature_PosFy, scene::TrafficScene, carind::Int) = scene.cars[carind].pos.y
function get(::Feature_PosFy, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		return gete(pdset, :posFy, frameind)
	end
	getc(pdset, "posFy", carind, validfind)
end

create_feature_basics("Speed", InherentFeature, "m/s", false, false, Inf, -Inf, false, :speed, L"\|v\|", "speed")
get(         ::Feature_Speed, scene::TrafficScene, carind::Int) = scene.cars[carind].spd
function get(::Feature_Speed, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		velx = gete(pdset, :velFx, frameind)
		vely = gete(pdset, :velFy, frameind)
		return hypot(velx, vely)
	end
	velx = getc(pdset, "velFx", carind, validfind)
	vely = getc(pdset, "velFy", carind, validfind)
	return hypot(velx, vely)
end

create_feature_basics("Delta_Speed_Limit", InherentFeature, "m/s", false, false, Inf, -Inf, false, :delta_speed_limit, L"Δv_{\text{limit}}", "difference between current speed and speed limit")
get(         ::Feature_Delta_Speed_Limit, scene::TrafficScene, carind::Int) = SPEED_LIMIT - scene.cars[carind].spd
function get(::Feature_Delta_Speed_Limit, pdset::PrimaryDataset, carind::Int, validfind::Int)
	speed = get(SPEED, pdset, carind, validfind)
	SPEED_LIMIT - speed
end

create_feature_basics("VelFx", InherentFeature, "m/s", false, false, Inf, -Inf, false, :velFx, L"v^F_x", "velocity along the lane centerline")
get(         ::Feature_VelFx, scene::TrafficScene, carind::Int) = scene.cars[carind].spd * cos(scene.cars[carind].pos.yaw)
function get(::Feature_VelFx, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		return gete(pdset, :velFx, frameind)
	end
	getc(pdset, "velFx", carind, validfind)
end

create_feature_basics("VelFy", InherentFeature, "m/s", false, false, Inf, -Inf, false, :velFy, L"v^F_y", "velocity perpendicular to the lane centerline")
get(         ::Feature_VelFy, scene::TrafficScene, carind::Int) = scene.cars[carind].spd * sin(scene.cars[carind].pos.yaw)
function get(::Feature_VelFy, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		return gete(pdset, :velFy, frameind)
	end
	getc(pdset, "velFy", carind, validfind)
end

create_feature_basics("IsEgo", InherentFeature, "-", true, true, 1.0, 0.0, false, :isego, L"1_{ego}", "whether the car is the ego car")
get(         ::Feature_IsEgo, scene::TrafficScene, carind::Int) = carind == CARIND_EGO
get(         ::Feature_IsEgo, pdset::PrimaryDataset, carind::Int, validfind::Int) = carind == CARIND_EGO

# ----------------------------------
# markov features

create_feature_basics( "CL", MarkovFeature, "-", true, false, Inf, 1.0, false, :cl, L"cl", "index of the closest centerline. 1 is rightmost")
function _get(::Feature_CL, scene::TrafficScene, carind::Int)
	car = scene.cars[carind]
	# generate the lane centerlines
	posFy = car.pos.y
	n_lanes = scene.road.n_lanes
	lane_centers = [0:(n_lanes-1)].*scene.road.lane_width # [0,lw,2lw,...]
	indmin(abs(lane_centers .- posFy))
end
function _get(::Feature_CL, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		return gete(pdset, :lane, frameind)
	end
	getc(pdset, "lane", carind, validfind)
end
		  get(::Feature_CL, ::PrimaryDataset, ::StreetNetwork, ::Int, ::Int) = error("CL not implemented with StreetNetworks")

create_feature_basics( "D_CL", MarkovFeature, "m", false, false, Inf, -Inf, false, :d_cl, L"d_{cl}", "distance to the closest centerline")
function _get(::Feature_D_CL, scene::TrafficScene, carind::Int)
	car = scene.cars[carind]
	lane = get(CL, scene, carind)
	car.pos.y - (lane-1)*scene.road.lane_width
end
function  get(::Feature_D_CL, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		return gete(pdset, :d_cl, frameind)
	end
	getc(pdset, "d_cl", carind, validfind)
end
		  get(::Feature_D_CL, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int) = get(D_CL, pdset, carind, validfind)

create_feature_basics( "D_ML", MarkovFeature, "m", false, false, Inf, 0.0, true, :d_ml, L"d_{ml}", "lateral distance between center of car and the left lane marker")
function _get(::Feature_D_ML, scene::TrafficScene, carind::Int)
 	d_cl = get(D_CL, scene, carind)
 	d_cl > 0.5LANE_WIDTH_STANDARD ? NA_ALIAS : 0.5LANE_WIDTH_STANDARD - d_cl
end
function  get(::Feature_D_ML, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		d_ml = gete(pdset, :d_ml, frameind)
		return isa(d_ml, NAtype) ? NA_ALIAS : d_ml
	end
	d_ml = getc(pdset, "d_ml", carind, validfind)
	return isa(d_ml, NAtype) ? NA_ALIAS : d_ml
end
	      get(::Feature_D_ML, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int) = get(D_ML, pdset, carind, validfind)

create_feature_basics( "D_MR", MarkovFeature, "m", false, false, Inf, 0.0, true, :d_mr, L"d_{mr}", "lateral distance (strictly positive) between center of car and the right lane marker")
function _get(::Feature_D_MR, scene::TrafficScene, carind::Int)
 	d_cl = get(D_CL, scene, carind)
 	d_cl < 0.5LANE_WIDTH_STANDARD ? NA_ALIAS : 0.5LANE_WIDTH_STANDARD + d_cl
end
function  get(::Feature_D_MR, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		d_mr = gete(pdset, :d_mr, frameind)
		return isa(d_mr, NAtype) ? NA_ALIAS : d_mr
	end
	d_mr = getc(pdset, "d_mr", carind, validfind)
	return isa(d_mr, NAtype) ? NA_ALIAS : d_mr
end
	      get(::Feature_D_MR, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int) = get(D_MR, pdset, carind, validfind)

create_feature_basics( "TimeToCrossing_Right", MarkovFeature, "s", false, false, Inf, 0.0, true, :ttcr_mr, L"ttcr^{mr}_y", "time to cross the right marking of assigned lane")
function  get(::Feature_TimeToCrossing_Right, scene::TrafficScene, carind::Int)
	d_mr = get(D_MR, scene, carind)
	velFy = get(VELFY, scene, carind)
	d_mr > 0.0 && velFy < 0.0 ? min(d_mr / velFy,THRESHOLD_TIME_TO_CROSSING) : NA_ALIAS
end
function  get(::Feature_TimeToCrossing_Right, pdset::PrimaryDataset, carind::Int, validfind::Int)
	d_mr = get(D_MR, pdset, carind, validfind)
	velFy = get(VELFY, pdset, carind, validfind)
	d_mr > 0.0 && velFy < 0.0 ? min(-d_mr / velFy,THRESHOLD_TIME_TO_CROSSING) : NA_ALIAS
end
	      get(::Feature_TimeToCrossing_Right, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int) = get(TIMETOCROSSING_RIGHT, pdset, carind, validfind)

create_feature_basics( "TimeToCrossing_Left", MarkovFeature, "s", false, false, Inf, 0.0, true, :ttcr_ml, L"ttcr^{ml}_y", "time to cross the left marking of assigned lane")
function  get(::Feature_TimeToCrossing_Left, scene::TrafficScene, carind::Int)
	d_ml = get(D_MR, scene, carind)
	velFy = get(VELFY, scene, carind)
	d_ml > 0.0 && velFy > 0.0 ? min(d_ml / velFy,THRESHOLD_TIME_TO_CROSSING) : NA_ALIAS
end
function  get(::Feature_TimeToCrossing_Left, pdset::PrimaryDataset, carind::Int, validfind::Int)
	d_ml = get(D_ML, pdset, carind, validfind)
	velFy = get(VELFY, pdset, carind, validfind)
	d_ml > 0.0 && velFy > 0.0 ? min(d_ml / velFy,THRESHOLD_TIME_TO_CROSSING) : NA_ALIAS
end
function  get(::Feature_TimeToCrossing_Left, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	d_ml = get(D_ML, pdset, sn, carind, validfind)
	velFy = get(VELFY, pdset, sn, carind, validfind)
	d_ml < 0.0 && velFy > 0.0 ? min(-d_ml / velFy,THRESHOLD_TIME_TO_CROSSING) : NA_ALIAS
end

create_feature_basics( "A_REQ_StayInLane", MarkovFeature, "m/s2", false, false, Inf, -Inf, true, :a_req_stayinlane, L"a^{req}_y", "acceleration required to stay in the current lane")
function get(::Feature_A_REQ_StayInLane, scene::TrafficScene, carind::Int)
	velFy = get(VELFY, scene, carind)
	d_mr = get(D_MR, scene, carind)
	if d_mr > 0.0
		return -min(0.5velFy*velFy / d_mr, THRESHOLD_A_REQ)
	end
	d_ml = get(D_ML, scene, carind)
	d_ml > 0.0 ?  -min(0.5velFy*velFy / d_ml, THRESHOLD_A_REQ) : NA_ALIAS
end
function get(::Feature_A_REQ_StayInLane, pdset::PrimaryDataset, carind::Int, validfind::Int)
	velFy = get(VELFY, pdset, carind, validfind)
	d_mr = get(D_MR, pdset, carind, validfind)
	if d_mr > 0.0
		return -min(0.5velFy*velFy / d_mr, THRESHOLD_A_REQ)
	end
	d_ml = get(D_ML, pdset, carind, validfind)
	d_ml > 0.0 ?  -min(0.5velFy*velFy / d_ml, THRESHOLD_A_REQ) : NA_ALIAS
end
function get(::Feature_A_REQ_StayInLane, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	velFy = get(VELFY, pdset, sn, carind, validfind)

	if velFy > 0.0
		d_mr = get(D_MR, pdset, sn, carind, validfind)
		return d_mr > 0.0 ? min( 0.5velFy*velFy / d_mr, THRESHOLD_A_REQ) : NA_ALIAS
	else
		d_ml = get(D_ML, pdset, sn, carind, validfind)
		return d_ml < 0.0 ? min(-0.5velFy*velFy / d_ml, THRESHOLD_A_REQ) : NA_ALIAS
	end
end

create_feature_basics( "N_LANE_L", MarkovFeature, "-", true, false, 10.0, 0.0, false, :n_lane_left, L"nlane_l", "number of lanes on the left side of this vehicle")
function get(::Feature_N_LANE_L, scene::TrafficScene, carind::Int)
	n_lanes = scene.road.n_lanes
	lane    = get(CL, scene, carind)
	float64(n_lanes - lane)
end
function get(::Feature_N_LANE_L, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		nll = float64(gete(pdset, :nll, frameind))
		@assert(0.0 ≤ nll ≤ 6.0)
		return nll
	end
	float64(getc(pdset, "nll", carind, validfind))
end
		 get(::Feature_N_LANE_L, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int) = get(N_LANE_L, pdset, carind, validfind)

create_feature_basics( "N_LANE_R", MarkovFeature, "-", true, false, 10.0, 0.0, false, :n_lane_right, L"nlane_r", "number of lanes on the right side of this vehicle")
		 get(::Feature_N_LANE_R, scene::TrafficScene, carind::Int) = get(CL, scene, carind) + 1.0
function get(::Feature_N_LANE_R, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		return float64(gete(pdset, :nlr, frameind))
	end
	float64(getc(pdset, "nlr", carind, validfind))
end
	     get(::Feature_N_LANE_R, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int) = get(N_LANE_R, pdset, carind, validfind)

create_feature_basics( "HAS_LANE_R", MarkovFeature, "-", true, true, 1.0, 0.0, false, :has_lane_right, L"\exists_{\text{lane}}^\text{r}", "whether at least one lane exists to right")
function get(::Feature_HAS_LANE_R, scene::TrafficScene, carind::Int)

	float64(get(N_LANE_R, scene, carind) > 0.0)
end
function get(::Feature_HAS_LANE_R, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	float64(get(N_LANE_R, pdset, carind, validfind) > 0.0)
end
	     get(::Feature_HAS_LANE_R, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int) = get(HAS_LANE_R, pdset, carind, validfind)

create_feature_basics( "HAS_LANE_L", MarkovFeature, "-", true, true, 1.0, 0.0, false, :has_lane_left, L"\exists_{\text{lane}}^\text{l}", "whether at least one lane exists to the left")
function get(::Feature_HAS_LANE_L, scene::TrafficScene, carind::Int)
	
	float64(get(N_LANE_L, scene, carind) > 0.0)
end
function get(::Feature_HAS_LANE_L, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	float64(get(N_LANE_L, pdset, carind, validfind) > 0.0)
end
		 get(::Feature_HAS_LANE_L, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int) = get(HAS_LANE_L, pdset, carind, validfind)


create_feature_basics( "IndFront", MarkovFeature, "-", true, false, Inf, -1.0, true, :ind_front, L"i_{front}", "index of the closest car in front")
function _get(::Feature_IndFront, scene::TrafficScene, carind::Int)
	car = scene.cars[carind]

	mylane = get(CL, scene, carind)
	myFy   = car.pos.y
	myFx   = car.pos.x
	myVy   = car.spd*sin(car.pos.yaw)
	myVx   = car.spd*cos(car.pos.yaw)

	min_dx = Inf
	ind_front = -2

	d_y_front = 0.0
	v_x_front = 0.0
	v_y_front = 0.0
	yaw_front = 0.0

	for (i,car2) in enumerate(scene.cars)
		if i == carind
			continue
		end

		dy = car2.pos.y - myFy
		dlane = get(CL, scene, i)
		if isapprox(dlane, mylane) || abs(dy) < THRESHOLD_DY_CONSIDERED_IN_FRONT
			dx = car2.pos.x - myFx
			if dx < min_dx && dx > 0.0
				min_dx = dx
				d_y_front = car2.pos.y - myFy
				v_x_front = car2.spd*cos(car2.pos.yaw) - myVx
				v_y_front = car2.spd*sin(car2.pos.yaw) - myVy
				yaw_front = car2.pos.yaw
				ind_front = i
			end
		end
	end

	if isinf(min_dx)
		car.meta[:d_x_front] = NA_ALIAS
		car.meta[:d_y_front] = NA_ALIAS
		car.meta[:v_x_front] = NA_ALIAS
		car.meta[:v_y_front] = NA_ALIAS
		car.meta[:yaw_front] = NA_ALIAS
		return NA_ALIAS
	end

	car.meta[:d_x_front] = min_dx
	car.meta[:d_y_front] = d_y_front
	car.meta[:v_x_front] = v_x_front
	car.meta[:v_y_front] = v_y_front
	car.meta[:yaw_front] = yaw_front
	ind_front
end
function _get(::Feature_IndFront, pdset::PrimaryDataset, carind::Int, validfind::Int)

	frameind = validfind2frameind(pdset, validfind)
	curlane = gete(pdset, :lane,  frameind)
	egoFx   = gete(pdset, :posFx, frameind)
	egoFy   = gete(pdset, :posFy, frameind)
	ncarsinframe = ncars_in_frame(pdset, validfind)

	searchcarinds = [0:ncarsinframe-1]
	if carind == CARIND_EGO

		if ncarsinframe == 0 # no other cars!
			meta[CARIND_EGO][:d_x_front] = NA_ALIAS
			meta[CARIND_EGO][:d_y_front] = NA_ALIAS
			meta[CARIND_EGO][:v_x_front] = NA_ALIAS
			meta[CARIND_EGO][:v_y_front] = NA_ALIAS
			meta[CARIND_EGO][:yaw_front] = NA_ALIAS
			meta[CARIND_EGO][:turnrate_front] = NA_ALIAS
			return NA_ALIAS
		end

		d_x_front, frontcar_ind = findmin(map(ci->begin
				dy    = getc(pdset, "posFy", ci, validfind) - egoFy
				dlane = getc(pdset, "lane",  ci, validfind) - curlane
				if isapprox(dlane, 0.0) || abs(dy) < THRESHOLD_DY_CONSIDERED_IN_FRONT
					dx = getc(pdset, "posFx", ci, validfind) - egoFx
					return dx > 0 ? dx : Inf
				end
				return Inf
			end, searchcarinds))

		if isinf(d_x_front)
			meta[CARIND_EGO][:d_x_front] = NA_ALIAS
			meta[CARIND_EGO][:d_y_front] = NA_ALIAS
			meta[CARIND_EGO][:v_x_front] = NA_ALIAS
			meta[CARIND_EGO][:v_y_front] = NA_ALIAS
			meta[CARIND_EGO][:yaw_front] = NA_ALIAS
			meta[CARIND_EGO][:turnrate_front] = NA_ALIAS
			return NA_ALIAS
		end

		@assert(!isa(d_x_front, NAtype))
		meta[CARIND_EGO][:d_x_front] = d_x_front
		meta[CARIND_EGO][:d_y_front] = getc(pdset, "posFy", searchcarinds[frontcar_ind], validfind) - egoFy
		meta[CARIND_EGO][:v_x_front] = getc(pdset, "velFx", searchcarinds[frontcar_ind], validfind) - gete(pdset, :velFx, frameind)
		meta[CARIND_EGO][:v_y_front] = getc(pdset, "velFy", searchcarinds[frontcar_ind], validfind) - gete(pdset, :velFy, frameind)
		meta[CARIND_EGO][:yaw_front] = getc(pdset, "posFyaw", searchcarinds[frontcar_ind], validfind)
		meta[CARIND_EGO][:turnrate_front] = get(TURNRATE, pdset, searchcarinds[frontcar_ind], validfind)
		return frontcar_ind
	end
	
	myFx   = getc(pdset, "posFx", carind, validfind)
	myFy   = getc(pdset, "posFy", carind, validfind)
	mylane = getc(pdset, "lane",  carind, validfind)
	frontcar_dist, frontcar_ind = findmin(map(carind2->begin
			if carind2 == carind
				return Inf
			end
			dy    = getc(pdset, "posFy", carind2, validfind) - myFy
			dlane = getc(pdset, "lane", carind2, validfind) -  mylane
			if isapprox(dlane, 0.0) || abs(dy) < THRESHOLD_DY_CONSIDERED_IN_FRONT
				dx = getc(pdset, "posFx", carind2, validfind) - myFx
				return dx > 0 ? dx : Inf
			end
			return Inf
		end, searchcarinds))
	dy_ego = myFy - egoFy
	dlane =  mylane - curlane
	if (isapprox(dlane, 0.0) || abs(dy_ego) < THRESHOLD_DY_CONSIDERED_IN_FRONT) &&
		egoFx - myFx > 0.0 &&
		egoFx - myFx < frontcar_dist

		# ego is better
		meta[carind][:d_x_front] = egoFx - myFx
		meta[carind][:d_y_front] = egoFy - myFy
		meta[carind][:v_x_front] = gete(pdset, :velFx, frameind) - getc(pdset, "velFx", carind, validfind)
		meta[carind][:v_y_front] = gete(pdset, :velFy, frameind) - getc(pdset, "velFy", carind, validfind)
		meta[carind][:yaw_front] = gete(pdset, :posFyaw, frameind)
		meta[carind][:turnrate_front] = get(TURNRATE, pdset, CARIND_EGO, validfind)
		return frontcar_ind
	end

	if isinf(frontcar_dist)
		meta[carind][:d_x_front] = NA_ALIAS
		meta[carind][:d_y_front] = NA_ALIAS
		meta[carind][:v_x_front] = NA_ALIAS
		meta[carind][:v_y_front] = NA_ALIAS
		meta[carind][:yaw_front] = NA_ALIAS
		meta[carind][:turnrate_front] = NA_ALIAS
		return NA_ALIAS
	end

	# other car is better
	meta[carind][:d_x_front] = frontcar_dist
	meta[carind][:d_y_front] = getc(pdset, "posFy", searchcarinds[frontcar_ind], validfind) - myFy
	meta[carind][:v_x_front] = getc(pdset, "velFx", searchcarinds[frontcar_ind], validfind) - getc(pdset, "velFx", carind, validfind)
	meta[carind][:v_y_front] = getc(pdset, "velFy", searchcarinds[frontcar_ind], validfind) - getc(pdset, "velFy", carind, validfind)
	meta[carind][:yaw_front] = getc(pdset, "posFyaw", searchcarinds[frontcar_ind], validfind)
	meta[carind][:turnrate_front] = get(TURNRATE, pdset, searchcarinds[frontcar_ind], validfind)
	frontcar_ind
end
function _get(::Feature_IndFront, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)

	best_carind = -2
	best_dist   = 100.0 # start at max
	best_ΔpFy   = 0.0

	frameind = validfind2frameind(pdset, validfind)

	ncarsinframe = ncars_in_frame(pdset, validfind)
	cars_to_check = Set([-1 : (ncarsinframe-1)])

	lanetags = Array(LaneTag, ncarsinframe+1)
	for cind in cars_to_check
		lanetags[cind+2] = get(pdset, :lanetag, "lanetag", cind, frameind, validfind)::LaneTag
	end

	posFx = get(pdset, :posFx, "posFx",  carind, frameind, validfind)::Float64
	posFy = get(pdset, :posFy, "posFy",  carind, frameind, validfind)::Float64
	active_lanetag = lanetags[carind+2]
	active_lane = get_lane(sn, active_lanetag)
	search_dist = 0.0

	delete!(cars_to_check, carind)

	done = false
	while !done
		to_remove = Set{Int}()
		for target_carind in cars_to_check
			if active_lanetag == lanetags[target_carind+2]
				target_posFx = get(pdset, :posFx, "posFx", target_carind, frameind, validfind)

				# only accept cars that are in front of us
				# and better than what we already have

				target_dist = target_posFx - posFx + search_dist
				if 0.0 < target_dist < best_dist

					target_posFy = get(pdset, :posFy, "posFy", target_carind, frameind, validfind)
					ΔpFy = target_posFy - posFy
					if abs(ΔpFy) < THRESHOLD_DY_CONSIDERED_IN_FRONT
						best_carind, best_dist = target_carind, target_dist
						best_ΔpFy = ΔpFy
					end
				end

				push!(to_remove, target_carind)
			end
		end
		if best_carind != -2
			break
		end

		for target_carind in to_remove
			delete!(cars_to_check, target_carind)
		end
		if isempty(cars_to_check)
			break
		end

		# move to next lane
		if has_next_lane(sn, active_lane)
			search_dist += active_lane.curve.s[end]
			active_lane = next_lane(sn, active_lane)
			active_lanetag = LaneTag(sn, active_lane)
			done = search_dist > best_dist
		else
			done = true
		end
	end

	if best_carind != -2
		meta[carind][:d_x_front] = best_dist
		meta[carind][:d_y_front] = best_ΔpFy
		meta[carind][:v_x_front] = get(pdset, :velFx,   "velFx",   best_carind, frameind, validfind) - get(pdset, :velFx,   "velFx", carind, frameind, validfind)
		meta[carind][:v_y_front] = get(pdset, :velFy,   "velFy",   best_carind, frameind, validfind) - get(pdset, :velFy,   "velFy", carind, frameind, validfind)
		meta[carind][:yaw_front] = get(pdset, :posFyaw, "posFyaw", best_carind, frameind, validfind)
		meta[carind][:turnrate_front] = get(TURNRATE, pdset, sn, best_carind, validfind)
		return float64(best_carind)
	else
		meta[carind][:d_x_front] = NA_ALIAS
		meta[carind][:d_y_front] = NA_ALIAS
		meta[carind][:v_x_front] = NA_ALIAS
		meta[carind][:v_y_front] = NA_ALIAS
		meta[carind][:yaw_front] = NA_ALIAS
		meta[carind][:turnrate_front] = NA_ALIAS
		return NA_ALIAS
	end
end

create_feature_basics( "HAS_FRONT", MarkovFeature, "-", true, true, 1.0, 0.0, false, :has_front, L"\exists_{fo}", "whether there is a car in front")
function get( ::Feature_HAS_FRONT, scene::TrafficScene, carind::Int)

	float64(get(INDFRONT, scene, carind) != NA_ALIAS)
end
function get( ::Feature_HAS_FRONT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	float64(get(INDFRONT, pdset, carind, validfind) != NA_ALIAS)
end
function get( ::Feature_HAS_FRONT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	
	float64(get(INDFRONT, pdset, sn, carind, validfind) != NA_ALIAS)
end

create_feature_basics( "D_X_FRONT", MarkovFeature, "m", false, false, Inf, 0.0, true, :d_x_front, L"d_{x,fo}", "longitudinal distance to the closest vehicle in the same lane in front")
function get( ::Feature_D_X_FRONT, scene::TrafficScene, carind::Int)
	get(INDFRONT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:d_x_front] # pull the processed result
end
function get( ::Feature_D_X_FRONT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDFRONT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:d_x_front] # pull the processed result
end
function get( ::Feature_D_X_FRONT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDFRONT, pdset, sn, carind, validfind) # call to get it to do the calculation
	meta[carind][:d_x_front] # pull the processed result
end

create_feature_basics( "D_Y_FRONT", MarkovFeature, "m", false, false, Inf, -Inf, true, :d_y_front, L"d_{y,fo}", "lateral distance to the closest vehicle in the same lane in front")
function get( ::Feature_D_Y_FRONT, scene::TrafficScene, carind::Int)
	get(INDFRONT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:d_y_front] # pull the processed result
end
function get( ::Feature_D_Y_FRONT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDFRONT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:d_y_front] # pull the processed result
end
function get( ::Feature_D_Y_FRONT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDFRONT, pdset, sn, carind, validfind) # call to get it to do the calculation
	meta[carind][:d_y_front] # pull the processed result
end

create_feature_basics( "V_X_FRONT", MarkovFeature, "m/s", false, false, Inf, -Inf, true, :v_x_front, L"v^{rel}_{x,fo}", "relative x velocity of the vehicle in front of you")
function get( ::Feature_V_X_FRONT, scene::TrafficScene, carind::Int)
	get(INDFRONT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:v_x_front] # pull the processed result
end
function get( ::Feature_V_X_FRONT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDFRONT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:v_x_front] # pull the processed result
end
function get( ::Feature_V_X_FRONT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDFRONT, pdset, sn, carind, validfind) # call to get it to do the calculation
	meta[carind][:v_x_front] # pull the processed result
end

create_feature_basics( "V_Y_FRONT", MarkovFeature, "m/s", false, false, Inf, -Inf, true, :v_y_front, L"v^{rel}_{y,fo}", "relative y velocity of the vehicle in front of you")
function get( ::Feature_V_Y_FRONT, scene::TrafficScene, carind::Int)
	get(INDFRONT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:v_y_front] # pull the processed result
end
function get( ::Feature_V_Y_FRONT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDFRONT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:v_y_front] # pull the processed result
end
function get( ::Feature_V_Y_FRONT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDFRONT, pdset, sn, carind, validfind) # call to get it to do the calculation
	meta[carind][:v_y_front] # pull the processed result
end

create_feature_basics( "YAW_FRONT", MarkovFeature, "rad", false, false, float64(pi), float64(-pi), true, :yaw_front, L"\psi_{fo}", "yaw of the vehicle in front of you")
function get( ::Feature_YAW_FRONT, scene::TrafficScene, carind::Int)
	get(INDFRONT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:yaw_front] # pull the processed result
end
function get( ::Feature_YAW_FRONT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDFRONT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:yaw_front] # pull the processed result
end
function get( ::Feature_YAW_FRONT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDFRONT, pdset, sn, carind, validfind) # call to get it to do the calculation
	meta[carind][:yaw_front] # pull the processed result
end

create_feature_basics( "TURNRATE_FRONT", MarkovFeature, "rad/s", false, false, Inf, -Inf, true, :turnrate_front, L"\dot{\psi}_{fo}", "turnrate of the vehicle in front of you")
function get( ::Feature_TURNRATE_FRONT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDFRONT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:turnrate_front] # pull the processed result
end
function get( ::Feature_TURNRATE_FRONT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDFRONT, pdset, sn, carind, validfind) # call to get it to do the calculation
	meta[carind][:turnrate_front] # pull the processed result
end
function get( ::Feature_TURNRATE_FRONT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDFRONT, pdset, sn, carind, validfind) # call to get it to do the calculation
	meta[carind][:turnrate_front] # pull the processed result
end

create_feature_basics( "A_REQ_FRONT", MarkovFeature, "m/s2", false, false, 0.0, -Inf, true, :a_req_front, L"a^{req}_{x,fo}", "const acceleration required to prevent collision with car in front assuming constant velocity")
function get( ::Feature_A_REQ_FRONT, scene::TrafficScene, carind::Int)

	# the constant acceleration required so as not to collide with the car in front
	# assumes both cars maintain constant vel
	# note: this assumes cars have zero-length, which is NOT the case
	# note: finds the acceleration such that we have zero velocity when reaching other car

	ind_front = get(INDFRONT, scene, carind)
	if ind_front == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_FRONT, scene, carind) # distance between cars
	dv = get(V_X_FRONT, scene, carind) # v_front - v_back

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	-min(dv*dv / (2dx),THRESHOLD_A_REQ)
end
function get( ::Feature_A_REQ_FRONT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	ind_front = get(INDFRONT, pdset, carind, validfind)
	if ind_front == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_FRONT, pdset, carind, validfind) # distance between cars
	dv = get(V_X_FRONT, pdset, carind, validfind) # v_front - v_back

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	-min(dv*dv / (2dx), THRESHOLD_A_REQ)
end
function get( ::Feature_A_REQ_FRONT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	
	ind_front = get(INDFRONT, pdset, sn, carind, validfind)
	if ind_front == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_FRONT, pdset, sn, carind, validfind) # distance between cars
	dv = get(V_X_FRONT, pdset, sn, carind, validfind) # v_front - v_back

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	-min(dv*dv / (2dx), THRESHOLD_A_REQ)
end

create_feature_basics( "Gaining_On_Front", MarkovFeature, "-", true, true, 1.0, 0.0, true, :gaining_on_front, L"1\{v_\text{ego} > v_\text{front}\}", "whether the car will collide with front if no action is taken and both const. vel")
function get( ::Feature_Gaining_On_Front, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDFRONT, pdset, sn, carind, validfind) # call to get it to do the calculation
	ΔV = meta[carind][:v_x_front] # pull the processed result
	float64(ΔV < 0.0)
end

create_feature_basics( "TTC_X_FRONT", MarkovFeature, "s", false, false, Inf, 0.0, true, :ttc_x_front, L"ttc_{x,fo}", "time to collision with car in front assuming constant velocities")
function get( ::Feature_TTC_X_FRONT, scene::TrafficScene, carind::Int)

	# assumes both cars maintain constant vel
	# note: this assumes cars have zero-length, which is NOT the case

	ind_front = get(INDFRONT, scene, carind)
	if ind_front == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_FRONT, scene, carind) # distance between cars
	dv = get(V_X_FRONT, scene, carind) # v_front - v_back

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	min(-dx / dv, THRESHOLD_TIME_TO_COLLISION)
end
function get( ::Feature_TTC_X_FRONT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	ind_front = get(INDFRONT, pdset, carind, validfind)
	if ind_front == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_FRONT, pdset, carind, validfind) # distance between cars
	dv = get(V_X_FRONT, pdset, carind, validfind) # v_front - v_back

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	min(-dx / dv, THRESHOLD_TIME_TO_COLLISION)
end
function get( ::Feature_TTC_X_FRONT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	
	ind_front = get(INDFRONT, pdset, sn, carind, validfind)
	if ind_front == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_FRONT, pdset, sn, carind, validfind) # distance between cars
	dv = get(V_X_FRONT, pdset, sn, carind, validfind) # v_front - v_back

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	min(-dx / dv, THRESHOLD_TIME_TO_COLLISION)
end

create_feature_basics( "TimeGap_X_FRONT", MarkovFeature, "s", false, false, Inf, 0.0, false, :timegap_x_front, L"\tau_{x,fo}", "timegap between cars")
function get( ::Feature_TimeGap_X_FRONT, scene::TrafficScene, carind::Int)

	# assumes we maintain constant velocity
	# computes time it will take for us to pass same location that the other car is at now

	ind_front = get(INDFRONT, scene, carind)
	if ind_front == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_FRONT, scene, carind) # distance between cars
	v  = get(VELFX,     scene, carind) # our current velocity

	if v <= 0.0
		return NA_ALIAS
	end

	dx / v
end
function get( ::Feature_TimeGap_X_FRONT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	ind_front = get(INDFRONT, pdset, carind, validfind)
	if ind_front == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_FRONT, pdset, carind, validfind) # distance between cars
	v  = get(VELFX,     pdset, carind, validfind) # our current velocity

	if v <= 0.0
		return NA_ALIAS
	end

	dx / v
end
function get( ::Feature_TimeGap_X_FRONT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	
	ind_front = get(INDFRONT, pdset, sn, carind, validfind)
	if ind_front == NA_ALIAS
		return Features.THRESHOLD_TIMEGAP
	end

	dx = get(D_X_FRONT, pdset, sn, carind, validfind) # distance between cars
	v  = get(VELFX,     pdset, sn, carind, validfind) # our current velocity

	if v <= 0.0
		return Features.THRESHOLD_TIMEGAP
	end

	min(dx / v, THRESHOLD_TIMEGAP)
end

create_feature_basics( "IndRear", MarkovFeature, "-", true, false, Inf, -2.0, true, :ind_rear, L"i_{rear}", "index of the closest car behind")
function _get(::Feature_IndRear, scene::TrafficScene, carind::Int)
	car = scene.cars[carind]

	mylane = get(CL, scene, carind)
	myFy   = car.pos.y
	myFx   = car.pos.x
	myVy   = car.spd*sin(car.pos.yaw)
	myVx   = car.spd*cos(car.pos.yaw)

	min_dx = Inf
	ind_rear = -2

	d_y_rear = 0.0
	v_x_rear = 0.0
	v_y_rear = 0.0
	yaw_rate = 0.0

	for (i,car2) in enumerate(scene.cars)
		if i == carind
			continue
		end

		dy = myFy - car2.pos.y
		dlane = get(CL, scene, i)
		if isapprox(dlane, mylane) || abs(dy) < THRESHOLD_DY_CONSIDERED_IN_FRONT
			dx = myFx - car2.pos.x
			if dx < min_dx && dx > 0.0
				min_dx = dx
				d_y_rear = myFy - car2.pos.y
				v_x_rear = myVx - car2.spd*cos(car2.pos.yaw)
				v_y_rear = myVy - car2.spd*sin(car2.pos.yaw)
				yaw_rear = car2.pos.yaw
				ind_rear = i
			end
		end
	end

	if isinf(min_dx)
		car.meta[:d_x_rear] = NA_ALIAS
		car.meta[:d_y_rear] = NA_ALIAS
		car.meta[:v_x_rear] = NA_ALIAS
		car.meta[:v_y_rear] = NA_ALIAS
		car.meta[:yaw_rear] = NA_ALIAS
		return NA_ALIAS
	end

	car.meta[:d_x_rear] = min_dx
	car.meta[:d_y_rear] = d_y_rear
	car.meta[:v_x_rear] = v_x_rear
	car.meta[:v_y_rear] = v_y_rear
	car.meta[:yaw_rear] = yaw_rear
	ind_rear
end
function _get(::Feature_IndRear, pdset::PrimaryDataset, carind::Int, validfind::Int)

	frameind = validfind2frameind(pdset, validfind)
	curlane = gete(pdset, :lane,  frameind)
	egoFx   = gete(pdset, :posFx, frameind)
	egoFy   = gete(pdset, :posFy, frameind)
	ncarsinframe = ncars_in_frame(pdset, validfind)

	searchcarinds = [0:ncarsinframe-1]
	if carind == CARIND_EGO

		if ncarsinframe == 0 # no other cars!
			meta[CARIND_EGO][:d_x_rear] = NA_ALIAS
			meta[CARIND_EGO][:d_y_rear] = NA_ALIAS
			meta[CARIND_EGO][:v_x_rear] = NA_ALIAS
			meta[CARIND_EGO][:v_y_rear] = NA_ALIAS
			meta[CARIND_EGO][:yaw_rear] = NA_ALIAS
			meta[CARIND_EGO][:turnrate_rear] = NA_ALIAS
			return NA_ALIAS
		end

		d_x_rear, rearcar_ind = findmin(map(ci->begin
				dy    = getc(pdset, "posFy", ci, validfind) - egoFy
				dlane = getc(pdset, "lane",  ci, validfind) - curlane
				if isapprox(dlane, 0.0) || abs(dy) < THRESHOLD_DY_CONSIDERED_IN_FRONT
					dx = egoFx - getc(pdset, "posFx", ci, validfind)
					return dx > 0.0 ? dx : Inf
				end
				return Inf
			end, searchcarinds))

		if isinf(d_x_rear)
			meta[CARIND_EGO][:d_x_rear] = NA_ALIAS
			meta[CARIND_EGO][:d_y_rear] = NA_ALIAS
			meta[CARIND_EGO][:v_x_rear] = NA_ALIAS
			meta[CARIND_EGO][:v_y_rear] = NA_ALIAS
			meta[CARIND_EGO][:yaw_rear] = NA_ALIAS
			meta[CARIND_EGO][:turnrate_rear] = NA_ALIAS
			return NA_ALIAS
		end

		@assert(!isa(d_x_rear, NAtype))
		meta[CARIND_EGO][:d_x_rear] = d_x_rear
		meta[CARIND_EGO][:d_y_rear] = egoFy - getc(pdset, "posFy", searchcarinds[rearcar_ind], validfind)
		meta[CARIND_EGO][:v_x_rear] = gete(pdset, :velFx, frameind) - getc(pdset, "velFx", searchcarinds[rearcar_ind], validfind)
		meta[CARIND_EGO][:v_y_rear] = gete(pdset, :velFy, frameind) - getc(pdset, "velFy", searchcarinds[rearcar_ind], validfind)
		meta[CARIND_EGO][:yaw_rear] = getc(pdset, "posFyaw", searchcarinds[rearcar_ind], validfind)
		meta[CARIND_EGO][:turnrate_rear] = get(TURNRATE, pdset, searchcarinds[rearcar_ind], validfind)
		return rearcar_ind
	end
	
	myFx   = getc(pdset, "posFx", carind, validfind)
	myFy   = getc(pdset, "posFy", carind, validfind)
	mylane = getc(pdset, "lane",  carind, validfind)
	rearcar_dist, rearcar_ind = findmin(map(carind2->begin
			if carind2 == carind
				return Inf
			end
			dy    = getc(pdset, "posFy", carind2, validfind) - myFy
			dlane = getc(pdset, "lane", carind2, validfind) -  mylane
			if isapprox(dlane, 0.0) || dy < THRESHOLD_DY_CONSIDERED_IN_FRONT
				dx = myFx - getc(pdset, "posFx", carind2, validfind)
				return dx > 0 ? dx : Inf
			end
			return Inf
		end, searchcarinds))
	dy_ego = myFy - egoFy
	dlane =  mylane - curlane
	if (isapprox(dlane, 0.0) || abs(dy_ego) < THRESHOLD_DY_CONSIDERED_IN_FRONT) &&
		myFx - egoFx > 0.0 &&
		myFx - egoFx < rearcar_dist

		# ego is better
		meta[carind][:d_x_rear] = myFx - egoFx
		meta[carind][:d_y_rear] = myFy - egoFy
		meta[carind][:v_x_rear] = getc(pdset, "velFx", carind, validfind) - gete(pdset, :velFx, frameind)
		meta[carind][:v_y_rear] = getc(pdset, "velFy", carind, validfind) - gete(pdset, :velFy, frameind)
		meta[carind][:yaw_rear] = gete(pdset, :posFyaw, frameind)
		meta[carind][:turnrate_rear] = get(TURNRATE, pdset, CARIND_EGO, validfind)
		return rearcar_ind
	end

	if isinf(rearcar_dist)
		meta[carind][:d_x_rear] = NA_ALIAS
		meta[carind][:d_y_rear] = NA_ALIAS
		meta[carind][:v_x_rear] = NA_ALIAS
		meta[carind][:v_y_rear] = NA_ALIAS
		meta[carind][:yaw_rear] = NA_ALIAS
		meta[carind][:turnrate_rear] = NA_ALIAS
		return NA_ALIAS
	end

	# other car is better
	meta[carind][:d_x_rear] = rearcar_dist
	meta[carind][:d_y_rear] = myFy                                    - getc(pdset, "posFy", searchcarinds[rearcar_ind], validfind)
	meta[carind][:v_x_rear] = getc(pdset, "velFx", carind, validfind) - getc(pdset, "velFx", searchcarinds[rearcar_ind], validfind)
	meta[carind][:v_y_rear] = getc(pdset, "velFy", carind, validfind) - getc(pdset, "velFy", searchcarinds[rearcar_ind], validfind)
	meta[carind][:yaw_rear] = getc(pdset, "posFyaw", searchcarinds[rearcar_ind], validfind)
	meta[carind][:turnrate_rear] = get(TURNRATE, pdset, searchcarinds[rearcar_ind], validfind)
	rearcar_ind
end
function _get(::Feature_IndRear, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)

	best_carind = -2
	best_dist   = 100.0 # start at max
	best_ΔpFy   = 0.0

	frameind = validfind2frameind(pdset, validfind)

	ncarsinframe = ncars_in_frame(pdset, validfind)
	cars_to_check = Set([-1 : (ncarsinframe-1)])

	lanetags = Array(LaneTag, ncarsinframe+1)
	for cind in cars_to_check
		lanetags[cind+2] = get(pdset, :lanetag, "lanetag", cind, frameind, validfind)::LaneTag
	end

	posFx = get(pdset, :posFx, "posFx",  carind, frameind, validfind)::Float64
	posFy = get(pdset, :posFy, "posFy",  carind, frameind, validfind)::Float64
	active_lanetag = lanetags[carind+2]
	active_lane = get_lane(sn, active_lanetag)
	search_dist = 0.0

	delete!(cars_to_check, carind)

	done = false
	while !done
		to_remove = Set{Int}()
		for target_carind in cars_to_check
			if active_lanetag == lanetags[target_carind+2]
				target_posFx = get(pdset, :posFx, "posFx", target_carind, frameind, validfind)

				target_dist = posFx - target_posFx + search_dist

				if 0.0 < target_dist < best_dist
					target_posFy = get(pdset, :posFy, "posFy", target_carind, frameind, validfind)
					ΔpFy = target_posFy - posFy
					if abs(ΔpFy) < THRESHOLD_DY_CONSIDERED_IN_FRONT
						best_carind, best_dist = target_carind, target_dist
						best_ΔpFy = ΔpFy
					end
				end

				push!(to_remove, target_carind)
			end
		end
		if best_carind != -2 || search_dist > best_dist
			break
		end

		for target_carind in to_remove
			delete!(cars_to_check, target_carind)
		end
		if isempty(cars_to_check)
			break
		end

		if has_prev_lane(sn, active_lane)
			active_lane = prev_lane(sn, active_lane)
			active_lanetag = LaneTag(sn, active_lane)
			search_dist += active_lane.curve.s[end]
		else
			done = true
		end
	end

	if best_carind != -2
		meta[carind][:d_x_rear] = best_dist
		meta[carind][:d_y_rear] = best_ΔpFy
		meta[carind][:v_x_rear] = get(pdset, :velFx,   "velFx",   best_carind, frameind, validfind) - get(pdset, :velFx,   "velFx", carind, frameind, validfind)
		meta[carind][:v_y_rear] = get(pdset, :velFy,   "velFy",   best_carind, frameind, validfind) - get(pdset, :velFy,   "velFy", carind, frameind, validfind)
		meta[carind][:yaw_rear] = get(pdset, :posFyaw, "posFyaw", best_carind, frameind, validfind)
		meta[carind][:turnrate_rear] = get(TURNRATE, pdset, sn, best_carind, validfind)
		return float64(best_carind)
	else
		meta[carind][:d_x_rear] = NA_ALIAS
		meta[carind][:d_y_rear] = NA_ALIAS
		meta[carind][:v_x_rear] = NA_ALIAS
		meta[carind][:v_y_rear] = NA_ALIAS
		meta[carind][:yaw_rear] = NA_ALIAS
		meta[carind][:turnrate_rear] = NA_ALIAS
		return NA_ALIAS
	end
end

create_feature_basics( "HAS_REAR", MarkovFeature, "-", true, true, 1.0, 0.0, false, :has_rear, L"\exists_{re}", "whether there is a car behind")
function get( ::Feature_HAS_REAR, scene::TrafficScene, carind::Int)

	float64(get(INDREAR, scene, carind) != NA_ALIAS)
end
function get( ::Feature_HAS_REAR, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	float64(get(INDREAR, pdset, carind, validfind) != NA_ALIAS)
end
function get( ::Feature_HAS_REAR, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	
	float64(get(INDREAR, pdset, sn, carind, validfind) != NA_ALIAS)
end

create_feature_basics( "D_X_REAR", MarkovFeature, "m", false, false, Inf, 0.0, true, :d_x_rear, L"d_{x,re}", "longitudinal distance to the closest vehicle in the same lane in rear")
function get( ::Feature_D_X_REAR, scene::TrafficScene, carind::Int)
	get(INDREAR, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:d_x_rear] # pull the processed result
end
function get( ::Feature_D_X_REAR, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDREAR, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:d_x_rear] # pull the processed result
end
function get( ::Feature_D_X_REAR, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDREAR, pdset, sn, carind, validfind)
	meta[carind][:d_x_rear]
end

create_feature_basics( "D_Y_REAR", MarkovFeature, "m", false, false, Inf, -Inf, true, :d_y_rear, L"d_{y,re}", "lateral distance to the closest vehicle in the same lane in rear")
function get( ::Feature_D_Y_REAR, scene::TrafficScene, carind::Int)
	get(INDREAR, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:d_y_rear] # pull the processed result
end
function get( ::Feature_D_Y_REAR, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDREAR, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:d_y_rear] # pull the processed result
end
function get( ::Feature_D_Y_REAR, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDREAR, pdset, sn, carind, validfind)
	meta[carind][:d_y_rear]
end

create_feature_basics( "V_X_REAR", MarkovFeature, "m/s", false, false, Inf, -Inf, true, :v_x_rear, L"v^{rel}_{x,re}", "relative x velocity of the vehicle behind you")
function get( ::Feature_V_X_REAR, scene::TrafficScene, carind::Int)
	get(INDREAR, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:v_x_rear] # pull the processed result
end
function get( ::Feature_V_X_REAR, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDREAR, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:v_x_rear] # pull the processed result
end
function get( ::Feature_V_X_REAR, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDREAR, pdset, sn, carind, validfind)
	meta[carind][:v_x_rear]
end

create_feature_basics( "V_Y_REAR", MarkovFeature, "m/s", false, false, Inf, -Inf, true, :v_y_rear, L"v^{rel}_{y,re}", "relative y velocity of the vehicle behind you")
function get( ::Feature_V_Y_REAR, scene::TrafficScene, carind::Int)
	get(INDREAR, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:v_y_rear] # pull the processed result
end
function get( ::Feature_V_Y_REAR, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDREAR, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:v_y_rear] # pull the processed result
end
function get( ::Feature_V_Y_REAR, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDREAR, pdset, sn, carind, validfind)
	meta[carind][:v_y_rear]
end

create_feature_basics( "YAW_REAR", MarkovFeature, "rad", false, false, float64(pi), float64(-pi), true, :yaw_rear, L"\psi^{rel}_{re}", "yaw of the vehicle behind you")
function get( ::Feature_YAW_REAR, scene::TrafficScene, carind::Int)
	get(INDREAR, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:yaw_rear] # pull the processed result
end
function get( ::Feature_YAW_REAR, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDREAR, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:yaw_rear] # pull the processed result
end
function get( ::Feature_YAW_REAR, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDREAR, pdset, sn, carind, validfind)
	meta[carind][:yaw_rear]
end

create_feature_basics( "TURNRATE_REAR", MarkovFeature, "rad/s", false, false, Inf, -Inf, true, :turnrate_rear, L"\dot{\psi}^{rel}_{re}", "turnrate of the vehicle behind you")
function get( ::Feature_TURNRATE_REAR, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDREAR, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:turnrate_rear] # pull the processed result
end
function get( ::Feature_TURNRATE_REAR, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDREAR, pdset, sn, carind, validfind)
	meta[carind][:turnrate_rear]
end

create_feature_basics( "A_REQ_REAR", MarkovFeature, "m/s2", false, false, Inf, 0.0, true, :a_req_rear, L"a^{req}_{x,re}", "const acceleration required to prevent collision with car behind assuming constant velocity")
function get( ::Feature_A_REQ_REAR, scene::TrafficScene, carind::Int)

	# the constant acceleration required so as not to collide with the other car
	# assumes both cars maintain constant vel
	# note: this assumes cars have zero-length, which is NOT the case

	ind_rear = get(INDREAR, scene, carind)
	if ind_rear == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_REAR, scene, carind) # distance between cars
	dv = get(V_X_REAR, scene, carind) # v_front - v_back

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	-min(dv*dv / (2dx), THRESHOLD_A_REQ)
end
function get( ::Feature_A_REQ_REAR, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	ind_rear = get(INDREAR, pdset, carind, validfind)
	if ind_rear == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_REAR, pdset, carind, validfind) # distance between cars
	dv = get(V_X_REAR, pdset, carind, validfind) # v_front - v_back

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	-(dv*dv / (2dx), THRESHOLD_A_REQ)
end
function get( ::Feature_A_REQ_REAR, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	
	ind_rear = get(INDREAR, pdset, sn, carind, validfind)
	if ind_rear == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_REAR, pdset, sn, carind, validfind) # distance between cars
	dv = get(V_X_REAR, pdset, sn, carind, validfind) # v_front - v_back

	if dv <= 0.0
		return NA_ALIAS
	end

	min(dv*dv / (2dx), THRESHOLD_A_REQ)
end

create_feature_basics( "TTC_X_REAR", MarkovFeature, "s", false, false, Inf, 0.0, true, :ttc_x_rear, L"ttc_{x,re}", "time to collision with rear car assuming constant velocities")
function get( ::Feature_TTC_X_REAR, scene::TrafficScene, carind::Int)

	# assumes both cars maintain constant vel
	# note: this assumes cars have zero-length, which is NOT the case

	ind_rear = get(INDREAR, scene, carind)
	if ind_rear == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_REAR, scene, carind) # distance between cars
	dv = get(V_X_REAR, scene, carind) # v_front - v_back

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	min(-dx / dv, THRESHOLD_TIME_TO_COLLISION)
end
function get( ::Feature_TTC_X_REAR, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	ind_rear = get(INDREAR, pdset, carind, validfind)
	if ind_rear == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_REAR, pdset, carind, validfind) # distance between cars
	dv = get(V_X_REAR, pdset, carind, validfind) # v_front - v_back

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	min(-dx / dv, THRESHOLD_TIME_TO_COLLISION)
end
function get( ::Feature_TTC_X_REAR, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	
	ind_rear = get(INDREAR, pdset, sn, carind, validfind)
	if ind_rear == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_REAR, pdset, sn, carind, validfind) # distance between cars
	dv = get(V_X_REAR, pdset, sn, carind, validfind) # v_them - v_us

	if dv <= 0.0
		return NA_ALIAS
	end

	min(dx / dv, THRESHOLD_TIME_TO_COLLISION)
end

create_feature_basics( "Timegap_X_REAR", MarkovFeature, "s", false, false, Inf, 0.0, true, :timegap_x_rear, L"\tau_{x,re}", "timegap with rear car")
function get( ::Feature_Timegap_X_REAR, scene::TrafficScene, carind::Int)

	ind_rear = get(INDREAR, scene, carind)
	if ind_rear == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_REAR, scene, carind) # distance between cars
	v  = get(VELFX,    scene, carind)

	if v <= 0.0
		return NA_ALIAS
	end

	dx / v
end
function get( ::Feature_Timegap_X_REAR, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	ind_rear = get(INDREAR, pdset, carind, validfind)
	if ind_rear == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_REAR, pdset, carind, validfind) # distance between cars
	v  = get(VELFX,    pdset, carind, validfind)

	if v <= 0.0
		return NA_ALIAS
	end

	dx / v
end
function get( ::Feature_Timegap_X_REAR, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	
	ind_rear = get(INDREAR, pdset, sn, carind, validfind)
	if ind_rear == NA_ALIAS
		return THRESHOLD_TIMEGAP
	end

	dx = get(D_X_REAR, pdset, sn, carind, validfind) # distance between cars
	v  = get(VELFX,    pdset, sn, carind, validfind)

	if v <= 0.0
		return THRESHOLD_TIMEGAP
	end

	min(dx / v, THRESHOLD_TIMEGAP)
end

create_feature_basics( "IndLeft", MarkovFeature, "-", true, false, Inf, -2.0, true, :ind_left, L"i_{\text{left}}", "index of the closest car in the left-hand lane")
function _get(::Feature_IndLeft, scene::TrafficScene, carind::Int)
	car = scene.cars[carind]

	mylane = get(CL, scene, carind)
	myFy   = car.pos.y
	myFx   = car.pos.x
	myVy   = car.spd*sin(car.pos.yaw)
	myVx   = car.spd*cos(car.pos.yaw)

	min_dx = Inf
	ind_left = -2

	d_y_left = 0.0
	v_x_left = 0.0
	v_y_left = 0.0
	yaw_left = 0.0

	for (i,car2) in enumerate(scene.cars)
		if i == carind
			continue
		end

		dy = car2.pos.y - myFy
		their_lane = get(CL, scene, i)
		if isapprox(their_lane, mylane+1)
			dx = car2.pos.x - myFx
			if abs(dx) < abs(min_dx)
				min_dx = dx
				d_y_left = car2.pos.y - myFy
				v_x_left = car2.spd*cos(car2.pos.yaw) - myVx
				v_y_left = car2.spd*sin(car2.pos.yaw) - myVy
				yaw_left = car2.pos.yaw
				ind_left = i
			end
		end
	end

	if isinf(min_dx)
		car.meta[:d_x_left] = NA_ALIAS
		car.meta[:d_y_left] = NA_ALIAS
		car.meta[:v_x_left] = NA_ALIAS
		car.meta[:v_y_left] = NA_ALIAS
		car.meta[:yaw_left] = NA_ALIAS
		return NA_ALIAS
	end

	car.meta[:d_x_left] = min_dx
	car.meta[:d_y_left] = d_y_left
	car.meta[:v_x_left] = v_x_left
	car.meta[:v_y_left] = v_y_left
	car.meta[:yaw_left] = yaw_left
	ind_left
end
function _get(::Feature_IndLeft, pdset::PrimaryDataset, carind::Int, validfind::Int)

	frameind = validfind2frameind(pdset, validfind)
	curlane = gete(pdset, :lane,  frameind)
	egoFx   = gete(pdset, :posFx, frameind)
	egoFy   = gete(pdset, :posFy, frameind)
	ncarsinframe = ncars_in_frame(pdset, validfind)

	searchcarinds = [0:ncarsinframe-1]
	if carind == CARIND_EGO

		if ncarsinframe == 0 # no other cars!
			meta[CARIND_EGO][:d_x_left] = NA_ALIAS
			meta[CARIND_EGO][:d_y_left] = NA_ALIAS
			meta[CARIND_EGO][:v_x_left] = NA_ALIAS
			meta[CARIND_EGO][:v_y_left] = NA_ALIAS
			meta[CARIND_EGO][:yaw_left] = NA_ALIAS
			meta[CARIND_EGO][:turnrate_left] = NA_ALIAS
			return NA_ALIAS
		end

		d_x_left, leftcar_ind = findmin(map(ci->begin
				dlane = getc(pdset, "lane",  ci, validfind) - curlane
				return isapprox(dlane, 1.0) ? abs(getc(pdset, "posFx", ci, validfind) - egoFx) : Inf
			end, searchcarinds))

		if isinf(d_x_left)
			meta[CARIND_EGO][:d_x_left] = NA_ALIAS
			meta[CARIND_EGO][:d_y_left] = NA_ALIAS
			meta[CARIND_EGO][:v_x_left] = NA_ALIAS
			meta[CARIND_EGO][:v_y_left] = NA_ALIAS
			meta[CARIND_EGO][:yaw_left] = NA_ALIAS
			meta[CARIND_EGO][:turnrate_left] = NA_ALIAS
			return NA_ALIAS
		end

		@assert(!isa(d_x_left, NAtype))
		meta[CARIND_EGO][:d_x_left] = getc(pdset, "posFx", searchcarinds[leftcar_ind], validfind) - egoFx
		meta[CARIND_EGO][:d_y_left] = getc(pdset, "posFy", searchcarinds[leftcar_ind], validfind) - egoFy
		meta[CARIND_EGO][:v_x_left] = getc(pdset, "velFx", searchcarinds[leftcar_ind], validfind) - gete(pdset, :velFx, frameind)
		meta[CARIND_EGO][:v_y_left] = getc(pdset, "velFy", searchcarinds[leftcar_ind], validfind) - gete(pdset, :velFy, frameind)
		meta[CARIND_EGO][:yaw_left] = getc(pdset, "posFyaw", searchcarinds[leftcar_ind], validfind)
		meta[CARIND_EGO][:turnrate_left] = get(TURNRATE, pdset, searchcarinds[leftcar_ind], validfind)
		return leftcar_ind
	end
	
	myFx   = getc(pdset, "posFx", carind, validfind)
	myFy   = getc(pdset, "posFy", carind, validfind)
	mylane = getc(pdset, "lane",  carind, validfind)
	leftcar_dist, leftcar_ind = findmin(map(carind2->begin
			if carind2 == carind
				return Inf
			end
			dlane = getc(pdset, "lane", carind2, validfind) -  mylane
			return isapprox(dlane, 1.0) ? abs(getc(pdset, "posFx", carind2, validfind) - myFx) : Inf
		end, searchcarinds))

	dy_ego = myFy - egoFy
	dlane =  mylane - curlane
	if isapprox(dlane, 1.0) && abs(egoFx - myFx) < leftcar_dist

		# ego is better
		meta[carind][:d_x_left] = egoFx - myFx
		meta[carind][:d_y_left] = egoFy - myFy
		meta[carind][:v_x_left] = gete(pdset, :velFx, frameind) - getc(pdset, "velFx", carind, validfind)
		meta[carind][:v_y_left] = gete(pdset, :velFy, frameind) - getc(pdset, "velFy", carind, validfind)
		meta[carind][:yaw_left] = gete(pdset, :posFyaw, frameind)
		meta[carind][:turnrate_left] = get(TURNRATE, pdset, CARIND_EGO, validfind)
		return leftcar_ind
	end

	if isinf(leftcar_dist)
		meta[carind][:d_x_left] = NA_ALIAS
		meta[carind][:d_y_left] = NA_ALIAS
		meta[carind][:v_x_left] = NA_ALIAS
		meta[carind][:v_y_left] = NA_ALIAS
		meta[carind][:yaw_left] = NA_ALIAS
		meta[carind][:turnrate_left] = NA_ALIAS
		return NA_ALIAS
	end

	# other car is better
	meta[carind][:d_x_left] = getc(pdset, "posFx", searchcarinds[leftcar_ind], validfind) - myFx
	meta[carind][:d_y_left] = getc(pdset, "posFy", searchcarinds[leftcar_ind], validfind) - myFy
	meta[carind][:v_x_left] = getc(pdset, "velFx", searchcarinds[leftcar_ind], validfind) - getc(pdset, "velFx", carind, validfind)
	meta[carind][:v_y_left] = getc(pdset, "velFy", searchcarinds[leftcar_ind], validfind) - getc(pdset, "velFy", carind, validfind)
	meta[carind][:yaw_left] = getc(pdset, "posFyaw", searchcarinds[leftcar_ind], validfind)
	meta[carind][:turnrate_left] = get(TURNRATE, pdset, searchcarinds[leftcar_ind], validfind)
	leftcar_ind
end
function _get(::Feature_IndLeft, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)

	# Finds the car that is closest in physical distance to the car in the left lane

	# TODO(tim): may still have an issue with computing distances across curve boundaries

	best_carind = -2
	best_dist   = 50.0
	best_ΔpFy   = NA_ALIAS

	frameind = validfind2frameind(pdset, validfind)

	if get(pdset, :nll, "nll",  carind, frameind, validfind)::Int8 > 0

		ncarsinframe = ncars_in_frame(pdset, validfind)
		cars_to_check = Set([-1 : (ncarsinframe-1)])

		lanetags = Array(LaneTag, ncarsinframe+1)
		for cind in cars_to_check
			lanetags[cind+2] = get(pdset, :lanetag, "lanetag", cind, frameind, validfind)::LaneTag
		end

		current_lanetag = lanetags[carind+2]
		current_lane    = get_lane(sn, current_lanetag)
		current_lane_next = has_next_lane(sn, current_lane) ? next_lane(sn, current_lane) : current_lane
		current_lane_prev = has_prev_lane(sn, current_lane) ? prev_lane(sn, current_lane) : current_lane

		posGx = get(pdset, :posGx,   "posGx",   carind, frameind, validfind)::Float64
		posGy = get(pdset, :posGy,   "posGy",   carind, frameind, validfind)::Float64
		posGθ = get(pdset, :posGyaw, "posGyaw", carind, frameind, validfind)::Float64

		rayEx = posGx + cos(posGθ)
		rayEy = posGy + sin(posGθ)

		# project the current location to the tilemap, but accept only lanes to the left of current location & not the current lane
		function f_filter(curve_pt::Vector{Float64}, lane::StreetLane)
			x = curve_pt[XIND]
			y = curve_pt[YIND]
			is_pt_left_of_ray(x, y, posGx, posGy, rayEx, rayEy) && 
				!(current_lane      === lane) &&
				!(current_lane_next === lane) && 
				!(current_lane_prev === lane)
		end
		proj = project_point_to_streetmap(posGx, posGy, sn, f_filter)

		if proj.successful

			posFx, posFy = pt_to_frenet_xy(proj.curvept, posGx, posGy)

			left_lanetag = LaneTag(proj.tile, proj.laneid)
			left_lane = get_lane(sn, left_lanetag)
			
			delete!(cars_to_check, carind)

			pq = Collections.PriorityQueue()
			Collections.enqueue!(pq, (left_lane, left_lanetag, true, 0.0), 0.0)

			done = false
			while !done && !isempty(pq)
				to_remove = Set{Int}()

				active_lane, active_lanetag, is_forward, search_dist = Collections.dequeue!(pq)

				# if ( is_forward && search_dist > abs(best_dist)) ||
				#    (!is_forward && search_dist > (abs(best_dist) + active_lane.curve.s[end]))
				#    break
				# end

				for target_carind in cars_to_check
					if active_lanetag == lanetags[target_carind+2]
						target_posFx = get(pdset, :posFx, "posFx", target_carind, frameind, validfind)

						target_dist = is_forward ? target_posFx - posFx + search_dist :
						                           posFx - target_posFx + search_dist

						if abs(target_dist) < abs(best_dist)
							best_carind, best_dist = target_carind, target_dist
							target_posFy = get(pdset, :posFy, "posFy", target_carind, frameind, validfind)
							best_ΔpFy = target_posFy - posFy
						end

						push!(to_remove, target_carind)
					end
				end	

				for target_carind in to_remove
					delete!(cars_to_check, target_carind)
				end
				if isempty(cars_to_check)
					break
				end

				if is_forward && has_next_lane(sn, active_lane)
					next_search_dist = search_dist + active_lane.curve.s[end]
					next_active_lane = next_lane(sn, active_lane)
					Collections.enqueue!(pq, (next_active_lane, LaneTag(sn, next_active_lane), true, next_search_dist), next_search_dist)
				end
				if (!is_forward || isapprox(search_dist, 0.0)) && has_prev_lane(sn, active_lane)
					prev_active_lane = prev_lane(sn, active_lane)
					prev_search_dist = search_dist + prev_active_lane.curve.s[end]
					Collections.enqueue!(pq, (prev_active_lane, LaneTag(sn, prev_active_lane), false, prev_search_dist), prev_search_dist)
				end
			end
		end
	end

	if best_carind != -2
		meta[carind][:d_x_left] = best_dist # NOTE(tim): + if other car in front, - if other car behind
		meta[carind][:d_y_left] = best_ΔpFy
		meta[carind][:v_x_left] = get(pdset, :velFx,   "velFx",   best_carind, frameind, validfind) - get(pdset, :velFx,   "velFx", carind, frameind, validfind)
		meta[carind][:v_y_left] = get(pdset, :velFy,   "velFy",   best_carind, frameind, validfind) - get(pdset, :velFy,   "velFy", carind, frameind, validfind)
		meta[carind][:yaw_left] = get(pdset, :posFyaw, "posFyaw", best_carind, frameind, validfind)
		meta[carind][:turnrate_left] = get(TURNRATE, pdset, sn, best_carind, validfind)
		return float64(best_carind)
	else
		meta[carind][:d_x_left] = NA_ALIAS
		meta[carind][:d_y_left] = NA_ALIAS
		meta[carind][:v_x_left] = NA_ALIAS
		meta[carind][:v_y_left] = NA_ALIAS
		meta[carind][:yaw_left] = NA_ALIAS
		meta[carind][:turnrate_left] = NA_ALIAS
		return NA_ALIAS
	end
end

create_feature_basics( "D_X_LEFT", MarkovFeature, "m", false, false, Inf, -Inf, true, :d_x_left, L"d_{x,lf}", "longitudinal distance to the closest vehicle in the left lane")
function get( ::Feature_D_X_LEFT, scene::TrafficScene, carind::Int)
	get(INDLEFT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:d_x_left] # pull the processed result
end
function get( ::Feature_D_X_LEFT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDLEFT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:d_x_left] # pull the processed result
end
function get( ::Feature_D_X_LEFT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDLEFT, pdset, sn, carind, validfind)
	meta[carind][:d_x_left]
end

create_feature_basics( "D_Y_LEFT", MarkovFeature, "m", false, false, Inf, -Inf, true, :d_y_left, L"d_{y,lf}", "lateral distance to the closest vehicle in the left lane")
function get( ::Feature_D_Y_LEFT, scene::TrafficScene, carind::Int)
	get(INDLEFT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:d_y_left] # pull the processed result
end
function get( ::Feature_D_Y_LEFT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDLEFT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:d_y_left] # pull the processed result
end
function get( ::Feature_D_Y_LEFT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDLEFT, pdset, sn, carind, validfind)
	meta[carind][:d_y_left]
end

create_feature_basics( "V_X_LEFT", MarkovFeature, "m/s", false, false, Inf, -Inf, true, :v_x_left, L"v^{rel}_{x,lf}", "relative x velocity of the closest vehicle in left lane")
function get( ::Feature_V_X_LEFT, scene::TrafficScene, carind::Int)
	get(INDLEFT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:v_x_left] # pull the processed result
end
function get( ::Feature_V_X_LEFT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDLEFT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:v_x_left] # pull the processed result
end
function get( ::Feature_V_X_LEFT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDLEFT, pdset, sn, carind, validfind)
	meta[carind][:v_x_left]
end

create_feature_basics( "V_Y_LEFT", MarkovFeature, "m/s", false, false, Inf, -Inf, true, :v_y_left, L"v^{rel}_{y,le}", "relative y velocity of the closest vehicle in left lane")
function get( ::Feature_V_Y_LEFT, scene::TrafficScene, carind::Int)
	get(INDLEFT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:v_y_left] # pull the processed result
end
function get( ::Feature_V_Y_LEFT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDLEFT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:v_y_left] # pull the processed result
end
function get( ::Feature_V_Y_LEFT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDLEFT, pdset, sn, carind, validfind)
	meta[carind][:v_y_left]
end

create_feature_basics( "YAW_LEFT", MarkovFeature, "rad", false, false, float64(pi), float64(-pi), true, :yaw_left, L"\psi_{le}", "yaw of the closest vehicle in left lane")
function get( ::Feature_YAW_LEFT, scene::TrafficScene, carind::Int)
	get(INDLEFT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:yaw_left] # pull the processed result
end
function get( ::Feature_YAW_LEFT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDLEFT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:yaw_left] # pull the processed result
end
function get( ::Feature_YAW_LEFT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDLEFT, pdset, sn, carind, validfind)
	meta[carind][:yaw_left]
end

create_feature_basics( "TURNRATE_LEFT", MarkovFeature, "rad/s", false, false, Inf, -Inf, true, :turnrate_left, L"\dot{\psi}_{le}", "turnrate of the closest vehicle in left lane")
function get( ::Feature_TURNRATE_LEFT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDLEFT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:turnrate_left] # pull the processed result
end
function get( ::Feature_TURNRATE_LEFT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDLEFT, pdset, sn, carind, validfind)
	meta[carind][:turnrate_left]
end

create_feature_basics( "A_REQ_LEFT", MarkovFeature, "m/s2", false, false, 0.0, -Inf, true, :a_req_left, L"a^{req}_{x,le}", "const acceleration (+ to left) required to prevent collision with car to left assuming constant velocity")
function get( ::Feature_A_REQ_LEFT, scene::TrafficScene, carind::Int)

	# the constant acceleration required so as not to collide with the other car
	# assumes both cars maintain constant vel
	# note: this assumes cars have zero-length, which is NOT the case

	ind_left = get(INDLEFT, scene, carind)
	if ind_left == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_LEFT, scene, carind) # distance between cars
	dv = get(V_X_LEFT, scene, carind) # v_other - v_me (+ if I am going to run into them)

	if dv <= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	-min(dv*dv / (2dx), THRESHOLD_A_REQ)
end
function get( ::Feature_A_REQ_LEFT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	ind_left = get(INDLEFT, pdset, carind, validfind)
	if ind_left == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_LEFT, pdset, carind, validfind) # distance between cars
	dv = get(V_X_LEFT, pdset, carind, validfind) # v_other - v_me

	if dv <= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	-min(dv*dv / (2dx), THRESHOLD_A_REQ)
end
function get( ::Feature_A_REQ_LEFT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	
	ind_left = get(INDLEFT, pdset, sn, carind, validfind)
	if ind_left == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_LEFT, pdset, sn, carind, validfind) # distance between cars
	dv = get(V_X_LEFT, pdset, sn, carind, validfind) # v_other - v_me

	if (dx > 0.0 && dv > 0.0) || (dx < 0.0 && dv < 0.0)
		return NA_ALIAS
	end

	if dx > 0.0
		-min(dv*dv / (2*dx), THRESHOLD_A_REQ)
	else
		min(dv*dv / (2*abs(dx)), THRESHOLD_A_REQ)
	end
end

create_feature_basics( "TTC_X_LEFT", MarkovFeature, "s", false, false, Inf, 0.0, true, :ttc_x_left, L"ttc_{x,le}", "time to collision with left car assuming constant velocities")
function get( ::Feature_TTC_X_LEFT, scene::TrafficScene, carind::Int)

	# assumes both cars maintain constant vel
	# note: this assumes cars have zero-length, which is NOT the case

	ind_left = get(INDLEFT, scene, carind)
	if ind_left == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_LEFT, scene, carind) # distance between cars
	dv = get(V_X_LEFT, scene, carind) # v_other - v_me (+ if I am going to run into them)

	if dv <= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	min(dx / dv, THRESHOLD_TIME_TO_COLLISION)
end
function get( ::Feature_TTC_X_LEFT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	ind_left = get(INDLEFT, pdset, carind, validfind)
	if ind_left == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_LEFT, pdset, carind, validfind) # distance between cars
	dv = get(V_X_LEFT, pdset, carind, validfind) # v_other - v_me

	if dv <= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	min(dx / dv, THRESHOLD_TIME_TO_COLLISION)
end
function get( ::Feature_TTC_X_LEFT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	
	ind_left = get(INDLEFT, pdset, sn, carind, validfind)
	if ind_left == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_LEFT, pdset, sn, carind, validfind) # distance between cars
	dv = get(V_X_LEFT, pdset, sn, carind, validfind) # v_other - v_me

	if (dx > 0.0 && dv > 0.0) || (dx < 0.0 && dv < 0.0)
		return NA_ALIAS
	end

	min(abs(dx / dv), THRESHOLD_TIME_TO_COLLISION)
end

create_feature_basics( "Timegap_X_LEFT", MarkovFeature, "s", false, false, Inf, 0.0, true, :timegap_x_left, L"\tau_{x,le}", "timegap with left car")
function get( ::Feature_Timegap_X_LEFT, scene::TrafficScene, carind::Int)

	ind_left = get(INDLEFT, scene, carind)
	if ind_left == NA_ALIAS
		return THRESHOLD_TIMEGAP
	end

	dx = get(D_X_LEFT, scene, carind) # distance between cars
	v  = get(VELFX,    scene, carind)

	if v <= 0.0
		return THRESHOLD_TIMEGAP
	end

	dx / v
end
function get( ::Feature_Timegap_X_LEFT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	ind_left = get(INDLEFT, pdset, carind, validfind)
	if ind_left == NA_ALIAS
		return THRESHOLD_TIMEGAP
	end

	dx = get(D_X_LEFT, pdset, carind, validfind) # distance between cars
	v  = get(VELFX,    pdset, carind, validfind)

	if v <= 0.0
		return THRESHOLD_TIMEGAP
	end

	dx / v
end
function get( ::Feature_Timegap_X_LEFT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)

	ind_left = get(INDLEFT, pdset, sn, carind, validfind)
	if ind_left == NA_ALIAS
		return THRESHOLD_TIMEGAP
	end

	dx = get(D_X_LEFT, pdset, sn, carind, validfind) # distance between cars
	v  = get(VELFX,    pdset, sn, carind, validfind)

	if (dx > 0.0 && v > 0.0) || (dx < 0.0 && v < 0.0)
		return THRESHOLD_TIMEGAP
	end

	min(abs(dx / v), THRESHOLD_TIMEGAP)
end

create_feature_basics( "IndRight", MarkovFeature, "-", true, false, Inf, -2.0, true, :ind_right, L"i_{right}", "index of the closest car in the right-hand lane")
function _get(::Feature_IndRight, scene::TrafficScene, carind::Int)
	car = scene.cars[carind]

	mylane = get(CL, scene, carind)
	myFy   = car.pos.y
	myFx   = car.pos.x
	myVy   = car.spd*sin(car.pos.yaw)
	myVx   = car.spd*cos(car.pos.yaw)

	min_dx = Inf
	ind_right = -2

	d_y_right = 0.0
	v_x_right = 0.0
	v_y_right = 0.0
	yaw_right = 0.0

	for (i,car2) in enumerate(scene.cars)
		if i == carind
			continue
		end

		dy = car2.pos.y - myFy
		dlane = get(CL, scene, i)
		if isapprox(dlane, mylane-1)
			dx = car2.pos.x - myFx
			if abs(dx) < abs(min_dx)
				min_dx = dx
				d_y_right = car2.pos.y - myFy
				v_x_right = car2.spd*cos(car2.pos.yaw) - myVx
				v_y_right = car2.spd*sin(car2.pos.yaw) - myVy
				yaw_right = car2.pos.yaw
				ind_right = i
			end
		end
	end

	if isinf(min_dx)
		car.meta[:d_x_right] = NA_ALIAS
		car.meta[:d_y_right] = NA_ALIAS
		car.meta[:v_x_right] = NA_ALIAS
		car.meta[:v_y_right] = NA_ALIAS
		car.meta[:yaw_right] = NA_ALIAS
		return NA_ALIAS
	end

	car.meta[:d_x_right] = min_dx
	car.meta[:d_y_right] = d_y_right
	car.meta[:v_x_right] = v_x_right
	car.meta[:v_y_right] = v_y_right
	car.meta[:yaw_right] = yaw_right
	ind_right
end
function _get(::Feature_IndRight, pdset::PrimaryDataset, carind::Int, validfind::Int)

	frameind = validfind2frameind(pdset, validfind)
	curlane = gete(pdset, :lane,  frameind)
	egoFx   = gete(pdset, :posFx, frameind)
	egoFy   = gete(pdset, :posFy, frameind)
	ncarsinframe = ncars_in_frame(pdset, validfind)

	searchcarinds = [0:ncarsinframe-1]
	if carind == CARIND_EGO

		if ncarsinframe == 0 # no other cars!
			meta[CARIND_EGO][:d_x_right] = NA_ALIAS
			meta[CARIND_EGO][:d_y_right] = NA_ALIAS
			meta[CARIND_EGO][:v_x_right] = NA_ALIAS
			meta[CARIND_EGO][:v_y_right] = NA_ALIAS
			meta[CARIND_EGO][:yaw_right] = NA_ALIAS
			meta[CARIND_EGO][:turnrate_right] = NA_ALIAS
			return NA_ALIAS
		end

		d_x_right, rightcar_ind = findmin(map(ci->begin
				dlane = getc(pdset, "lane",  ci, validfind) - curlane
				return isapprox(dlane, -1.0) ? abs(getc(pdset, "posFx", ci, validfind) - egoFx) : Inf
			end, searchcarinds))

		if isinf(d_x_right)
			meta[CARIND_EGO][:d_x_right] = NA_ALIAS
			meta[CARIND_EGO][:d_y_right] = NA_ALIAS
			meta[CARIND_EGO][:v_x_right] = NA_ALIAS
			meta[CARIND_EGO][:v_y_right] = NA_ALIAS
			meta[CARIND_EGO][:yaw_right] = NA_ALIAS
			meta[CARIND_EGO][:turnrate_right] = NA_ALIAS
			return NA_ALIAS
		end

		@assert(!isa(d_x_right, NAtype))
		meta[CARIND_EGO][:d_x_right] = getc(pdset, "posFx", searchcarinds[rightcar_ind], validfind) - egoFx
		meta[CARIND_EGO][:d_y_right] = getc(pdset, "posFy", searchcarinds[rightcar_ind], validfind) - egoFy
		meta[CARIND_EGO][:v_x_right] = getc(pdset, "velFx", searchcarinds[rightcar_ind], validfind) - gete(pdset, :velFx, frameind)
		meta[CARIND_EGO][:v_y_right] = getc(pdset, "velFy", searchcarinds[rightcar_ind], validfind) - gete(pdset, :velFy, frameind)
		meta[CARIND_EGO][:yaw_right] = getc(pdset, "posFyaw", searchcarinds[rightcar_ind], validfind)
		meta[CARIND_EGO][:turnrate_right] = get(TURNRATE, pdset, searchcarinds[rightcar_ind], validfind)
		return rightcar_ind
	end
	
	myFx   = getc(pdset, "posFx", carind, validfind)
	myFy   = getc(pdset, "posFy", carind, validfind)
	mylane = getc(pdset, "lane",  carind, validfind)
	rightcar_dist, rightcar_ind = findmin(map(carind2->begin
			if carind2 == carind
				return Inf
			end
			dlane = getc(pdset, "lane", carind2, validfind) -  mylane
			return isapprox(dlane, -1.0) ? abs(getc(pdset, "posFx", carind2, validfind) - myFx) : Inf
		end, searchcarinds))

	dy_ego = myFy - egoFy
	dlane =  mylane - curlane
	if isapprox(dlane, -1.0) && abs(egoFx - myFx) < rightcar_dist

		# ego is better
		meta[carind][:d_x_right] = egoFx - myFx
		meta[carind][:d_y_right] = egoFy - myFy
		meta[carind][:v_x_right] = gete(pdset, :velFx, frameind) - getc(pdset, "velFx", carind, validfind)
		meta[carind][:v_y_right] = gete(pdset, :velFy, frameind) - getc(pdset, "velFy", carind, validfind)
		meta[carind][:yaw_right] = gete(pdset, :posFyaw, frameind)
		meta[carind][:turnrate_right] = get(TURNRATE, pdset, CARIND_EGO, validfind)
		return rightcar_ind
	end

	if isinf(rightcar_dist)
		meta[carind][:d_x_right] = NA_ALIAS
		meta[carind][:d_y_right] = NA_ALIAS
		meta[carind][:v_x_right] = NA_ALIAS
		meta[carind][:v_y_right] = NA_ALIAS
		meta[carind][:yaw_right] = NA_ALIAS
		meta[carind][:turnrate_right] = NA_ALIAS
		return NA_ALIAS
	end

	# other car is better
	meta[carind][:d_x_right] = getc(pdset, "posFx", searchcarinds[rightcar_ind], validfind) - myFx
	meta[carind][:d_y_right] = getc(pdset, "posFy", searchcarinds[rightcar_ind], validfind) - myFy
	meta[carind][:v_x_right] = getc(pdset, "velFx", searchcarinds[rightcar_ind], validfind) - getc(pdset, "velFx", carind, validfind)
	meta[carind][:v_y_right] = getc(pdset, "velFy", searchcarinds[rightcar_ind], validfind) - getc(pdset, "velFy", carind, validfind)
	meta[carind][:yaw_right] = getc(pdset, "posFyaw", searchcarinds[rightcar_ind], validfind)
	meta[carind][:turnrate_right] = get(TURNRATE, pdset, searchcarinds[rightcar_ind], validfind)
	rightcar_ind
end
function _get(::Feature_IndRight, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)

	# TODO(tim): may still have an issue with computing distances across curve boundaries

	best_carind = -2
	best_dist   = 50.0
	best_ΔpFy   = NA_ALIAS

	frameind = validfind2frameind(pdset, validfind)

	if get(pdset, :nlr, "nlr",  carind, frameind, validfind)::Int8 > 0

		ncarsinframe = ncars_in_frame(pdset, validfind)
		cars_to_check = Set([-1 : (ncarsinframe-1)])

		lanetags = Array(LaneTag, ncarsinframe+1)
		for cind in cars_to_check
			lanetags[cind+2] = get(pdset, :lanetag, "lanetag", cind, frameind, validfind)::LaneTag
		end

		current_lanetag = lanetags[carind+2]
		current_lane    = get_lane(sn, current_lanetag)
		current_lane_next = has_next_lane(sn, current_lane) ? next_lane(sn, current_lane) : current_lane
		current_lane_prev = has_prev_lane(sn, current_lane) ? prev_lane(sn, current_lane) : current_lane

		posGx = get(pdset, :posGx,   "posGx",   carind, frameind, validfind)::Float64
		posGy = get(pdset, :posGy,   "posGy",   carind, frameind, validfind)::Float64
		posGθ = get(pdset, :posGyaw, "posGyaw", carind, frameind, validfind)::Float64

		rayEx = posGx + cos(posGθ)
		rayEy = posGy + sin(posGθ)

		function f_filter(curve_pt::Vector{Float64}, lane::StreetLane)
			x = curve_pt[XIND]
			y = curve_pt[YIND]
			!is_pt_left_of_ray(x, y, posGx, posGy, rayEx, rayEy) && 
				!(current_lane      === lane) &&
				!(current_lane_next === lane) && 
				!(current_lane_prev === lane)
		end
		proj = project_point_to_streetmap(posGx, posGy, sn, f_filter)

		if proj.successful

			posFx, posFy = pt_to_frenet_xy(proj.curvept, posGx, posGy)

			right_lanetag = LaneTag(proj.tile, proj.laneid)
			right_lane = get_lane(sn, right_lanetag)
			
			delete!(cars_to_check, carind)

			pq = Collections.PriorityQueue()
			Collections.enqueue!(pq, (right_lane, right_lanetag, true, 0.0), 0.0)

			done = false
			while !done && !isempty(pq)
				to_remove = Set{Int}()

				active_lane, active_lanetag, is_forward, search_dist = Collections.dequeue!(pq)

				for target_carind in cars_to_check
					if active_lanetag == lanetags[target_carind+2]
						target_posFx = get(pdset, :posFx, "posFx", target_carind, frameind, validfind)

						target_dist = is_forward ? target_posFx - posFx + search_dist :
						                           posFx - target_posFx + search_dist

						if abs(target_dist) < abs(best_dist)
							best_carind, best_dist = target_carind, target_dist
							target_posFy = get(pdset, :posFy, "posFy", target_carind, frameind, validfind)
							best_ΔpFy = target_posFy - posFy
						end

						push!(to_remove, target_carind)
					end
				end	

				for target_carind in to_remove
					delete!(cars_to_check, target_carind)
				end
				if isempty(cars_to_check)
					break
				end

				if is_forward && has_next_lane(sn, active_lane)
					next_search_dist = search_dist + active_lane.curve.s[end]
					next_active_lane = next_lane(sn, active_lane)
					Collections.enqueue!(pq, (next_active_lane, LaneTag(sn, next_active_lane), true, next_search_dist), next_search_dist)
				end
				if (!is_forward || isapprox(search_dist, 0.0)) && has_prev_lane(sn, active_lane)
					prev_active_lane = prev_lane(sn, active_lane)
					prev_search_dist = search_dist + prev_active_lane.curve.s[end]
					Collections.enqueue!(pq, (prev_active_lane, LaneTag(sn, prev_active_lane), false, prev_search_dist), prev_search_dist)
				end
			end
		end
	end

	if best_carind != -2
		meta[carind][:d_x_right] = best_dist # NOTE(tim): + if other car in front, - if other car behind
		meta[carind][:d_y_right] = best_ΔpFy
		meta[carind][:v_x_right] = get(pdset, :velFx,   "velFx",   best_carind, frameind, validfind) - get(pdset, :velFx,   "velFx", carind, frameind, validfind)
		meta[carind][:v_y_right] = get(pdset, :velFy,   "velFy",   best_carind, frameind, validfind) - get(pdset, :velFy,   "velFy", carind, frameind, validfind)
		meta[carind][:yaw_right] = get(pdset, :posFyaw, "posFyaw", best_carind, frameind, validfind)
		meta[carind][:turnrate_right] = get(TURNRATE, pdset, sn, best_carind, validfind)
		return float64(best_carind)
	else
		meta[carind][:d_x_right] = NA_ALIAS
		meta[carind][:d_y_right] = NA_ALIAS
		meta[carind][:v_x_right] = NA_ALIAS
		meta[carind][:v_y_right] = NA_ALIAS
		meta[carind][:yaw_right] = NA_ALIAS
		meta[carind][:turnrate_right] = NA_ALIAS
		return NA_ALIAS
	end
end

create_feature_basics( "D_X_RIGHT", MarkovFeature, "m", false, false, Inf, -Inf, true, :d_x_right, L"d_{x,ri}", "longitudinal distance to the closest vehicle in the right lane")
function get( ::Feature_D_X_RIGHT, scene::TrafficScene, carind::Int)
	get(INDRIGHT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:d_x_right] # pull the processed result
end
function get( ::Feature_D_X_RIGHT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDRIGHT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:d_x_right] # pull the processed result
end
function get( ::Feature_D_X_RIGHT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDRIGHT, pdset, sn, carind, validfind)
	meta[carind][:d_x_right]
end

create_feature_basics( "D_Y_RIGHT", MarkovFeature, "m", false, false, Inf, -Inf, true, :d_y_right, L"d_{y,ri}", "lateral distance to the closest vehicle in the right lane")
function get( ::Feature_D_Y_RIGHT, scene::TrafficScene, carind::Int)
	get(INDRIGHT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:d_y_right] # pull the processed result
end
function get( ::Feature_D_Y_RIGHT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDRIGHT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:d_y_right] # pull the processed result
end
function get( ::Feature_D_Y_RIGHT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDRIGHT, pdset, sn, carind, validfind)
	meta[carind][:d_y_right]
end

create_feature_basics( "V_X_RIGHT", MarkovFeature, "m/s", false, false, Inf, -Inf, true, :v_x_right, L"v^{rel}_{x,ri}", "relative x velocity of the closest vehicle in right lane")
function get( ::Feature_V_X_RIGHT, scene::TrafficScene, carind::Int)
	get(INDRIGHT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:v_x_right] # pull the processed result
end
function get( ::Feature_V_X_RIGHT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDRIGHT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:v_x_right] # pull the processed result
end
function get( ::Feature_V_X_RIGHT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDRIGHT, pdset, sn, carind, validfind)
	meta[carind][:v_x_right]
end

create_feature_basics( "V_Y_RIGHT", MarkovFeature, "m/s", false, false, Inf, -Inf, true, :v_y_right, L"v^{rel}_{y,ri}", "relative y velocity of the closest vehicle in right lane")
function get( ::Feature_V_Y_RIGHT, scene::TrafficScene, carind::Int)
	get(INDRIGHT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:v_y_right] # pull the processed result
end
function get( ::Feature_V_Y_RIGHT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDRIGHT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:v_y_right] # pull the processed result
end
function get( ::Feature_V_Y_RIGHT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDRIGHT, pdset, sn, carind, validfind)
	meta[carind][:v_y_right]
end

create_feature_basics( "YAW_RIGHT", MarkovFeature, "rad", false, false, float64(pi), float64(-pi), true, :yaw_right, L"\psi_{ri}", "yaw of the closest vehicle in right lane")
function get( ::Feature_YAW_RIGHT, scene::TrafficScene, carind::Int)
	get(INDRIGHT, scene, carind) # call to get it to do the calculation
	scene.cars[carind].meta[:yaw_right] # pull the processed result
end
function get( ::Feature_YAW_RIGHT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDRIGHT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:yaw_right] # pull the processed result
end
function get( ::Feature_YAW_RIGHT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDRIGHT, pdset, sn, carind, validfind)
	meta[carind][:yaw_right]
end

create_feature_basics( "TURNRATE_RIGHT", MarkovFeature, "rad", false, false, Inf, -Inf, true, :turnrate_right, L"\dot{\psi}_{ri}", "turnrate of the closest vehicle in right lane")
function get( ::Feature_TURNRATE_RIGHT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	get(INDRIGHT, pdset, carind, validfind) # call to get it to do the calculation
	meta[carind][:turnrate_right] # pull the processed result
end
function get( ::Feature_TURNRATE_RIGHT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	get(INDRIGHT, pdset, sn, carind, validfind)
	meta[carind][:turnrate_right]
end

create_feature_basics( "A_REQ_RIGHT", MarkovFeature, "m/s2", false, false, Inf, 0.0, true, :a_req_right, L"a^{req}_{x,ri}", "const acceleration (+ to left) required to prevent collision with car to right assuming constant velocity")
function get( ::Feature_A_REQ_RIGHT, scene::TrafficScene, carind::Int)

	# the constant acceleration required so as not to collide with the other car
	# assumes both cars maintain constant vel
	# note: this assumes cars have zero-length, which is NOT the case

	ind_right = get(INDRIGHT, scene, carind)
	if ind_right == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_RIGHT, scene, carind) # distance between cars
	dv = get(V_X_RIGHT, scene, carind) # v_other - v_me (- if I am going to run into them)

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	-min(dv*dv / (2dx), THRESHOLD_A_REQ)
end
function get( ::Feature_A_REQ_RIGHT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	ind_right = get(INDRIGHT, pdset, carind, validfind)
	if ind_right == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_RIGHT, pdset, carind, validfind) # distance between cars
	dv = get(V_X_RIGHT, pdset, carind, validfind) # v_other - v_me

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	min(dv*dv / (2dx), THRESHOLD_A_REQ)
end
function get( ::Feature_A_REQ_RIGHT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	
	ind_right = get(INDRIGHT, pdset, sn, carind, validfind)
	if ind_right == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_RIGHT, pdset, sn, carind, validfind) # distance between cars
	dv = get(V_X_RIGHT, pdset, sn, carind, validfind) # v_other - v_me

	if (dx > 0.0 && dv > 0.0) || (dx < 0.0 && dv < 0.0)
		return NA_ALIAS
	end

	if dx > 0.0
		-min(dv*dv / (2*dx), THRESHOLD_A_REQ)
	else
		min(dv*dv / (2*abs(dx)), THRESHOLD_A_REQ)
	end
end

create_feature_basics( "TTC_X_RIGHT", MarkovFeature, "s", false, false, Inf, 0.0, true, :ttc_x_right, L"ttc_{x,ri}", "time to collision with right car assuming constant velocities")
function get( ::Feature_TTC_X_RIGHT, scene::TrafficScene, carind::Int)

	# assumes both cars maintain constant vel
	# note: this assumes cars have zero-length, which is NOT the case

	ind_right = get(INDRIGHT, scene, carind)
	if ind_right == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_RIGHT, scene, carind) # distance between cars
	dv = get(V_X_RIGHT, scene, carind) # v_other - v_me (- if I am going to run into them)

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	min(-dx / dv, THRESHOLD_TIME_TO_COLLISION)
end
function get( ::Feature_TTC_X_RIGHT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	ind_right = get(INDRIGHT, pdset, carind, validfind)
	if ind_right == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_RIGHT, pdset, carind, validfind) # distance between cars
	dv = get(V_X_RIGHT, pdset, carind, validfind) # v_other - v_me

	if dv >= 0.0 # they are pulling away; we are good
		return NA_ALIAS
	end

	min(-dx / dv, THRESHOLD_TIME_TO_COLLISION)
end
function get( ::Feature_TTC_X_RIGHT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	
	ind_right = get(INDRIGHT, pdset, sn, carind, validfind)
	if ind_right == NA_ALIAS
		return NA_ALIAS
	end

	dx = get(D_X_RIGHT, pdset, sn, carind, validfind) # distance between cars
	dv = get(V_X_RIGHT, pdset, sn, carind, validfind) # v_other - v_me

	if (dx > 0.0 && dv > 0.0) || (dx < 0.0 && dv < 0.0)
		return NA_ALIAS
	end

	min(abs(dx / dv), THRESHOLD_TIME_TO_COLLISION)
end

create_feature_basics( "Timegap_X_RIGHT", MarkovFeature, "s", false, false, Inf, 0.0, true, :timegap_x_right, L"\tau_{x,ri}", "timegap with right car")
function get( ::Feature_Timegap_X_RIGHT, scene::TrafficScene, carind::Int)

	ind_right = get(INDRIGHT, scene, carind)
	if ind_right == NA_ALIAS
		return THRESHOLD_TIMEGAP
	end

	dx = get(D_X_RIGHT, scene, carind) # distance between cars
	v  = get(VELFX,     scene, carind)

	if v <= 0.0
		return THRESHOLD_TIMEGAP
	end

	dx / v
end
function get( ::Feature_Timegap_X_RIGHT, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	ind_right = get(INDRIGHT, pdset, carind, validfind)
	if ind_right == NA_ALIAS
		return THRESHOLD_TIMEGAP
	end

	dx = get(D_X_RIGHT, pdset, carind, validfind) # distance between cars
	v  = get(VELFX,     pdset, carind, validfind)

	if v <= 0.0
		return THRESHOLD_TIMEGAP
	end

	dx / v
end
function get( ::Feature_Timegap_X_RIGHT, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	
	ind_right = get(INDRIGHT, pdset, sn, carind, validfind)
	if ind_right == NA_ALIAS
		return THRESHOLD_TIMEGAP
	end

	dx = get(D_X_RIGHT, pdset, sn, carind, validfind) # distance between cars
	v  = get(VELFX,     pdset, sn, carind, validfind)

	if (dx > 0.0 && v > 0.0) || (dx < 0.0 && v < 0.0)
		return THRESHOLD_TIMEGAP
	end

	min(abs(dx / v), THRESHOLD_TIMEGAP)
end

create_feature_basics( "SceneVelFx", MarkovFeature, "m", false, false, Inf, -Inf, false, :scene_velFx, L"\|v\|_{scene}", "average velocity along the lane across all cars in the scene")
function _get( ::Feature_SceneVelFx, scene::TrafficScene, carind::Int)
	val = 0.0
	count = 0
	for car in scene.cars
		val += car.spd * cos(car.pos.yaw)
		count += 1
	end
	if count == 0
		return NA_ALIAS
	end
	val / count
end
function _get( ::Feature_SceneVelFx, pdset::PrimaryDataset, carind::Int, validfind::Int)
	count = 1
	val   = gete(pdset, :velFx, validfind2frameind(pdset, validfind))
	for carind = 0 : ncars_in_frame(pdset, validfind)-1
		val += getc(pdset, "velFx", carind, validfind)
		count += 1
	end
	val / count
end
		 _get( ::Feature_SceneVelFx, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int) =
		 	_get(SCENEVELFX, pdset, carind, validfind)

# ----------------------------------
# iterated features

create_feature_basics( "TurnRate", UnextractableFeature, "rad/s", false, false, Inf, -Inf, false, :turnrate, L"\dot{\psi}", "turn rate")
function _get(::Feature_TurnRate, pdset::PrimaryDataset, carind::Int, validfind::Int)
	# NOTE(tim): defaults to 0.0

	past_validfind = jumpframe(pdset, validfind, -1) # look back one frame
	if past_validfind == 0 # Does not exist
		return 0.0
	end

	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		curr = gete(pdset, :posFyaw, frameind)
		past = gete(pdset, :posFyaw, frameind-1)
		return deltaangle(curr, past) / SEC_PER_FRAME
	end

	carid = carind2id(pdset, carind, validfind)
	if idinframe(pdset, carid, past_validfind)
		farind = carid2ind(pdset, carid, past_validfind)
		curr = getc(pdset, "posFyaw", carind, validfind)
		past = getc(pdset, "posFyaw", farind, past_validfind)
		return deltaangle(curr, past) / SEC_PER_FRAME
	end
	0.0
end

create_feature_basics( "TurnRate_Global", UnextractableFeature, "rad/s", false, false, Inf, -Inf, false, :turnrate_global, L"\dot{\psi}^G", "turn rate in the global frame")
function _get(::Feature_TurnRate_Global, pdset::PrimaryDataset, carind::Int, validfind::Int)
	# NOTE(tim): defaults to 0.0

	past_validfind = jumpframe(pdset, validfind, -1) # look back one frame
	if past_validfind == 0 # Does not exist
		return 0.0
	end

	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		curr = gete(pdset, :posGyaw, frameind)
		past = gete(pdset, :posGyaw, frameind-1)
		return deltaangle(curr, past) / SEC_PER_FRAME
	end

	carid = carind2id(pdset, carind, validfind)
	if idinframe(pdset, carid, past_validfind)
		farind = carid2ind(pdset, carid, past_validfind)
		curr = getc(pdset, "posGyaw", carind, validfind)
		past = getc(pdset, "posGyaw", farind, past_validfind)
		return deltaangle(curr, past) / SEC_PER_FRAME
	end
	0.0
end

create_feature_basics( "AccFx", UnextractableFeature, "m/s2", false, false, Inf, -Inf, false, :accFx, L"a^F_x", "instantaneous acceleration along the lane")
function _get(::Feature_AccFx, pdset::PrimaryDataset, carind::Int, validfind::Int)
	# NOTE(tim): defaults to 0.0

	past_validfind = jumpframe(pdset, validfind, -1) # look back one frame
	if past_validfind == 0 # Does not exist
		return 0.0
	end

	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		curr = gete(pdset, :velFx, frameind)
		past = gete(pdset, :velFx, frameind-1)
		return (curr - past) / SEC_PER_FRAME
	end

	carid = carind2id(pdset, carind, validfind)
	if idinframe(pdset, carid, past_validfind)
		farind = carid2ind(pdset, carid, past_validfind)
		curr = getc(pdset, "velFx", carind, validfind)
		past = getc(pdset, "velFx", farind, past_validfind)
		return (curr - past) / SEC_PER_FRAME
	end
	0.0
end

create_feature_basics( "AccFy", UnextractableFeature, "m/s2", false, false, Inf, -Inf, false, :accFy, L"a^F_y", "instantaneous acceleration perpendicular to the lane")
function _get(::Feature_AccFy, pdset::PrimaryDataset, carind::Int, validfind::Int)
	# NOTE(tim): defaults to 0.0

	past_validfind = jumpframe(pdset, validfind, -1) # look back one frame
	if past_validfind == 0 # Does not exist
		return 0.0
	end

	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		curr = gete(pdset, :velFy, frameind)
		past = gete(pdset, :velFy, frameind-1)
		return (curr - past) / SEC_PER_FRAME
	end

	carid = carind2id(pdset, carind, validfind)
	if idinframe(pdset, carid, past_validfind)
		farind = carid2ind(pdset, carid, past_validfind)
		curr = getc(pdset, "velFy", carind, validfind)
		past = getc(pdset, "velFy", farind, past_validfind)
		return (curr - past) / SEC_PER_FRAME
	end
	0.0
end

create_feature_basics( "Acc", UnextractableFeature, "m/s2", false, false, Inf, -Inf, false, :acc, L"a", "instantaneous ego longitudinal acceleration")
function _get(::Feature_Acc, pdset::PrimaryDataset, carind::Int, validfind::Int)
	# NOTE(tim): defaults to 0.0

	past_validfind = int(jumpframe(pdset, validfind, -1)) # look back one frame
	if past_validfind == 0 # Does not exist
		return 0.0
	end

	if carind == CARIND_EGO
		curr = get(SPEED, pdset, carind, validfind)
		past = get(SPEED, pdset, carind, past_validfind)
		return (curr - past) / SEC_PER_FRAME
	end

	carid = carind2id(pdset, carind, validfind)
	if idinframe(pdset, carid, past_validfind)
		farind = int(carid2ind(pdset, carid, past_validfind))
		curr = get(SPEED, pdset, carind, validfind)
		past = get(SPEED, pdset, farind, past_validfind)
		return (curr - past) / SEC_PER_FRAME
	end
	0.0
end

# ----------------------------------
# unextractable features

create_feature_basics( "D_OnRamp", UnextractableFeature, "m", false, false, Inf, 0.0, true, :d_onramp, L"d_{onramp}", "the distance along the RHS until the next onramp")
function  get(::Feature_D_OnRamp, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		d_onramp = gete(pdset, :d_onramp, frameind)
		return isa(d_onramp, NAtype) ? NA_ALIAS : d_onramp
	end
	d_onramp = getc(pdset, "d_onramp", carind, validfind)
	return isa(d_onramp, NAtype) ? NA_ALIAS : d_onramp
end
get(::Feature_D_OnRamp, ::PrimaryDataset, ::StreetNetwork, ::Int, ::Int) = error("d_onramp unsupported for StreetNetwork use d_merge instead")

create_feature_basics( "D_OffRamp", UnextractableFeature, "m", false, false, Inf, 0.0, true, :d_offramp, L"d_{offramp}", "the distance along the RHS until the next offramp")
function  get(::Feature_D_OffRamp, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		d_offramp = gete(pdset, :d_offramp, frameind)
		return isa(d_offramp, NAtype) ? NA_ALIAS : d_offramp
	end
	d_offramp = getc(pdset, "d_offramp", carind, validfind)
	return isa(d_offramp, NAtype) ? NA_ALIAS : d_offramp
end
get(::Feature_D_OffRamp, ::PrimaryDataset, ::StreetNetwork, ::Int, ::Int) = error("d_onramp unsupported for StreetNetwork use d_split instead")

create_feature_basics( "D_Merge", UnextractableFeature, "m", false, false, Inf, 0.0, true, :d_merge, L"d_{merge}", "the distance along the lane until it merges")
get(::Feature_D_Merge, ::PrimaryDataset, ::Int, ::Int) = error("d_merge unsupported without StreetNetwork; use d_offramp instead")
function  get(::Feature_D_Merge, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		d_merge = gete(pdset, :d_merge, frameind)
		return isa(d_merge, NAtype) ? NA_ALIAS : d_merge
	end
	d_merge = getc(pdset, "d_merge", carind, validfind)
	return isa(d_merge, NAtype) ? NA_ALIAS : d_merge
end

create_feature_basics( "D_Split", UnextractableFeature, "m", false, false, Inf, 0.0, true, :d_split, L"d_{split}", "the distance along the lane until it splits")
get(::Feature_D_Split, ::PrimaryDataset, ::Int, ::Int) = error("d_split unsupported without StreetNetwork; use d_onramp instead")
function  get(::Feature_D_Split, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		d_split = gete(pdset, :d_split, frameind)
		return isa(d_split, NAtype) ? NA_ALIAS : d_split
	end
	d_split = getc(pdset, "d_split", carind, validfind)
	return isa(d_onramp, NAtype) ? NA_ALIAS : d_split
end

create_feature_basics( "FutureTurnRate_250ms", UnextractableFeature, "rad/s", false, false, Inf, -Inf, true, :f_turnrate_250ms, L"\dot{\psi}^{\text{fut}}_{250ms}", "the average rate of heading change over the next quarter second")
function _get(::Feature_FutureTurnRate_250ms, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	const lookahead = 5
	timestep = lookahead * SEC_PER_FRAME

	futrvfind = jumpframe(pdset, validfind, lookahead)
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		cur = gete(pdset, :posFyaw, frameind)
		fut = gete(pdset, :posFyaw, frameind+lookahead)
		return deltaangle(fut, cur) / timestep
	end

	carid = carind2id(pdset, carind, validfind)
	if idinframe(pdset, carid, futrvfind)
		farind = carid2ind(pdset, carid, futrvfind)
		cur = getc(pdset, "posFyaw", carind, validfind)
		fut = getc(pdset, "posFyaw", farind, futrvfind)
		return deltaangle(fut, cur) / timestep
	end
	NA_ALIAS
end

create_feature_basics( "FutureTurnRate_500ms", UnextractableFeature, "rad/s", false, false, Inf, -Inf, true, :f_turnrate_500ms, L"\dot{\psi}^{\text{fut}}_{500ms}", "the average rate of heading change over the next half second")
function _get(::Feature_FutureTurnRate_500ms, pdset::PrimaryDataset, carind::Int, validfind::Int)
	
	const lookahead = 10
	timestep = lookahead * SEC_PER_FRAME

	futrvfind = jumpframe(pdset, validfind, lookahead)
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		cur = gete(pdset, :posFyaw, frameind)
		fut = gete(pdset, :posFyaw, frameind+lookahead)
		return deltaangle(fut, cur) / timestep
	end

	carid = carind2id(pdset, carind, validfind)
	if idinframe(pdset, carid, futrvfind)
		farind = carid2ind(pdset, carid, futrvfind)
		cur = getc(pdset, "posFyaw", carind, validfind)
		fut = getc(pdset, "posFyaw", farind, futrvfind)
		return deltaangle(fut, cur) / timestep
	end
	NA_ALIAS
end

create_feature_basics( "FutureAcceleration_250ms", UnextractableFeature, "m/s2", false, false, Inf, -Inf, true, :f_accel_250ms, L"a^{\text{fut}}_{250ms}", "the average rate of speed change over the next quarter second")
function _get(::Feature_FutureAcceleration_250ms, pdset::PrimaryDataset, carind::Int, validfind::Int)

	const lookahead = 5

	futrvfind = jumpframe(pdset, validfind, lookahead)
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	timestep  = lookahead * SEC_PER_FRAME
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		curx = gete(pdset, :velFx, frameind)
		cury = gete(pdset, :velFy, frameind)
		cur  = hypot(curx, cury)
		futx = gete(pdset, :velFx, frameind + lookahead)
		futy = gete(pdset, :velFy, frameind + lookahead)
		fut  = hypot(futx, futy)
		return (fut - cur) / timestep
	end

	carid = carind2id(pdset, carind, validfind)
	if idinframe(pdset, carid, futrvfind)
		farind = carid2ind(pdset, carid, futrvfind)
		curx = getc(pdset, "velFx", carind, validfind)
		cury = getc(pdset, "velFy", carind, validfind)
		cur  = hypot(curx, cury)
		futx = getc(pdset, "velFx", farind, futrvfind)
		futy = getc(pdset, "velFy", farind, futrvfind)
		fut  = hypot(futx, futy)
		return (fut - cur) / timestep
	end
	NA_ALIAS
end

create_feature_basics( "FutureAcceleration_500ms", UnextractableFeature, "m/s2", false, false, Inf, -Inf, true, :f_accel_500ms, L"a^{\text{fut}}_{500ms}", "the average rate of speed change over the next half second")
function _get(::Feature_FutureAcceleration_500ms, pdset::PrimaryDataset, carind::Int, validfind::Int)

	const lookahead = 10

	futrvfind = jumpframe(pdset, validfind, lookahead)
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	timestep  = lookahead * SEC_PER_FRAME
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		curx = gete(pdset, :velFx, frameind)
		cury = gete(pdset, :velFy, frameind)
		cur  = hypot(curx, cury)
		futx = gete(pdset, :velFx, frameind + lookahead)
		futy = gete(pdset, :velFy, frameind + lookahead)
		fut  = hypot(futx, futy)
		return (fut - cur) / timestep
	end

	carid = carind2id(pdset, carind, validfind)
	if idinframe(pdset, carid, futrvfind)
		farind = carid2ind(pdset, carid, futrvfind)
		curx = getc(pdset, "velFx", carind, validfind)
		cury = getc(pdset, "velFy", carind, validfind)
		cur  = hypot(curx, cury)
		futx = getc(pdset, "velFx", farind, futrvfind)
		futy = getc(pdset, "velFy", farind, futrvfind)
		fut  = hypot(futx, futy)
		return (fut - cur) / timestep
	end
	NA_ALIAS
end

create_feature_basics( "FutureDesiredAngle_250ms", UnextractableFeature, "rad", false, false, float64(π), -float64(π), true, :f_des_angle_250ms, L"\phi^{\text{des}}_{250ms}", "the inferred desired heading angle over the next quarter second")
function _get(::Feature_FutureDesiredAngle_250ms, pdset::PrimaryDataset, carind::Int, validfind::Int)

	const lookahead = 5

	futrvfind = jumpframe(pdset, validfind, lookahead)
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	curϕ = futϕ = 0.0
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		curϕ = gete(pdset, :posFyaw, frameind)
		futϕ = gete(pdset, :posFyaw, frameind + lookahead)
	else
		carid = carind2id(pdset, carind, validfind)
		if idinframe(pdset, carid, futrvfind)
			farind = carid2ind(pdset, carid, futrvfind)
			curϕ = getc(pdset, "posFyaw", carid, validfind)
			futϕ = getc(pdset, "posFyaw", farind, futrvfind)
		else
			return NA_ALIAS
		end
	end

	T  = lookahead * SEC_PER_FRAME
	e = exp(-KP_DESIRED_ANGLE*T)

	(futϕ - curϕ*e) / (1.0 - e)
end

create_feature_basics( "FutureDesiredAngle_500ms", UnextractableFeature, "rad", false, false, float64(π), -float64(π), true, :f_des_angle_500ms, L"\phi^{\text{des}}_{500ms}", "the inferred desired heading angle over the next half second")
function _get(::Feature_FutureDesiredAngle_500ms, pdset::PrimaryDataset, carind::Int, validfind::Int)

	const lookahead = 10

	futrvfind = jumpframe(pdset, validfind, lookahead)
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	curϕ = futϕ = 0.0
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		curϕ = gete(pdset, :posFyaw, frameind)
		futϕ = gete(pdset, :posFyaw, frameind + lookahead)
	else
		carid = carind2id(pdset, carind, validfind)
		if idinframe(pdset, carid, futrvfind)
			farind = carid2ind(pdset, carid, futrvfind)
			curϕ = getc(pdset, "posFyaw", carid, validfind)
			futϕ = getc(pdset, "posFyaw", farind, futrvfind)
		else
			return NA_ALIAS
		end
	end

	T  = lookahead * SEC_PER_FRAME
	e = exp(-KP_DESIRED_ANGLE*T)

	(futϕ - curϕ*e) / (1.0 - e)
end

function _get_futuredesired_speed(pdset::PrimaryDataset, carind::Int, validfind::Int, lookahead::Int)

	futrvfind = int(jumpframe(pdset, validfind, lookahead))
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	curv = get(SPEED, pdset, carind, validfind)
	futv = get(SPEED, pdset, carind, futrvfind)
	if isinf(futv) || isinf(curv)
		return NA_ALIAS
	end

	T = lookahead * SEC_PER_FRAME
	e = exp(-KP_DESIRED_SPEED*T)

	(futv - curv*e) / (1.0 - e) # - curv
end
create_feature_basics( "FutureDesiredSpeed_250ms", UnextractableFeature, "m/s", false, false, Inf, -Inf, true, :f_des_speed_250ms, L"|v|^{\text{des}}_{250ms}", "the inferred desired speed over the next quarter second")
_get(::Feature_FutureDesiredSpeed_250ms, pdset::PrimaryDataset, carind::Int, validfind::Int) =
	_get_futuredesired_speed(pdset, carind, validfind,  5)
create_feature_basics( "FutureDesiredSpeed_500ms", UnextractableFeature, "m/s", false, false, Inf, -Inf, true, :f_des_speed_500ms, L"|v|^{\text{des}}_{500ms}", "the inferred desired speed over the next half second")
_get(::Feature_FutureDesiredSpeed_500ms, pdset::PrimaryDataset, carind::Int, validfind::Int) =
	_get_futuredesired_speed(pdset, carind, validfind, 10)

function _get_futureaccel_control(pdset::PrimaryDataset, carind::Int, validfind::Int, lookahead::Int)

	futrvfind = int(jumpframe(pdset, validfind, lookahead))
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	curv = get(SPEED, pdset, carind, validfind)
	futv = get(SPEED, pdset, carind, futrvfind)
	if isinf(futv) || isinf(curv)
		return NA_ALIAS
	end

	T = lookahead * SEC_PER_FRAME

	a = (futv - curv) / T
	a / (SPEED_LIMIT - curv) # Kp = a / (v_des - v)
end
create_feature_basics( "FutureAccelControl_250ms", UnextractableFeature, "m/s", false, false, Inf, -Inf, true, :f_acc_control_250ms, L"K^a_{250ms}", "the inferred speed control constant over the next quarter second")
_get(::Feature_FutureAccelControl_250ms, pdset::PrimaryDataset, carind::Int, validfind::Int) =
	_get_futureaccel_control(pdset, carind, validfind,  5)
create_feature_basics( "FutureAccelControl_500ms", UnextractableFeature, "m/s", false, false, Inf, -Inf, true, :f_acc_control_500ms, L"K^a_{500ms}", "the inferred speed control constant over the next half second")
_get(::Feature_FutureAccelControl_500ms, pdset::PrimaryDataset, carind::Int, validfind::Int) =
	_get_futureaccel_control(pdset, carind, validfind, 10)

create_feature_basics( "FutureDeltaY2s", UnextractableFeature, "m", false, false, Inf, -Inf, true, :f_deltaY2s, L"{\Delta Y}^{\text{fut}}_{2s}", "the change in frenet frame y coordinate in 2 seconds")
function _get(::Feature_FutureDeltaY2s, pdset::PrimaryDataset, carind::Int, validfind::Int)

	timestep = 2.0
	jump = iround(timestep / SEC_PER_FRAME) # number of frames to look ahead

	futrvfind = jumpframe(pdset, validfind, jump)
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	cur,fut = 0.0,0.0
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		cur = gete(pdset, :posFy, frameind)
		fut = gete(pdset, :posFy, frameind + jump)
	elseif idinframe(pdset, carind2id(pdset, carind, validfind), futrvfind)
		carid = carind2id(pdset, carind, validfind)
		farind = carid2ind(pdset, carid, futrvfind)
		cur = getc(pdset, "posFy", carind, validfind)
		fut = getc(pdset, "posFy", farind, futrvfind)
	else
		return NA_ALIAS
	end

	fut-cur # positive is to the car's left
end
get(::Feature_FutureDeltaY2s, ::PrimaryDataset, ::StreetNetwork, ::Int, ::Int) = error("future delta Y 2s unsupported for StreetNetwork")

create_feature_basics( "FutureDeltaY1s", UnextractableFeature, "m", false, false, Inf, -Inf, true, :f_deltaY1s, L"{\Delta Y}^{\text{fut}}_{1s}", "the change in frenet frame y coordinate in 1 second")
function _get(::Feature_FutureDeltaY1s, pdset::PrimaryDataset, carind::Int, validfind::Int)

	timestep = 1.0
	jump = iround(timestep / SEC_PER_FRAME) # number of frames to look ahead

	futrvfind = jumpframe(pdset, validfind, jump)
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	cur,fut = 0.0,0.0
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		cur = gete(pdset, :posFy, frameind)
		fut = gete(pdset, :posFy, frameind + jump)
	elseif idinframe(pdset, carind2id(pdset, carind, validfind), futrvfind)
		carid = carind2id(pdset, carind, validfind)
		farind = carid2ind(pdset, carid, futrvfind)
		cur = getc(pdset, "posFy", carind, validfind)
		fut = getc(pdset, "posFy", farind, futrvfind)
	else
		return NA_ALIAS
	end

	fut-cur # positive is to the car's left
end
get(::Feature_FutureDeltaY1s, ::PrimaryDataset, ::StreetNetwork, ::Int, ::Int) = error("future delta Y 1s unsupported for StreetNetwork")

create_feature_basics( "FutureDeltaY_250ms", UnextractableFeature, "m", false, false, Inf, -Inf, true, :f_deltaY_250ms, L"{\Delta Y}^{\text{fut}}_{250ms}", "the change in frenet y over the next quarter second")
function _get(::Feature_FutureDeltaY_250ms, pdset::PrimaryDataset, carind::Int, validfind::Int)

	timestep = 0.25
	jump = iround(timestep / SEC_PER_FRAME) # number of frames to look ahead

	futrvfind = jumpframe(pdset, validfind, jump)
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	cur,fut = 0.0,0.0
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		cur = gete(pdset, :posFy, frameind)
		fut = gete(pdset, :posFy, frameind + jump)
	elseif idinframe(pdset, carind2id(pdset, carind, validfind), futrvfind)
		carid = carind2id(pdset, carind, validfind)
		farind = carid2ind(pdset, carid, futrvfind)
		cur = getc(pdset, "posFy", carind, validfind)
		fut = getc(pdset, "posFy", farind, futrvfind)
	else
		return NA_ALIAS
	end

	fut-cur # positive is to the car's left
end
get(::Feature_FutureDeltaY1s, ::PrimaryDataset, ::StreetNetwork, ::Int, ::Int) = error("future delta Y 250ms unsupported for StreetNetwork")

create_feature_basics( "ID", UnextractableFeature, "-", true, false, Inf, -1.0, false, :id, L"\text{id}", "the corresponding carid, -1 for the ego car")
get(::Feature_ID, pdset::PrimaryDataset, carind::Int, validfind::Int) = float64(carind2id(pdset, carind, validfind))

create_feature_basics( "LaneCurvature", UnextractableFeature, "1/m", false, false, Inf, -Inf, false, :curvature, L"\kappa", "the local lane curvature")
function get(::Feature_LaneCurvature, pdset::PrimaryDataset, carind::Int, validfind::Int)
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		return gete(pdset, :curvature, frameind)
	end
	getc(pdset, "curvature", carind, validfind)
end
	     get(::Feature_LaneCurvature, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int) = get(LANECURVATURE, pdset, carind, validfind)

create_feature_basics( "TimeToLaneCrossing", UnextractableFeature, "s", false, false, THRESHOLD_TIMETOLANECROSSING, 0.0, false, :timetolanecrossing, L"t_{\text{crossing}}^{+}", "the time until the next lane crossing")
function _get(::Feature_TimeToLaneCrossing, pdset::PrimaryDataset, carind::Int, validfind::Int)

	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		cur = gete(pdset, :lane, frameind)
		jump = 0
		fut = cur
		while fut == cur && frameind_inbounds(pdset, frameind+jump+1)
			jump += 1
			fut = gete(pdset, :lane, frameind+jump)
			if isa(fut, NAtype)
				fut = cur
			end
		end
		return fut != cur ? jump*SEC_PER_FRAME : NA_ALIAS
	else
		carid = carind2id(pdset, carind, validfind)
		cur = getc(pdset, "lane", carind, validfind)
		jump = 0
		fut = cur
		t_inview = getc(pdset, "t_inview", carind, validfind) # time in view
		while fut == cur && jump*SEC_PER_FRAME < t_inview
			jump += 1
			futrvfind = jumpframe(pdset, validfind, jump)
			if futrvfind != 0 && idinframe(pdset, carid, futrvfind)
				farind = carid2ind(pdset, carid, futrvfind)
				fut = getc(pdset, "lane", farind, futrvfind)
			end
		end
		return fut != cur ? jump*SEC_PER_FRAME : NA_ALIAS
	end
end
function _get(::Feature_TimeToLaneCrossing, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	# scan forward until the car is no longer in the same lane

	frameind = validfind2frameind(pdset, validfind)
	t_orig = gete(pdset, :time, frameind)::Float64

	cur_lanetag = get(pdset, :lanetag, "lanetag", carind, frameind, validfind)::LaneTag
	cur_lane = get_lane(sn, cur_lanetag)
	cur_frameind = frameind

	done = false
	while !done
		cur_frameind += 1
		cur_validfind = frameind2validfind(pdset, cur_frameind)
		if cur_validfind != 0
			Δt = gete(pdset, :time, cur_frameind) - t_orig
			if Δt > THRESHOLD_TIMETOLANECROSSING
				return THRESHOLD_TIMETOLANECROSSING
			end
			fut_lanetag = get(pdset, :lanetag, "lanetag", carind, cur_frameind, cur_validfind)::LaneTag
			if fut_lanetag != cur_lanetag					
				if same_tile(cur_lanetag, fut_lanetag) || !has_next_lane(sn, cur_lane)
					return Δt
				else
					cur_lane = next_lane(sn, cur_lane)
					cur_lanetag = LaneTag(sn, cur_lane)
					if fut_lanetag != cur_lanetag
						return Δt
					end
				end
			end
		else
			return THRESHOLD_TIMETOLANECROSSING
		end
	end

	error("INVALID CODEPATH")
	return THRESHOLD_TIMETOLANECROSSING
end

create_feature_basics( "TimeSinceLaneCrossing", UnextractableFeature, "s", false, false, THRESHOLD_TIMESINCELANECROSSING, 0.0, false, :timesincelanecrossing, L"t_{\text{crossing}}^{-}", "the time since the last lane crossing")
function _get(::Feature_TimeSinceLaneCrossing, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	# scan backward until the car is no longer in the same lane
	# NOTE(tim): returns positive time values

	frameind = validfind2frameind(pdset, validfind)
	t_orig = gete(pdset, :time, frameind)::Float64

	cur_lanetag = get(pdset, :lanetag, "lanetag", carind, frameind, validfind)::LaneTag
	cur_lane = get_lane(sn, cur_lanetag)
	cur_frameind = frameind

	done = false
	while !done
		cur_frameind -= 1
		cur_validfind = frameind2validfind(pdset, cur_frameind)
		if cur_validfind != 0
			Δt = t_orig - gete(pdset, :time, cur_frameind)
			if Δt > THRESHOLD_TIMESINCELANECROSSING
				return THRESHOLD_TIMESINCELANECROSSING
			end
			past_lanetag = get(pdset, :lanetag, "lanetag", carind, cur_frameind, cur_validfind)::LaneTag
			if past_lanetag != cur_lanetag					
				if same_tile(cur_lanetag, past_lanetag) || !has_prev_lane(sn, cur_lane)
					return Δt
				else
					cur_lane = prev_lane(sn, cur_lane)
					cur_lanetag = LaneTag(sn, cur_lane)
					if past_lanetag != cur_lanetag
						return Δt
					end
				end
			end
		else
			return THRESHOLD_TIMESINCELANECROSSING
		end
	end

	error("INVALID CODEPATH")
	return THRESHOLD_TIMESINCELANECROSSING
end

create_feature_basics( "Time_Consecutive_Brake", UnextractableFeature, "s", false, false, THRESHOLD_TIMECONSECUTIVEBRAKE, 0.0, false, :time_consecutive_brake, L"t_\text{brake}", "the consecutive time during which the car has been decelerating")
function _get(::Feature_Time_Consecutive_Brake, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	# scan backward until the car is no longer braking
	# NOTE(tim): returns positive time values

	const THRESHOLD_BRAKING = -0.05 # [m/s²]

	frameind = validfind2frameind(pdset, validfind)
	t_orig = gete(pdset, :time, frameind)::Float64

	past_validfind = int(jumpframe(pdset, validfind, -1))
	if past_validfind == 0
		meta[carind][:time_consecutive_accel] = 0.0
		return 0.0 # default
	end

	cur_accel = get(ACC, pdset, sn, carind, validfind)
	if cur_accel > THRESHOLD_BRAKING
		return 0.0
	end

	meta[carind][:time_consecutive_accel] = 0.0

	past_frameind = frameind
	done = false
	while !done
		past_frameind -= 1
		past_validfind = int(jumpframe(pdset, past_validfind, -1))
		if past_validfind != 0
			past_accel = get(ACC, pdset, sn, carind, past_validfind)
			if past_accel > THRESHOLD_BRAKING
				return t_orig - gete(pdset, :time, past_frameind+1)
			end

			Δt = t_orig - gete(pdset, :time, past_frameind)	
			if Δt > THRESHOLD_TIMECONSECUTIVEBRAKE
				return THRESHOLD_TIMECONSECUTIVEBRAKE
			end
		else
			return t_orig - gete(pdset, :time, past_frameind+1)
		end
	end

	error("INVALID CODEPATH")
	return 0.0
end

create_feature_basics( "Time_Consecutive_Accel", UnextractableFeature, "s", false, false, THRESHOLD_TIMECONSECUTIVEACCEL, 0.0, false, :time_consecutive_accel, L"t_\text{accel}", "the consecutive time during which the car has been accelerating")
function _get(::Feature_Time_Consecutive_Accel, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)
	# scan backward until the car is no longer braking
	# NOTE(tim): returns positive time values

	const THRESHOLD_ACCELERATION = 0.05 # [m/s²]

	frameind = validfind2frameind(pdset, validfind)
	t_orig = gete(pdset, :time, frameind)::Float64

	past_validfind = int(jumpframe(pdset, validfind, -1))
	if past_validfind == 0
		meta[carind][:time_consecutive_brake] = 0.0
		return 0.0 # default
	end

	cur_accel = get(ACC, pdset, sn, carind, validfind)
	if cur_accel < THRESHOLD_ACCELERATION
		return 0.0
	end

	meta[carind][:time_consecutive_brake] = 0.0

	past_frameind = frameind
	done = false
	while !done
		past_frameind -= 1
		past_validfind = int(jumpframe(pdset, past_validfind, -1))
		if past_validfind != 0

			past_accel = get(ACC, pdset, sn, carind, past_validfind)
			if past_accel < THRESHOLD_ACCELERATION
				return t_orig - gete(pdset, :time, past_frameind+1)
			end

			Δt = t_orig - gete(pdset, :time, past_frameind)	
			if Δt > THRESHOLD_TIMECONSECUTIVEBRAKE
				return THRESHOLD_TIMECONSECUTIVEBRAKE
			end
		else
			return t_orig - gete(pdset, :time, past_frameind+1)
		end
	end

	error("INVALID CODEPATH")
	return 0.0
end

create_feature_basics( "LaneChangeDir", UnextractableFeature, "s", true, false, 1.0, -1.0, true, :lanechangedir, L"\delta_{lc}", "the direction of the next lane change. Positive is to the left")
function _get(::Feature_LaneChangeDir, pdset::PrimaryDataset, carind::Int, validfind::Int)

	ttlc = get(TIMETOLANECROSSING, pdset, carind, validfind)
	if ttlc == NA_ALIAS
		return NA_ALIAS
	end

	jump = iround(ttlc / SEC_PER_FRAME)

	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		cur = gete(pdset, :lane, frameind)
		fut = gete(pdset, :lane, frameind + jump)
		return fut > cur ? 1.0 : -1.0
	else
		carid = carind2id(pdset, carind, validfind)
		futrvfind = jumpframe(pdset, validfind, jump)
		farind = carid2ind(pdset, carid, futrvfind)
		cur = getc(pdset, "lane", carind, validfind)
		fut = getc(pdset, "lane", farind, futrvfind)
		return fut > cur ? 1.0 : -1.0
	end
	error("IMPOSSIBLE LOCATION REACHED")
end
get(::Feature_LaneChangeDir, ::PrimaryDataset, ::StreetNetwork, ::Int, ::Int) = error("lane change dir not supported for StreetNetwork")

create_feature_basics( "LaneChange2s", UnextractableFeature, "-", true, true, 2.0, 0.0, true, :lanechange2s, L"\Delta L_{2s}", "whether the car is in another lane two seconds in the future")
function _get(::Feature_LaneChange2s, pdset::PrimaryDataset, carind::Int, validfind::Int)

	timestep = 2.0
	jump = iround(timestep / SEC_PER_FRAME) # number of frames to look ahead

	futrvfind = jumpframe(pdset, validfind, jump)
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	cur,fut = 0.0,0.0
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		cur = gete(pdset, :lane, frameind)
		fut = gete(pdset, :lane, frameind + jump)
	elseif idinframe(pdset, carind2id(pdset, carind, validfind), futrvfind)
		carid = carind2id(pdset, carind, validfind)
		farind = carid2ind(pdset, carid, futrvfind)
		cur = getc(pdset, "lane", carind, validfind)
		fut = getc(pdset, "lane", farind, futrvfind)
	else
		return NA_ALIAS
	end

	islanechange = iround(cur) != iround(fut)
	islanechange ? 1.0 : 0.0
end
get(::Feature_LaneChange2s, ::PrimaryDataset, ::StreetNetwork, ::Int, ::Int) = error("lane change 2s not supported for StreetNetwork")

create_feature_basics( "LaneChange1s", UnextractableFeature, "-", true, true, 2.0, 0.0, true, :lanechange1s, L"\Delta L_{1s}", "whether the car is in another lane one second in the future")
function _get(::Feature_LaneChange1s, pdset::PrimaryDataset, carind::Int, validfind::Int)

	timestep = 1.0
	jump = iround(timestep / SEC_PER_FRAME) # number of frames to look ahead

	futrvfind = jumpframe(pdset, validfind, jump)
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	cur,fut = 0.0,0.0
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		cur = gete(pdset, :lane, frameind)
		fut = gete(pdset, :lane, frameind + jump)
	elseif idinframe(pdset, carind2id(pdset, carind, validfind), futrvfind)
		carid = carind2id(pdset, carind, validfind)
		farind = carid2ind(pdset, carid, futrvfind)
		cur = getc(pdset, "lane", carind, validfind)
		fut = getc(pdset, "lane", farind, futrvfind)
	else
		return NA_ALIAS
	end

	islanechange = iround(cur) != iround(fut)
	islanechange ? 1.0 : 0.0
end
get(::Feature_LaneChange1s, ::PrimaryDataset, ::StreetNetwork, ::Int, ::Int) = error("lane change 1s not supported for StreetNetwork")

create_feature_basics( "LaneChange500ms", UnextractableFeature, "-", true, true, 1.0, 0.0, true, :lanechange500ms, L"\Delta L_{0.5s}", "whether the car is in another lane a half second in the future")
function _get(::Feature_LaneChange500ms, pdset::PrimaryDataset, carind::Int, validfind::Int)

	timestep = 0.5
	jump = iround(timestep / SEC_PER_FRAME) # number of frames to look ahead

	futrvfind = jumpframe(pdset, validfind, jump)
	if futrvfind == 0 # Does not exist
		return NA_ALIAS
	end

	cur,fut = 0.0,0.0
	if carind == CARIND_EGO
		frameind = validfind2frameind(pdset, validfind)
		cur = gete(pdset, :lane, frameind)
		fut = gete(pdset, :lane, frameind + jump)
	elseif idinframe(pdset, carind2id(pdset, carind, validfind), futrvfind)
		carid = carind2id(pdset, carind, validfind)
		farind = carid2ind(pdset, carid, futrvfind)
		cur = getc(pdset, "lane", carind, validfind)
		fut = getc(pdset, "lane", farind, futrvfind)
	else
		return NA_ALIAS
	end

	islanechange = iround(cur) != iround(fut)
	islanechange ? 1.0 : 0.0
end
get(::Feature_LaneChange500ms, ::PrimaryDataset, ::StreetNetwork, ::Int, ::Int) = error("lane change 500ms not supported for StreetNetwork")

# ----------------------------------
# submodels

create_feature_basics( "Subset_Emergency", UnextractableFeature, "-", true, true, 1.0, 0.0, false, :subset_emergency, L"\mathcal{D}_\text{emerg}", "subset of data for emergency behavior")
function _get(::Feature_Subset_Emergency, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int;
	threshold_acc :: Float64 = 2.0, 
	threshold_turnrate :: Float64 = 0.05
	)

	# emergency
	# - events close to a collision
	# - can be found when the driver does a big accel?
	# - definition
	#   - |acc| > threshold
	#           or
	#   - |turnrate| > threshold

	a = get(ACC,      pdset, sn, carind, validfind)
	ω = get(TURNRATE, pdset, sn, carind, validfind)

	retval = abs(a) > threshold_acc ||
	         abs(ω) < threshold_turnrate

	float64(retval)
end

create_feature_basics( "Subset_Free_Flow", UnextractableFeature, "-", true, true, 1.0, 0.0, false, :subset_free_flow, L"\mathcal{D}_\text{free}", "subset of data for free flow")
function _get(::Feature_Subset_Free_Flow, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int;
	threshold_timegap_front :: Float64 = 3.0, 
	threshold_d_v_front     :: Float64 = 0.5
	)

	ΔT = get(TIMEGAP_X_FRONT, pdset, sn, carind, validfind)
	dv = get(      V_X_FRONT, pdset, sn, carind, validfind)

	retval =     ΔT > threshold_timegap_front ||
	             dv > threshold_d_v_front

	float64(retval)
end

create_feature_basics( "Subset_Car_Following", UnextractableFeature, "-", true, true, 1.0, 0.0, false, :subset_car_following, L"\mathcal{D}_\text{follow}", "subset of data for car following")
function _get(::Feature_Subset_Car_Following, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int;
	threshold_timegap_front :: Float64 = 3.0,
	threshold_d_v_front     :: Float64 = 0.5
	)

	# following
	# - when there is a car in front of you
	# - you typically modulate your speed to match theirs and maintain a certain distance
	#   - distance is proportional to speed?
	#   - maintain a time gap?
	# - interesting:
	#   - what form does the timegap take?
	#   - how rigid is the tracking?
	#   - what sort of delays are there?
	# - definition
	#   - timegap_front < threshold
	#   - d_v_front < threshold

	ΔT = get(TIMEGAP_X_FRONT, pdset, sn, carind, validfind)
	dv = get(      V_X_FRONT, pdset, sn, carind, validfind)

	retval =    !(ΔT > threshold_timegap_front ||
	              dv > threshold_d_v_front)

	float64(retval)
end

create_feature_basics( "Subset_Lane_Crossing", UnextractableFeature, "-", true, true, 1.0, 0.0, false, :subset_lane_crossing, L"\mathcal{D}_\text{crossing}", "subset of data near lane crossing")
function _get(::Feature_Subset_Lane_Crossing, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int;
	max_time_to_lane_change::Float64 = 2.0, # [s]
	max_time_since_lane_change::Float64 = 2.0
	)

	# lane change
	# - when the driver changes lanes
	# - interesting:
	#   - what causes a driver to change lanes?
	# - definition
	#   - a lane crossing has occurred

	time_to_lane_change = get(TIMETOLANECROSSING, pdset, sn, carind, validfind)
	time_since_lane_change = get(TIMESINCELANECROSSING, pdset, sn, carind, validfind)

	retval = time_to_lane_change    < max_time_to_lane_change ||
	         time_since_lane_change < max_time_since_lane_change

	float64(retval)
end

create_feature_basics( "Subset_Sustained_Crossing", UnextractableFeature, "-", true, true, 1.0, 0.0, false, :subset_sustained_crossing, L"\mathcal{D}_\text{sustained}", "subset of data with a sustained lane crossing")
function _get(::Feature_Subset_Sustained_Crossing, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int;
	max_time_to_lane_change     :: Float64 = 2.0, # [s]
	max_time_since_lane_change  :: Float64 = 2.0, # [s]
	min_abs_d_cl_at_extrema     :: Float64 = 1.0, # [m]
	req_time_d_cl_threshold_met :: Float64 = 0.1, # [s] amount of time at extrema of segment that must meet d_cl requirements
	)

	had_lane_crossing = _get(SUBSET_LANE_CROSSING, pdset, sn, carind, validfind,
							max_time_to_lane_change = max_time_to_lane_change,
							max_time_since_lane_change = max_time_since_lane_change)

	if isapprox(had_lane_crossing, 0.0)
		return float64(false)
	end
	
	frameind = validfind2frameind(pdset, validfind)
	t = gete(pdset, :time, frameind)::Float64

	time_to_lane_crossing    = get(TIMETOLANECROSSING,    pdset, sn, carind, validfind)
	time_since_lane_crossing = get(TIMESINCELANECROSSING, pdset, sn, carind, validfind)

	if time_to_lane_crossing < time_since_lane_crossing
		Δt = time_to_lane_crossing
	else
		Δt = -time_since_lane_crossing
	end

	closest_validfind = Trajdata.closest_validfind(pdset, t + Δt)
	@assert(closest_validfind != 0)
	closest_frameind = validfind2frameind(pdset, closest_validfind)
	@assert(closest_frameind != 0)

	find, vind = closest_frameind, closest_validfind
	dir  = convert(typeof(closest_frameind), sign(Δt))
	t    = gete(pdset, :time, closest_frameind)::Float64
	dt   = 0.0
	while dt < req_time_d_cl_threshold_met
		d_cl = get(pdset, :d_cl, "d_cl", carind, find, vind)
		if abs(d_cl) > min_abs_d_cl_at_extrema
			return float64(false)
		end
		
		find += dir
		vind  = frameind2validfind(pdset, find)
		while vind == 0
			find += dir
			vind  = frameind2validfind(pdset, find)
		end

		dt = abs(gete(pdset, :time, find)::Float64 - t)
	end

	float64(true)
end

create_feature_basics( "Subset_At_SixtyFive", UnextractableFeature, "-", true, true, 1.0, 0.0, false, :subset_at_sixtyfive, L"\mathcal{D}_{v=65\text{mph}}", "subset of data where the car is close to 65 mph")
function _get(::Feature_Subset_At_SixtyFive, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int;
	speed_tolerance :: Float64 = 2.235 # [m/s²]
	)

	V = get(SPEED, pdset, sn, carind, validfind)
	float64(abs(V - SPEED_LIMIT) < speed_tolerance)
end

create_feature_basics( "Subset_Auto", UnextractableFeature, "-", true, true, 1.0, 0.0, false, :subset_auto, L"s_\text{auto}", "subset of data where the car is in autonomous mode")
function _get(::Feature_Subset_Auto, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)

	if carind != CARIND_EGO
		return float64(false)
	end

	status = gete_validfind(pdset, :control_status, validfind)::Int
	return float64(status == Trajdata.CONTROL_STATUS_AUTO)
end

# ----------------------------------
# misc

clearmeta!(car::RichVehicle) = empty!(car.meta) 
function clearmeta!(scene::TrafficScene)
	for car in scene.cars
		clearmeta!(car)
	end
end
clearmeta() = empty!(meta)
function checkmeta(carind::Int, validfind::Int, symb::Symbol)
	global meta_validfind
	global meta
	if meta_validfind != validfind
		clearmeta()
		meta_validfind = validfind
	end

	if !haskey(meta, carind)
		meta[carind] = Dict{Symbol,Float64}()
		return false	
	end

	haskey(meta[carind], symb)
end

# ----------------------------------

function ticks_to_time_string( ticks::Int )
	n_secs = ticks//FRAME_PER_SEC
	unit = "s"
	if den(n_secs) != 1
		unit = "ms"
		n_secs *= 1000
	end
	@assert(den(n_secs) == 1)
	(num(n_secs), unit)
end

# ----------------------------------

# signed delta angle
deltaangle( a::Real, b::Real ) = atan2(sin(a-b), cos(a-b))
# distance between two angles
angledist( a::Real, b::Real ) = abs(tan2(sin(a-b), cos(a-b)))

function is_pt_left_of_ray(x::Float64, y::Float64, raySx::Float64, raySy::Float64, rayEx::Float64, rayEy::Float64)

  	(y-raySy)*(rayEx-raySx) > (x-raySx)*(rayEy-raySy)
end
function is_pt_left_of_ray(x::Float64, y::Float64, rX::Float64, rY::Float64, rθ::Float64)

	is_pt_left_of_ray(x, y, rX, rY, rX + cos(θ), rY + sin(θ))
end

# ----------------------------------

function replace_featuresymbols_in_string_with_latex(str::String)

	strs = map(k->string(k), keys(sym2ftr))
	lstr = map(f->lsymbol(f), values(sym2ftr))

	p = sortperm(strs, by=length, rev=true)

	for i in p
		str = replace(str, strs[i], lstr[i])
	end
	str
end
function print_feature_table()
	for f in values(sym2ftr)
		println(string(f), " & ", lsymbol(f), " & ", units(f), " & ", description(f))
	end
end

# ----------------------------------

tick_list = [5,10,15,20,30,40,50,60,80] # histories used for accumulated features
tick_list_short = [5,10,15,20]          # histories used for past features

include("features/template_past_accel.jl")
include("features/template_past_turnrate.jl")
include("features/template_past_velFy.jl")
include("features/template_past_d_cl.jl")

include("features/template_max_accel.jl")
include("features/template_max_turnrate.jl")
include("features/template_mean_accel.jl")
include("features/template_mean_turnrate.jl")
include("features/template_std_accel.jl")
include("features/template_std_turnrate.jl")

end # end module