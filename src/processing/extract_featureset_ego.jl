push!(LOAD_PATH, "location_of_base_files");

using DataFrames
using Trajdata
using HDF5, JLD
using StreetMap
using Curves
using Features
using FilesystemUtils

function gen_featureset_ego{F <: AbstractFeature, G <: AbstractFeature}(
    pdset :: PrimaryDataset,
    validfind_regions :: Vector{Int},
    sn    :: StreetNetwork,
    features :: Vector{F};
    filters :: Vector{G} = AbstractFeature[] # list of boolean features; if true the car is kept
    )
    
    maxcarind = get_maxcarind(pdset)

    estimated_row_count = 0
    region_index_lo = 1
    while region_index_lo < length(validfind_regions)
        @assert(validfind_regions[region_index_lo+1] > validfind_regions[region_index_lo])
        estimated_row_count += validfind_regions[region_index_lo+1] - validfind_regions[region_index_lo] + 2
        region_index_lo += 2
    end

    df = DataFrame()
    for f in features
        df[symbol(f)] = DataArray(Float64, estimated_row_count)
    end

    row_index = 0
    region_index_lo = 1
    while region_index_lo < length(validfind_regions)
        for validfind in validfind_regions[region_index_lo] : validfind_regions[region_index_lo+1]

            violates_filter = false
            for f in filters
                if isapprox(get(f, pdset, sn, CARIND_EGO, validfind), 0.0)
                    violates_filter = true
                    break
                end
            end

            if !violates_filter
                row_index += 1
                for f in features
                    value::Float64 = get(f, pdset, sn, CARIND_EGO, validfind)
                    if isa(f, Features.Feature_N_LANE_L) && (value < 0.0 || value > 5.0)
                        println(f, "  ", value)
                    end
                    df[row_index, symbol(f)] = value
                end
            end
        end

        region_index_lo += 2
    end

    # NOTE(tim): drop any extra rows
    df[1:row_index,:]
end

type CSVFileSet
    csvfile :: String
    streetmapbasename :: String

    lanechanges_normal    :: Vector{Int} # lane changes made to pass another vehicle
    lanechanges_postpass  :: Vector{Int} # lane changes made to get back into a lane after passing
    lanechanges_arbitrary :: Vector{Int} # discretionary lane changes
    carfollow             :: Vector{Int}
    freeflow              :: Vector{Int}
end

const PRIMARYDATA_DIR = "location_to_import_primarydata"
const FEATUREMATRIX_DIR = "location_to_output_featurematrices"
const STREETMAP_BASE = "location_to_input_streetmaps"
const CSVFILESETS = (
            CSVFileSet("example", "some_streetmap", [128,297], [298,300], Int[], Int[], Int[])
            #=
            ACTUAL CSVFileSets CENSORED
            =#
        )
const FEATURES = [
        YAW, POSFX, POSFY, SPEED, VELFX, VELFY, DELTA_SPEED_LIMIT,
        D_CL, D_ML, D_MR, D_MERGE, D_SPLIT, 
        TIMETOCROSSING_RIGHT, TIMETOCROSSING_LEFT,
        N_LANE_L, N_LANE_R, HAS_LANE_L, HAS_LANE_R,
        TURNRATE, TURNRATE_GLOBAL, ACC, ACCFX, ACCFY, A_REQ_STAYINLANE, LANECURVATURE,

        HAS_FRONT, D_X_FRONT, D_Y_FRONT, V_X_FRONT, V_Y_FRONT, YAW_FRONT, TURNRATE_FRONT,
        HAS_REAR,  D_X_REAR,  D_Y_REAR,  V_X_REAR,  V_Y_REAR,  YAW_REAR,  TURNRATE_REAR,
                   D_X_LEFT,  D_Y_LEFT,  V_X_LEFT,  V_Y_LEFT,  YAW_LEFT,  TURNRATE_LEFT,
                   D_X_RIGHT, D_Y_RIGHT, V_X_RIGHT, V_Y_RIGHT, YAW_RIGHT, TURNRATE_RIGHT,
        A_REQ_FRONT, TTC_X_FRONT, TIMEGAP_X_FRONT, GAINING_ON_FRONT,
        A_REQ_REAR,  TTC_X_REAR,  TIMEGAP_X_REAR,
        A_REQ_LEFT,  TTC_X_LEFT,  TIMEGAP_X_LEFT,
        A_REQ_RIGHT, TTC_X_RIGHT, TIMEGAP_X_RIGHT,

        SCENEVELFX,

        TIME_CONSECUTIVE_BRAKE, TIME_CONSECUTIVE_ACCEL,
             PASTACC250MS,      PASTACC500MS,      PASTACC750MS,      PASTACC1S,
        PASTTURNRATE250MS, PASTTURNRATE500MS, PASTTURNRATE750MS, PASTTURNRATE1S,
           PASTVELFY250MS,    PASTVELFY500MS,    PASTVELFY750MS,    PASTVELFY1S,
            PASTD_CL250MS,     PASTD_CL500MS,     PASTD_CL750MS,     PASTD_CL1S,

            MAXACCFX250MS,     MAXACCFX500MS,     MAXACCFX750MS,     MAXACCFX1S,     MAXACCFX1500MS,     MAXACCFX2S,     MAXACCFX2500MS,     MAXACCFX3S,     MAXACCFX4S,
            MAXACCFY250MS,     MAXACCFY500MS,     MAXACCFY750MS,     MAXACCFY1S,     MAXACCFY1500MS,     MAXACCFY2S,     MAXACCFY2500MS,     MAXACCFY3S,     MAXACCFY4S,
         MAXTURNRATE250MS,  MAXTURNRATE500MS,  MAXTURNRATE750MS,  MAXTURNRATE1S,  MAXTURNRATE1500MS,  MAXTURNRATE2S,  MAXTURNRATE2500MS,  MAXTURNRATE3S,  MAXTURNRATE4S,
           MEANACCFX250MS,    MEANACCFX500MS,    MEANACCFX750MS,    MEANACCFX1S,    MEANACCFX1500MS,    MEANACCFX2S,    MEANACCFX2500MS,    MEANACCFX3S,    MEANACCFX4S,
           MEANACCFY250MS,    MEANACCFY500MS,    MEANACCFY750MS,    MEANACCFY1S,    MEANACCFY1500MS,    MEANACCFY2S,    MEANACCFY2500MS,    MEANACCFY3S,    MEANACCFY4S,
        MEANTURNRATE250MS, MEANTURNRATE500MS, MEANTURNRATE750MS, MEANTURNRATE1S, MEANTURNRATE1500MS, MEANTURNRATE2S, MEANTURNRATE2500MS, MEANTURNRATE3S, MEANTURNRATE4S,
            STDACCFX250MS,     STDACCFX500MS,     STDACCFX750MS,     STDACCFX1S,     STDACCFX1500MS,     STDACCFX2S,     STDACCFX2500MS,     STDACCFX3S,     STDACCFX4S,
            STDACCFY250MS,     STDACCFY500MS,     STDACCFY750MS,     STDACCFY1S,     STDACCFY1500MS,     STDACCFY2S,     STDACCFY2500MS,     STDACCFY3S,     STDACCFY4S,
         STDTURNRATE250MS,  STDTURNRATE500MS,  STDTURNRATE750MS,  STDTURNRATE1S,  STDTURNRATE1500MS,  STDTURNRATE2S,  STDTURNRATE2500MS,  STDTURNRATE3S,  STDTURNRATE4S,

        FUTUREACCELERATION_250MS, FUTURETURNRATE_250MS, FUTUREDESIREDANGLE_250MS, FUTUREDESIREDSPEED_250MS, FUTUREACCELCONTROL_250MS,
        FUTUREACCELERATION_500MS, FUTURETURNRATE_500MS, FUTUREDESIREDANGLE_500MS, FUTUREDESIREDSPEED_500MS, FUTUREACCELCONTROL_500MS,
        TIMETOLANECROSSING, TIMESINCELANECROSSING,

        SUBSET_EMERGENCY, SUBSET_FREE_FLOW, SUBSET_CAR_FOLLOWING, SUBSET_LANE_CROSSING, SUBSET_SUSTAINED_CROSSING,
        SUBSET_AT_SIXTYFIVE,
    ]
const FILTERS = AbstractFeature[]

const FEATURE_SETS = [
        ("lanechange", csvfileset -> [csvfileset.lanechanges_normal, csvfileset.lanechanges_postpass, csvfileset.lanechanges_arbitrary]),
        ("carfollow", csvfileset -> [csvfileset.carfollow]),
        ("freeflow", csvfileset -> [csvfileset.freeflow]),
    ]

streetnet_cache = Dict{String, StreetNetwork}()

for (featureset_name, f_regions) in FEATURE_SETS

    output_name = "ego_" * featureset_name
    featureset_dir  = joinpath("dir_to_featureset_output", output_name*"/")
    if !isdir(featureset_dir)
        mkdir(featureset_dir)
    end

    aggregate_datafile = DataFrame()

    tic()
    for csvfileset in CSVFILESETS

        csvfile = csvfileset.csvfile
        streetmapbasename = csvfileset.streetmapbasename
        validfinds = f_regions(csvfileset)

        if streetmapbasename == "skip"
            continue
        end

        println("Loading ", streetmapbasename)

        if !haskey(streetnet_cache, streetmapbasename)
            streetnet_cache[streetmapbasename] = load(STREETMAP_BASE*streetmapbasename*".jld")["streetmap"]
        end
        sn = streetnet_cache[streetmapbasename]

        tic()

        pdsetfile = joinpath(PRIMARYDATA_DIR, toext("primarydata_" * csvfile, "jld"))
        pdset = load(pdsetfile)["pdset"]

        featureset = gen_featureset_ego(pdset, validfinds, sn, FEATURES, filters=FILTERS)
        if isempty(aggregate_datafile)
            aggregate_datafile = featureset
        else
            aggregate_datafile = [aggregate_datafile, featureset]
        end
        println(size(aggregate_datafile))

        outfile = joinpath(featureset_dir, "featureset_" * csvfile * ".jld")
        JLD.save(outfile, "featureset", featureset)

        toc()
    end

    outfile = joinpath(FEATUREMATRIX_DIR, output_name*".jld")
    save(outfile, "data", aggregate_datafile)
end

print("TOTAL TIME: "); toc()