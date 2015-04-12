push!(LOAD_PATH, "path_to_base_files")

using StreetMap
using HDF5, JLD
using Features
using FilesystemUtils
using Trajdata

@everywhere using StreamStats
@everywhere using CarEM
@everywhere using Discretizers

@everywhere type OrigHistobinExtractResults
    histobin            :: Matrix{Float64}

    sample_count        :: Int
    n_lanechanges_ego   :: Int
    agg_speed_ego       :: StreamStats.Var
    agg_d_cl_ego        :: StreamStats.Var

    start_var_d_cl      :: StreamStats.Var
    start_var_yaw       :: StreamStats.Var
    start_var_turnrate  :: StreamStats.Var
    start_var_speed     :: StreamStats.Var
    start_var_accel     :: StreamStats.Var
    start_var_d_x_front :: StreamStats.Var
    start_var_v_x_front :: StreamStats.Var

    initial_scenes      :: Vector{Vector{VehicleInitialConditions}}
end

type CSVFileSet
    csvfile :: String
    streetmapbasename :: String

    lanechanges_normal    :: Vector{Int}
    lanechanges_postpass  :: Vector{Int}
    lanechanges_arbitrary :: Vector{Int}
    carfollow             :: Vector{Int}
    freeflow              :: Vector{Int}
end

const INPUT_EMSTATS_FOLDER = "where_to_import_emstats"
const EMSTATS_ARR = [
    # "freeflow_desang250_acc250_medium_forced.jld",
    # "freeflow_desang250_acc250_large_forced.jld",
    # "freeflow_desang250_acc250_medium_edge_addition.jld",
    # "freeflow_desang250_acc250_large_edge_addition.jld",
    # "freeflow_desang500_acc500_medium_forced.jld",
    # "freeflow_desang500_acc500_medium_edge_addition.jld",
    "following_desang250_acc250_medium_random_init.jld",
    "following_desang250_acc250_large_random_init.jld",
    # "following_desang250_acc250_large_edge_addition.jld",
    # "sustcross_desang250_acc250_medium_edge_addition.jld",
    # "sustcross_desang250_acc250_large_edge_addition.jld",
    # "lanechange_desang250_acc250_medium_forced.jld",
    ]

const OUTPUT_FOLDER_DIRICHLET   = "where_to_export_dirichlet"
const OUTPUT_FOLDER_CATEGORICAL = "where_to_export_categorical"

const SEC_PER_FRAME = 0.25 # seconds per frame to use in simulation
const SIM_HORIZON   = 4.0  # length of horizon [s]
const N_FRAMES      = int(SIM_HORIZON / SEC_PER_FRAME) + 2 # note: add some buffer frames to allow for initial conditions
const N_EULER_STEPS = 10
@everywhere const INITIAL_SPEED = 29.06 # [m/s] (65 mph)

const ROAD  = CarEM.StraightRoadway(3, 3.7)

const PRIMARYDATA_DIR = "where_to_import_primary_data"
const STREETMAP_BASE = "where_to_get_streetmaps"
const CSVFILESETS = (
            ("sample", "some_streetmap"),
            #=
            ACTUAL SET OF CSVFILESETS CENSORED
            =#
        )

function load_em(empstats_file::String)

    emstats = load(empstats_file)

    binmapdict = emstats["binmaps"]
    targets    = emstats["targets"]
    indicators = emstats["indicators"]
    stats      = emstats["statistics"]
    adj        = emstats["adjacency"]

    use_250ms = in(:f_accel_250ms, targets)

    encounter_model(build_bn(targets, indicators, stats, adj), 
                        stats, binmapdict, targets, indicators)
end

const EM_FREEFLOW = load_em(INPUT_EMSTATS_FOLDER * "freeflow_desang250_acc250_medium_edge_addition.jld")

@everywhere type ParamsHistobin
    discx :: LinearDiscretizer
    discy :: LinearDiscretizer
end
@everywhere function calc_histobin(
    road            :: StraightRoadway,
    em              :: Any, # EM or ScenarioSelector
    target_histobin :: OrigHistobinExtractResults,
    sim_params      :: SimParams,
    histobin_params :: ParamsHistobin,
    drivelog        :: Matrix{Float64}
    )

    n_bins_x = nlabels(histobin_params.discx)
    n_bins_y = nlabels(histobin_params.discy)

    histobin = zeros(Float64, n_bins_x, n_bins_y)
    for scene in target_histobin.initial_scenes

        # scene = sample_initial_scene(target_histobin)
        initial_speed = scene[1].speed

        runlog = simulate!(scene, road, em, sim_params, frameind=2)

        Δx = (runlog[end,LOG_COL_X] - runlog[2,LOG_COL_X]) / initial_speed
        Δy =  runlog[end,LOG_COL_Y] - runlog[2,LOG_COL_Y]

        bin_x = encode(histobin_params.discx, Δx)
        bin_y = encode(histobin_params.discy, Δy)

        histobin[bin_x, bin_y] += 1.0
    end
    histobin
end
@everywhere function KL_divergence_categorical(A::Matrix{Float64}, B::Matrix{Float64})
    
    Ap = A ./ sum(A)
    Bp = B ./ sum(B)

    KL_divergence = 0.0
    for i = 1 : length(A)
        KL_divergence += Ap[i]*log(Ap[i]/Bp[i])
    end
    KL_divergence::Float64
end
@everywhere function KL_divergence_dirichlet(A::Matrix{Float64}, B::Matrix{Float64})
    α0 = sum(A)
    KL_divergence = 0.0
    KL_divergence += lgamma(α0)
    KL_divergence -= lgamma(sum(B))
    KL_divergence -= sum([lgamma(a) for a in A])
    KL_divergence += sum([lgamma(b) for b in B])
    for i = 1 : length(A)
        KL_divergence += (A[i] - B[i])*(digamma(A[i]) - digamma(α0))
    end
    KL_divergence::Float64
end

function coordinate_descent_histobin_fit(
    road            :: StraightRoadway,
    em              :: EM,
    target_histobin :: OrigHistobinExtractResults,
    histobin_params :: ParamsHistobin,
    sec_per_frame   :: Float64,
    n_frames        :: Int,
    n_euler_steps   :: Int,
    KLdiv_method    :: Symbol; # ∈ :Dirichlet, :Categorical
    verbosity       :: Int = 0
    )

    # TODO(tim): also optimize n_euler_steps?

    if KLdiv_method == :Dirichlet
        KL_div_function = KL_divergence_dirichlet
    elseif KLdiv_method == :Categorical
        KL_div_function = KL_divergence_categorical
    else
        error("unknown KL divergence method $KLdiv_method")
    end

    if verbosity > 0
        println("Coordinte descent $KLdiv_method")
        tic()
    end

    A = target_histobin.histobin .+ 1.0

    param_options = (
        [SampleUniform(), SampleUniformZeroBin()],
        [(:none, 0), (:SMA, 2), (:SMA, 3), (:SMA, 4), (:SMA, 5), (:SMA, 6), (:SMA, 7), (:SMA, 8), (:SMA, 9), (:SMA, 10), (:SMA, 11), (:SMA, 12),
                     (:WMA, 2), (:WMA, 3), (:WMA, 4), (:WMA, 5), (:WMA, 6), (:WMA, 7), (:WMA, 8), (:WMA, 9), (:WMA, 10), (:WMA, 11), (:WMA, 12)],
        [SampleUniform(), SampleUniformZeroBin()],
        [(:none, 0), (:SMA, 2), (:SMA, 3), (:SMA, 4), (:SMA, 5), (:SMA, 6), (:SMA, 7), (:SMA, 8), (:SMA, 9), (:SMA, 10), (:SMA, 11), (:SMA, 12),
                     (:WMA, 2), (:WMA, 3), (:WMA, 4), (:WMA, 5), (:WMA, 6), (:WMA, 7), (:WMA, 8), (:WMA, 9), (:WMA, 10), (:WMA, 11), (:WMA, 12)]
    )
    n_params = length(param_options)
    drivelog = create_log(1, n_frames)

    params_tried = Set{Vector{Int}}() # param_inds
    paraminds = ones(Int, n_params)

    symbol_lat = symbol(get_target_lat(em))
    symbol_lon = symbol(get_target_lon(em))

    function inds_to_simparams(inds)
        smoothing_lat = param_options[2][inds[2]]
        smoothing_lon = param_options[4][inds[4]]
        sampling_lat = SamplingParams(param_options[1][inds[1]], smoothing_lat[1], smoothing_lat[2])
        sampling_lon = SamplingParams(param_options[3][inds[3]], smoothing_lon[1], smoothing_lon[2])
        SimParams(symbol_lat, symbol_lon, sampling_lat, sampling_lon, sec_per_frame, n_euler_steps, n_frames)
    end

    converged = false
    best_KLdiv = Inf
    iter = 0
    while !converged
        iter += 1
        if verbosity > 0
            println("iteration ", iter)
            toc()
            tic()
        end
        converged = true
        for coordinate = 1 : n_params
            if verbosity > 1
                println("\tcoordinate ", coordinate)
            end

            to_try = SimParams[]
            ip_arr = Int[]
            for ip in 1 : length(param_options[coordinate])

                newparams = copy(paraminds)
                newparams[coordinate] = ip

                if !in(newparams, params_tried)
                    push!(params_tried, newparams)
                    push!(ip_arr, ip)
                    push!(to_try, inds_to_simparams(newparams))
                end
            end

            if !isempty(to_try)

                tic()
                KL_divs = pmap(to_try) do sim_params

                    B = calc_histobin(road, em, target_histobin, sim_params, histobin_params, drivelog)
                    B .+= 1.0
                    KL_div_function(A,B)
                end
                toc()

                if verbosity > 1
                    for (i, kl_div) in enumerate(KL_divs)
                        ip = ip_arr[i]
                        println("\t coordinate -> ", param_options[coordinate][ip], "  ", kl_div)
                    end
                end

                KL_divergence, best_ind = findmin(KL_divs)
                if KL_divergence < best_KLdiv
                    best_KLdiv = KL_divergence
                    paraminds[coordinate] = ip_arr[best_ind]
                    converged = false
                    if verbosity > 0
                        println("\tfound better: ", coordinate, " -> ", param_options[coordinate][ip_arr[best_ind]])
                        println("\t\tKL: ", best_KLdiv)
                    end
                end
            end
        end
    end

    if verbosity > 0
        toc()
        println("optimal params: ", paraminds)
        println("KL-divergence: ", best_KLdiv)
    end

    (inds_to_simparams(paraminds), best_KLdiv)
end
function coordinate_descent_histobin_fit(
    road            :: StraightRoadway,
    ss              :: ScenarioSelector,
    target_histobin :: OrigHistobinExtractResults,
    histobin_params :: ParamsHistobin,
    sec_per_frame   :: Float64,
    n_frames        :: Int,
    n_euler_steps   :: Int,
    KLdiv_method    :: Symbol; # ∈ :Dirichlet, :Categorical
    verbosity       :: Int = 0
    )

    # TODO(tim): also optimize n_euler_steps?

    if KLdiv_method == :Dirichlet
        KL_div_function = KL_divergence_dirichlet
    elseif KLdiv_method == :Categorical
        KL_div_function = KL_divergence_categorical
    else
        error("unknown KL divergence method $KLdiv_method")
    end

    if verbosity > 0
        println("Coordinte descent $KLdiv_method")
        tic()
    end

    A = target_histobin.histobin .+ 1.0

    param_options = (
        [SampleUniform(), SampleUniformZeroBin()],
        [(:none, 0), (:SMA, 2), (:SMA, 3), (:SMA, 4), (:SMA, 5), (:SMA, 6), (:SMA, 7), (:SMA, 8), (:SMA, 9), (:SMA, 10), (:SMA, 11), (:SMA, 12),
                     (:WMA, 2), (:WMA, 3), (:WMA, 4), (:WMA, 5), (:WMA, 6), (:WMA, 7), (:WMA, 8), (:WMA, 9), (:WMA, 10), (:WMA, 11), (:WMA, 12)],
        [SampleUniform(), SampleUniformZeroBin()],
        [(:none, 0), (:SMA, 2), (:SMA, 3), (:SMA, 4), (:SMA, 5), (:SMA, 6), (:SMA, 7), (:SMA, 8), (:SMA, 9), (:SMA, 10), (:SMA, 11), (:SMA, 12),
                     (:WMA, 2), (:WMA, 3), (:WMA, 4), (:WMA, 5), (:WMA, 6), (:WMA, 7), (:WMA, 8), (:WMA, 9), (:WMA, 10), (:WMA, 11), (:WMA, 12)]
    )
    n_params = length(param_options)
    drivelog = create_log(1, n_frames)

    params_tried = Set{Vector{Int}}() # param_inds
    paraminds = ones(Int, n_params)

    em = ss.carfollow
    symbol_lat = symbol(get_target_lat(em))
    symbol_lon = symbol(get_target_lon(em))

    function inds_to_simparams(inds)
        smoothing_lat = param_options[2][inds[2]]
        smoothing_lon = param_options[4][inds[4]]
        sampling_lat = SamplingParams(param_options[1][inds[1]], smoothing_lat[1], smoothing_lat[2])
        sampling_lon = SamplingParams(param_options[3][inds[3]], smoothing_lon[1], smoothing_lon[2])
        SimParams(symbol_lat, symbol_lon, sampling_lat, sampling_lon, sec_per_frame, n_euler_steps, n_frames)
    end

    converged = false
    best_KLdiv = Inf
    iter = 0
    while !converged
        iter += 1
        if verbosity > 0
            println("iteration ", iter)
            toc()
            tic()
        end
        converged = true
        for coordinate = 1 : n_params
            if verbosity > 1
                println("\tcoordinate ", coordinate)
            end

            to_try = SimParams[]
            ip_arr = Int[]
            for ip in 1 : length(param_options[coordinate])

                newparams = copy(paraminds)
                newparams[coordinate] = ip

                if !in(newparams, params_tried)
                    push!(params_tried, newparams)
                    push!(ip_arr, ip)
                    push!(to_try, inds_to_simparams(newparams))
                end
            end

            if !isempty(to_try)

                tic()
                KL_divs = pmap(to_try) do sim_params

                    B = calc_histobin(road, ss, target_histobin, sim_params, histobin_params, drivelog)
                    B .+= 1.0
                    KL_div_function(A,B)
                end
                toc()

                if verbosity > 1
                    for (i, kl_div) in enumerate(KL_divs)
                        ip = ip_arr[i]
                        println("\t coordinate -> ", param_options[coordinate][ip], "  ", kl_div)
                    end
                end

                KL_divergence, best_ind = findmin(KL_divs)
                if KL_divergence < best_KLdiv
                    best_KLdiv = KL_divergence
                    paraminds[coordinate] = ip_arr[best_ind]
                    converged = false
                    if verbosity > 0
                        println("\tfound better: ", coordinate, " -> ", param_options[coordinate][ip_arr[best_ind]])
                        println("\t\tKL: ", best_KLdiv)
                    end
                end
            end
        end
    end

    if verbosity > 0
        toc()
        println("optimal params: ", paraminds)
        println("KL-divergence: ", best_KLdiv)
    end

    (inds_to_simparams(paraminds), best_KLdiv)
end
function extract_original_histobin(extract_stats::OrigHistobinExtractStats, histobin_params::ParamsHistobin)

    const FRAME_SKIP = iround(extract_stats.sim_horizon / 0.05)

    const TOLERANCE_D_CL      = extract_stats.tol_d_cl
    const TOLERANCE_YAW       = extract_stats.tol_yaw
    const TOLERANCE_TURNRATE  = extract_stats.tol_turnrate
    const TOLERANCE_SPEED     = extract_stats.tol_speed
    const TOLERANCE_ACCEL     = extract_stats.tol_accel
    const TOLERANCE_D_X_FRONT = extract_stats.tol_d_x_front
    const TOLERANCE_D_V_FRONT = extract_stats.tol_d_v_front

    n_bins_x = nlabels(histobin_params.discx)
    n_bins_y = nlabels(histobin_params.discy)
    histobin = zeros(Float64, n_bins_x, n_bins_y)

    xrange = (histobin_params.discx.binedges[1], histobin_params.discx.binedges[end])
    yrange = (histobin_params.discy.binedges[1], histobin_params.discy.binedges[end])

    sample_count      = 0
    n_lanechanges_ego = 0
    agg_speed_ego     = StreamStats.Var()
    agg_d_cl_ego      = StreamStats.Var()
    agg_timegap       = StreamStats.Var()

    start_var_d_cl       = StreamStats.Var()
    start_var_yaw        = StreamStats.Var()
    start_var_turnrate   = StreamStats.Var()
    start_var_speed      = StreamStats.Var()
    start_var_accel      = StreamStats.Var()
    start_var_d_x_front  = StreamStats.Var()
    start_var_v_x_front  = StreamStats.Var()

    n_initial_scenes = 0
    initial_scenes = Array(Vector{VehicleInitialConditions}, 130_000)

    streetnet_cache = Dict{String, StreetNetwork}()

    for (csvfile, streetmapbasename) in CSVFILESETS

        if streetmapbasename == "skip"
            continue
        end

        if !haskey(streetnet_cache, streetmapbasename)
            streetnet_cache[streetmapbasename] = load(STREETMAP_BASE*streetmapbasename*".jld")["streetmap"]
        end
        sn = streetnet_cache[streetmapbasename]

        pdsetfile = joinpath(PRIMARYDATA_DIR, toext("primarydata_" * csvfile, "jld"))
        pdset = load(pdsetfile)["pdset"]

        for validfind = 1 : nvalidfinds(pdset)

            validfind_fut = jumpframe(pdset, validfind, FRAME_SKIP)
            if validfind_fut != 0

                passes_subset = true
                for subset in extract_stats.subsets
                    if !bool(get(subset, pdset, sn, CARIND_EGO, validfind)::Float64)
                        passes_subset = false
                        break
                    end
                end
                if !passes_subset
                    continue
                end

                findA = validfind2frameind(pdset, validfind)
                findB = validfind2frameind(pdset, validfind_fut)

                posFy  = gete(pdset, :posFy, findA)
                v_orig = gete(pdset, :velFx, findA)
                ϕ      = gete(pdset, :posFyaw, findA)
                d_cl   = gete(pdset, :d_cl, findA)
                speed  = get(SPEED,    pdset, sn, CARIND_EGO, validfind)::Float64
                ω      = get(TURNRATE, pdset, sn, CARIND_EGO, validfind)::Float64
                a      = get(ACC,      pdset, sn, CARIND_EGO, validfind)::Float64
                has_front = bool(get(HAS_FRONT, pdset, sn, CARIND_EGO, validfind)::Float64)
                d_x_front = get(D_X_FRONT, pdset, sn, CARIND_EGO, validfind)::Float64
                v_x_front = get(V_X_FRONT, pdset, sn, CARIND_EGO, validfind)::Float64
                nll    = get(N_LANE_L, pdset, sn, CARIND_EGO, validfind)::Float64
                nlr    = get(N_LANE_R, pdset, sn, CARIND_EGO, validfind)::Float64
                
                if  abs(d_cl) ≤ TOLERANCE_D_CL      && 
                    abs(ϕ)    ≤ TOLERANCE_YAW       &&
                    abs(ω)    ≤ TOLERANCE_TURNRATE  &&
                    abs(a)    ≤ TOLERANCE_ACCEL     &&
                    d_x_front ≤ TOLERANCE_D_X_FRONT &&
                    v_x_front ≤ TOLERANCE_D_V_FRONT &&
                    (nll + nlr) > 0

                    posGxA = gete(pdset, :posGx, findA)
                    posGyA = gete(pdset, :posGy, findA)
                    posGxB = gete(pdset, :posGx, findB)
                    posGyB = gete(pdset, :posGy, findB)

                    (Δx, Δy) = frenet_distance_between_points(sn, posGxA, posGyA, posGxB, posGyB)
                    Δx /= v_orig
                    Δy -= posFy

                    if !isnan(Δx) && xrange[1] ≤ Δx ≤ xrange[2] && yrange[1] ≤ Δy ≤ yrange[2]

                        framecount = 0
                        inner_var_speed    = StreamStats.Var()
                        inner_var_abs_d_cl = StreamStats.Var()
                        inner_var_timegap  = StreamStats.Var()
                        
                        last_d_cl = d_cl
                        had_lanechange = false
                        for i = 0 : FRAME_SKIP
                            v = jumpframe(pdset, validfind, i)
                            if v != 0
                                framecount += 1
                                next_d_cl = get(D_CL, pdset, sn, CARIND_EGO, int(v))::Float64
                                if abs(next_d_cl - last_d_cl) > 1.0
                                    had_lanechange = true
                                end

                                v_ego = get(SPEED, pdset, sn, CARIND_EGO, int(v))::Float64
                                dx_front = get(D_X_FRONT, pdset, sn, CARIND_EGO, int(v))::Float64

                                update!(inner_var_speed, v_ego)
                                update!(inner_var_abs_d_cl, next_d_cl)
                                if !isinf(dx_front) && dx_front > 0
                                    update!(inner_var_timegap, dx_front / v_ego)
                                end
                                last_d_cl = next_d_cl
                            end
                        end

                        sample_count += 1
                        n_lanechanges_ego += had_lanechange
                        update!(agg_speed_ego, mean(inner_var_speed))
                        update!(agg_d_cl_ego,  mean(inner_var_abs_d_cl))
                        update!(agg_timegap,   mean(inner_var_timegap))

                        update!(start_var_d_cl,       d_cl)
                        update!(start_var_yaw,        ϕ)
                        update!(start_var_turnrate,   ω)
                        update!(start_var_speed,      speed)
                        update!(start_var_accel,      a)
                        if !isnan(d_x_front) && !isinf(d_x_front)
                            update!(start_var_d_x_front,  d_x_front)
                            update!(start_var_v_x_front,  v_x_front)
                        end

                        bin_x = encode(histobin_params.discx, Δx)
                        bin_y = encode(histobin_params.discy, Δy)

                        histobin[bin_x, bin_y] += 1.0


                        if nll == 0 && nlr == 0
                            lanecenter = 3.7
                        elseif nll == 0
                            lanecenter = 3.7*2
                        elseif nlr == 0
                            lanecenter = 0.0
                        else
                            lanecenter = 3.7
                        end

                        ego_vehicle = VehicleInitialConditions(0.0,posFy+lanecenter,ϕ,speed,ω,a)

                        if !has_front
                            initial_scenes[n_initial_scenes+=1] = [ego_vehicle]
                        else
                            ind_front = int(get(Features.INDFRONT, pdset, sn, CARIND_EGO, validfind)::Float64)
                            front_y = get(D_CL, pdset, sn, ind_front, validfind)::Float64
                            front_ϕ = get(YAW, pdset, sn, ind_front, validfind)::Float64
                            front_v = get(SPEED, pdset, sn, ind_front, validfind)::Float64
                            front_ω = get(TURNRATE, pdset, sn, ind_front, validfind)::Float64
                            front_a = get(ACC, pdset, sn, ind_front, validfind)::Float64

                            front_ω = isinf(front_ω) ? 0.0 : front_ω
                            front_a = isinf(front_a) ? 0.0 : front_a

                            front_vehicle = VehicleInitialConditions(d_x_front,front_y+lanecenter,front_ϕ,
                                                                       front_v,front_ω,front_a)
                            initial_scenes[n_initial_scenes+=1] = [ego_vehicle, front_vehicle]
                        end
                    end
                end
            end
        end
    end

    println("agg stats:")    
    println("\tn_samples:         ", sample_count)
    println("\tn_lanechanges_ego: ", n_lanechanges_ego)
    println("\td_cl_ego:          ", mean(agg_d_cl_ego),  "  ", sqrt(agg_d_cl_ego .v_hat))
    println("\tspeed_ego:         ", mean(agg_speed_ego), "  ", sqrt(agg_speed_ego.v_hat))
    println("\ttimegap:           ", mean(agg_timegap),   "  ", sqrt(agg_timegap.v_hat))

    println("initial conditions:")
    println("\td_cl:      ", mean(start_var_d_cl     ), "  ", sqrt(start_var_d_cl.v_hat     ))
    println("\tyaw:       ", rad2deg(mean(start_var_yaw)),      "  ", rad2deg(sqrt(start_var_yaw.v_hat)))
    println("\tturnrate:  ", rad2deg(mean(start_var_turnrate)), "  ", rad2deg(sqrt(start_var_turnrate.v_hat)))
    println("\tspeed:     ", mean(start_var_speed    ), "  ", sqrt(start_var_speed.v_hat    ))
    println("\taccel:     ", mean(start_var_accel    ), "  ", sqrt(start_var_accel.v_hat    ))
    println("\td_x_front: ", mean(start_var_d_x_front), "  ", sqrt(start_var_d_x_front.v_hat))
    println("\tv_x_front: ", mean(start_var_v_x_front), "  ", sqrt(start_var_v_x_front.v_hat))

    OrigHistobinExtractResults(histobin, sample_count, n_lanechanges_ego, agg_d_cl_ego, agg_speed_ego,
                               start_var_d_cl, start_var_yaw, start_var_turnrate, start_var_speed,
                               start_var_accel, start_var_d_x_front, start_var_v_x_front, 
                               initial_scenes[1:n_initial_scenes])
end
function extract_original_histobin_carfollow(extract_stats::OrigHistobinExtractStats, histobin_params::ParamsHistobin)

    const FRAME_SKIP = iround(extract_stats.sim_horizon / 0.05)

    const TOLERANCE_D_CL      = extract_stats.tol_d_cl
    const TOLERANCE_YAW       = extract_stats.tol_yaw
    const TOLERANCE_TURNRATE  = extract_stats.tol_turnrate
    const TOLERANCE_SPEED     = extract_stats.tol_speed
    const TOLERANCE_ACCEL     = extract_stats.tol_accel
    const TOLERANCE_D_X_FRONT = extract_stats.tol_d_x_front
    const TOLERANCE_D_V_FRONT = extract_stats.tol_d_v_front

    n_bins_x = nlabels(histobin_params.discx)
    n_bins_y = nlabels(histobin_params.discy)
    histobin = zeros(Float64, n_bins_x, n_bins_y)

    xrange = (histobin_params.discx.binedges[1], histobin_params.discx.binedges[end])
    yrange = (histobin_params.discy.binedges[1], histobin_params.discy.binedges[end])

    sample_count      = 0
    n_lanechanges_ego = 0
    agg_speed_ego     = StreamStats.Var()
    agg_d_cl_ego      = StreamStats.Var()
    agg_timegap       = StreamStats.Var()

    start_var_d_cl       = StreamStats.Var()
    start_var_yaw        = StreamStats.Var()
    start_var_turnrate   = StreamStats.Var()
    start_var_speed      = StreamStats.Var()
    start_var_accel      = StreamStats.Var()
    start_var_d_x_front  = StreamStats.Var()
    start_var_v_x_front  = StreamStats.Var()

    n_initial_scenes = 0
    initial_scenes = Array(Vector{VehicleInitialConditions}, 130_000)

    streetnet_cache = Dict{String, StreetNetwork}()

    for csvfileset in (
        #= 
        ACTUAL LIST OF CSVFileSets CENSORED
        =#
        )

        csvfile = csvfileset.csvfile
        streetmapbasename = csvfileset.streetmapbasename

        if streetmapbasename == "skip"
            continue
        end

        if !haskey(streetnet_cache, streetmapbasename)
            streetnet_cache[streetmapbasename] = load(STREETMAP_BASE*streetmapbasename*".jld")["streetmap"]
        end
        sn = streetnet_cache[streetmapbasename]

        pdsetfile = joinpath(PRIMARYDATA_DIR, toext("primarydata_" * csvfile, "jld"))
        pdset = load(pdsetfile)["pdset"]

        validfinds = falses(nvalidfinds(pdset))
        validfind_regions = [csvfileset.carfollow]

        row_index = 0
        region_index_lo = 1
        while region_index_lo < length(validfind_regions)
            for validfind in validfind_regions[region_index_lo] : validfind_regions[region_index_lo+1]
                validfinds[validfind] = true
            end
            region_index_lo += 2
        end

        for validfind in  1: length(validfinds)
            if !validfinds[validfind]
                continue
            end

            validfind_fut = jumpframe(pdset, validfind, FRAME_SKIP)
            if validfind_fut != 0 && !validfinds[validfind_fut]

                findA = validfind2frameind(pdset, validfind)
                findB = validfind2frameind(pdset, validfind_fut)

                posFy  = gete(pdset, :posFy, findA)
                v_orig = gete(pdset, :velFx, findA)
                ϕ      = gete(pdset, :posFyaw, findA)
                d_cl   = gete(pdset, :d_cl, findA)
                speed  = get(SPEED,    pdset, sn, CARIND_EGO, validfind)::Float64
                ω      = get(TURNRATE, pdset, sn, CARIND_EGO, validfind)::Float64
                a      = get(ACC,      pdset, sn, CARIND_EGO, validfind)::Float64
                has_front = bool(get(HAS_FRONT, pdset, sn, CARIND_EGO, validfind)::Float64)
                d_x_front = get(D_X_FRONT, pdset, sn, CARIND_EGO, validfind)::Float64
                v_x_front = get(V_X_FRONT, pdset, sn, CARIND_EGO, validfind)::Float64
                nll    = get(N_LANE_L, pdset, sn, CARIND_EGO, validfind)::Float64
                nlr    = get(N_LANE_R, pdset, sn, CARIND_EGO, validfind)::Float64
                
                if  abs(d_cl) ≤ TOLERANCE_D_CL      && 
                    abs(ϕ)    ≤ TOLERANCE_YAW       &&
                    abs(ω)    ≤ TOLERANCE_TURNRATE  &&
                    abs(a)    ≤ TOLERANCE_ACCEL     &&
                    d_x_front ≤ TOLERANCE_D_X_FRONT &&
                    v_x_front ≤ TOLERANCE_D_V_FRONT &&
                    (nll + nlr) > 0

                    posGxA = gete(pdset, :posGx, findA)
                    posGyA = gete(pdset, :posGy, findA)
                    posGxB = gete(pdset, :posGx, findB)
                    posGyB = gete(pdset, :posGy, findB)

                    (Δx, Δy) = frenet_distance_between_points(sn, posGxA, posGyA, posGxB, posGyB)

                    Δx /= v_orig
                    Δy -= posFy

                    if !isnan(Δx) && xrange[1] ≤ Δx ≤ xrange[2] && yrange[1] ≤ Δy ≤ yrange[2]

                        framecount = 0
                        inner_var_speed    = StreamStats.Var()
                        inner_var_abs_d_cl = StreamStats.Var()
                        inner_var_timegap  = StreamStats.Var()
                        
                        last_d_cl = d_cl
                        had_lanechange = false
                        for i = 0 : FRAME_SKIP
                            v = jumpframe(pdset, validfind, i)
                            if v != 0
                                framecount += 1
                                next_d_cl = get(D_CL, pdset, sn, CARIND_EGO, int(v))::Float64
                                if abs(next_d_cl - last_d_cl) > 1.0
                                    had_lanechange = true
                                end

                                v_ego = get(SPEED, pdset, sn, CARIND_EGO, int(v))::Float64
                                dx_front = get(D_X_FRONT, pdset, sn, CARIND_EGO, int(v))::Float64

                                update!(inner_var_speed, v_ego)
                                update!(inner_var_abs_d_cl, next_d_cl)
                                if !isinf(dx_front) && dx_front > 0
                                    update!(inner_var_timegap, dx_front / v_ego)
                                end
                                last_d_cl = next_d_cl
                            end
                        end

                        if had_lanechange
                            continue
                        end

                        sample_count += 1
                        n_lanechanges_ego += had_lanechange
                        update!(agg_speed_ego, mean(inner_var_speed))
                        update!(agg_d_cl_ego,  mean(inner_var_abs_d_cl))
                        update!(agg_timegap,   mean(inner_var_timegap))

                        update!(start_var_d_cl,       d_cl)
                        update!(start_var_yaw,        ϕ)
                        update!(start_var_turnrate,   ω)
                        update!(start_var_speed,      speed)
                        update!(start_var_accel,      a)
                        if !isnan(d_x_front) && !isinf(d_x_front)
                            update!(start_var_d_x_front,  d_x_front)
                            update!(start_var_v_x_front,  v_x_front)
                        end

                        bin_x = encode(histobin_params.discx, Δx)
                        bin_y = encode(histobin_params.discy, Δy)

                        histobin[bin_x, bin_y] += 1.0


                        if nll == 0 && nlr == 0
                            lanecenter = 3.7
                        elseif nll == 0
                            lanecenter = 3.7*2
                        elseif nlr == 0
                            lanecenter = 0.0
                        else
                            lanecenter = 3.7
                        end

                        ego_vehicle = VehicleInitialConditions(0.0,posFy+lanecenter,ϕ,speed,ω,a)

                        if !has_front
                            initial_scenes[n_initial_scenes+=1] = [ego_vehicle]
                        else
                            ind_front = int(get(Features.INDFRONT, pdset, sn, CARIND_EGO, validfind)::Float64)
                            front_y = get(D_CL,     pdset, sn, ind_front, validfind)::Float64
                            front_ϕ = get(YAW,      pdset, sn, ind_front, validfind)::Float64
                            front_v = get(SPEED,    pdset, sn, ind_front, validfind)::Float64
                            front_ω = get(TURNRATE, pdset, sn, ind_front, validfind)::Float64
                            front_a = get(ACC,      pdset, sn, ind_front, validfind)::Float64

                            front_ω = isinf(front_ω) ? 0.0 : front_ω
                            front_a = isinf(front_a) ? 0.0 : front_a

                            front_vehicle = VehicleInitialConditions(d_x_front,front_y+lanecenter,front_ϕ,
                                                                       front_v,front_ω,front_a)
                            initial_scenes[n_initial_scenes+=1] = [ego_vehicle, front_vehicle]
                        end
                    end
                end
            end
        end
    end

    println("agg stats:")    
    println("\tn_samples:         ", sample_count)
    println("\tn_lanechanges_ego: ", n_lanechanges_ego)
    println("\td_cl_ego:          ", mean(agg_d_cl_ego),  "  ", sqrt(agg_d_cl_ego .v_hat))
    println("\tspeed_ego:         ", mean(agg_speed_ego), "  ", sqrt(agg_speed_ego.v_hat))
    println("\ttimegap:           ", mean(agg_timegap),   "  ", sqrt(agg_timegap.v_hat))

    println("initial conditions:")
    println("\td_cl:      ", mean(start_var_d_cl     ), "  ", sqrt(start_var_d_cl.v_hat     ))
    println("\tyaw:       ", rad2deg(mean(start_var_yaw)),      "  ", rad2deg(sqrt(start_var_yaw.v_hat)))
    println("\tturnrate:  ", rad2deg(mean(start_var_turnrate)), "  ", rad2deg(sqrt(start_var_turnrate.v_hat)))
    println("\tspeed:     ", mean(start_var_speed    ), "  ", sqrt(start_var_speed.v_hat    ))
    println("\taccel:     ", mean(start_var_accel    ), "  ", sqrt(start_var_accel.v_hat    ))
    println("\td_x_front: ", mean(start_var_d_x_front), "  ", sqrt(start_var_d_x_front.v_hat))
    println("\tv_x_front: ", mean(start_var_v_x_front), "  ", sqrt(start_var_v_x_front.v_hat))

    OrigHistobinExtractResults(histobin, sample_count, n_lanechanges_ego, agg_d_cl_ego, agg_speed_ego,
                               start_var_d_cl, start_var_yaw, start_var_turnrate, start_var_speed,
                               start_var_accel, start_var_d_x_front, start_var_v_x_front, 
                               initial_scenes[1:n_initial_scenes])
end
function extract_original_histobin_lanechange(extract_stats::OrigHistobinExtractStats, histobin_params::ParamsHistobin)

    const FRAME_SKIP = iround(extract_stats.sim_horizon / 0.05)

    const TOLERANCE_D_CL      = extract_stats.tol_d_cl
    const TOLERANCE_YAW       = extract_stats.tol_yaw
    const TOLERANCE_TURNRATE  = extract_stats.tol_turnrate
    const TOLERANCE_SPEED     = extract_stats.tol_speed
    const TOLERANCE_ACCEL     = extract_stats.tol_accel
    const TOLERANCE_D_X_FRONT = extract_stats.tol_d_x_front
    const TOLERANCE_D_V_FRONT = extract_stats.tol_d_v_front

    n_bins_x = nlabels(histobin_params.discx)
    n_bins_y = nlabels(histobin_params.discy)
    histobin = zeros(Float64, n_bins_x, n_bins_y)

    xrange = (histobin_params.discx.binedges[1], histobin_params.discx.binedges[end])
    yrange = (histobin_params.discy.binedges[1], histobin_params.discy.binedges[end])

    sample_count      = 0
    n_lanechanges_ego = 0
    agg_speed_ego     = StreamStats.Var()
    agg_d_cl_ego      = StreamStats.Var()
    agg_timegap       = StreamStats.Var()

    start_var_d_cl       = StreamStats.Var()
    start_var_yaw        = StreamStats.Var()
    start_var_turnrate   = StreamStats.Var()
    start_var_speed      = StreamStats.Var()
    start_var_accel      = StreamStats.Var()
    start_var_d_x_front  = StreamStats.Var()
    start_var_v_x_front  = StreamStats.Var()

    n_initial_scenes = 0
    initial_scenes = Array(Vector{VehicleInitialConditions}, 130_000)

    streetnet_cache = Dict{String, StreetNetwork}()

    for csvfileset in (
        #=
        ACTUAL LIST OF CSVFileSets CENSORED
        =#
        )

        csvfile = csvfileset.csvfile
        streetmapbasename = csvfileset.streetmapbasename

        if streetmapbasename == "skip"
            continue
        end

        if !haskey(streetnet_cache, streetmapbasename)
            streetnet_cache[streetmapbasename] = load(STREETMAP_BASE*streetmapbasename*".jld")["streetmap"]
        end
        sn = streetnet_cache[streetmapbasename]

        pdsetfile = joinpath(PRIMARYDATA_DIR, toext("primarydata_" * csvfile, "jld"))
        pdset = load(pdsetfile)["pdset"]

        validfinds = falses(nvalidfinds(pdset))
        validfind_regions = [csvfileset.lanechanges_normal, csvfileset.lanechanges_postpass, csvfileset.lanechanges_arbitrary]

        row_index = 0
        region_index_lo = 1
        while region_index_lo < length(validfind_regions)
            for validfind in validfind_regions[region_index_lo] : validfind_regions[region_index_lo+1]
                validfinds[validfind] = true
            end
            region_index_lo += 2
        end

        for validfind in  1: length(validfinds)
            if !validfinds[validfind]
                continue
            end

            validfind_fut = jumpframe(pdset, validfind, FRAME_SKIP)
            if validfind_fut != 0 && !validfinds[validfind_fut]

                findA = validfind2frameind(pdset, validfind)
                findB = validfind2frameind(pdset, validfind_fut)

                posFy  = gete(pdset, :posFy, findA)
                v_orig = gete(pdset, :velFx, findA)
                ϕ      = gete(pdset, :posFyaw, findA)
                d_cl   = gete(pdset, :d_cl, findA)
                speed  = get(SPEED,    pdset, sn, CARIND_EGO, validfind)::Float64
                ω      = get(TURNRATE, pdset, sn, CARIND_EGO, validfind)::Float64
                a      = get(ACC,      pdset, sn, CARIND_EGO, validfind)::Float64
                has_front = bool(get(HAS_FRONT, pdset, sn, CARIND_EGO, validfind)::Float64)
                d_x_front = get(D_X_FRONT, pdset, sn, CARIND_EGO, validfind)::Float64
                v_x_front = get(V_X_FRONT, pdset, sn, CARIND_EGO, validfind)::Float64
                nll    = get(N_LANE_L, pdset, sn, CARIND_EGO, validfind)::Float64
                nlr    = get(N_LANE_R, pdset, sn, CARIND_EGO, validfind)::Float64
                
                if  abs(d_cl) ≤ TOLERANCE_D_CL      && 
                    abs(ϕ)    ≤ TOLERANCE_YAW       &&
                    abs(ω)    ≤ TOLERANCE_TURNRATE  &&
                    abs(a)    ≤ TOLERANCE_ACCEL     &&
                    d_x_front ≤ TOLERANCE_D_X_FRONT &&
                    v_x_front ≤ TOLERANCE_D_V_FRONT &&
                    (nll + nlr) > 0

                    posGxA = gete(pdset, :posGx, findA)
                    posGyA = gete(pdset, :posGy, findA)
                    posGxB = gete(pdset, :posGx, findB)
                    posGyB = gete(pdset, :posGy, findB)

                    (Δx, Δy) = frenet_distance_between_points(sn, posGxA, posGyA, posGxB, posGyB)
                    Δx /= v_orig
                    Δy -= posFy

                    if !(abs(Δy - 3.7) < 0.5)
                        continue
                    end

                    if !isnan(Δx) && xrange[1] ≤ Δx ≤ xrange[2] && yrange[1] ≤ Δy ≤ yrange[2]

                        framecount = 0
                        inner_var_speed    = StreamStats.Var()
                        inner_var_abs_d_cl = StreamStats.Var()
                        inner_var_timegap  = StreamStats.Var()
                        
                        last_d_cl = d_cl
                        had_lanechange = false
                        for i = 0 : FRAME_SKIP
                            v = jumpframe(pdset, validfind, i)
                            if v != 0
                                framecount += 1
                                next_d_cl = get(D_CL, pdset, sn, CARIND_EGO, int(v))::Float64
                                if abs(next_d_cl - last_d_cl) > 1.0
                                    had_lanechange = true
                                end

                                v_ego = get(SPEED, pdset, sn, CARIND_EGO, int(v))::Float64
                                dx_front = get(D_X_FRONT, pdset, sn, CARIND_EGO, int(v))::Float64

                                update!(inner_var_speed, v_ego)
                                update!(inner_var_abs_d_cl, next_d_cl)
                                if !isinf(dx_front) && dx_front > 0
                                    update!(inner_var_timegap, dx_front / v_ego)
                                end
                                last_d_cl = next_d_cl
                            end
                        end

                        sample_count += 1
                        n_lanechanges_ego += had_lanechange
                        update!(agg_speed_ego, mean(inner_var_speed))
                        update!(agg_d_cl_ego,  mean(inner_var_abs_d_cl))
                        update!(agg_timegap,   mean(inner_var_timegap))

                        update!(start_var_d_cl,       d_cl)
                        update!(start_var_yaw,        ϕ)
                        update!(start_var_turnrate,   ω)
                        update!(start_var_speed,      speed)
                        update!(start_var_accel,      a)
                        update!(start_var_d_x_front,  d_x_front)
                        update!(start_var_v_x_front,  v_x_front)

                        bin_x = encode(histobin_params.discx, Δx)
                        bin_y = encode(histobin_params.discy, Δy)

                        histobin[bin_x, bin_y] += 1.0


                        if nll == 0 && nlr == 0
                            lanecenter = 3.7
                        elseif nll == 0
                            lanecenter = 3.7*2
                        elseif nlr == 0
                            lanecenter = 0.0
                        else
                            lanecenter = 3.7
                        end

                        ego_vehicle = VehicleInitialConditions(0.0,posFy+lanecenter,ϕ,speed,ω,a)

                        if !has_front
                            initial_scenes[n_initial_scenes+=1] = [ego_vehicle]
                        else
                            ind_front = int(get(Features.INDFRONT, pdset, sn, CARIND_EGO, validfind)::Float64)
                            front_y = get(D_CL, pdset, sn, ind_front, validfind)::Float64
                            front_ϕ = get(YAW, pdset, sn, ind_front, validfind)::Float64
                            front_v = get(SPEED, pdset, sn, ind_front, validfind)::Float64
                            front_ω = get(TURNRATE, pdset, sn, ind_front, validfind)::Float64
                            front_a = get(ACC, pdset, sn, ind_front, validfind)::Float64

                            front_ω = isinf(front_ω) ? 0.0 : front_ω
                            front_a = isinf(front_a) ? 0.0 : front_a

                            front_vehicle = VehicleInitialConditions(d_x_front,front_y+lanecenter,front_ϕ,
                                                                       front_v,front_ω,front_a)
                            initial_scenes[n_initial_scenes+=1] = [ego_vehicle, front_vehicle]
                        end
                    end
                end
            end
        end
    end

    println("agg stats:")    
    println("\tn_samples:         ", sample_count)
    println("\tn_lanechanges_ego: ", n_lanechanges_ego)
    println("\td_cl_ego:          ", mean(agg_d_cl_ego),  "  ", sqrt(agg_d_cl_ego .v_hat))
    println("\tspeed_ego:         ", mean(agg_speed_ego), "  ", sqrt(agg_speed_ego.v_hat))
    println("\ttimegap:           ", mean(agg_timegap),   "  ", sqrt(agg_timegap.v_hat))

    println("initial conditions:")
    println("\td_cl:      ", mean(start_var_d_cl     ), "  ", sqrt(start_var_d_cl.v_hat     ))
    println("\tyaw:       ", rad2deg(mean(start_var_yaw)),      "  ", rad2deg(sqrt(start_var_yaw.v_hat)))
    println("\tturnrate:  ", rad2deg(mean(start_var_turnrate)), "  ", rad2deg(sqrt(start_var_turnrate.v_hat)))
    println("\tspeed:     ", mean(start_var_speed    ), "  ", sqrt(start_var_speed.v_hat    ))
    println("\taccel:     ", mean(start_var_accel    ), "  ", sqrt(start_var_accel.v_hat    ))
    println("\td_x_front: ", mean(start_var_d_x_front), "  ", sqrt(start_var_d_x_front.v_hat))
    println("\tv_x_front: ", mean(start_var_v_x_front), "  ", sqrt(start_var_v_x_front.v_hat))

    OrigHistobinExtractResults(histobin, sample_count, n_lanechanges_ego, agg_d_cl_ego, agg_speed_ego,
                               start_var_d_cl, start_var_yaw, start_var_turnrate, start_var_speed,
                               start_var_accel, start_var_d_x_front, start_var_v_x_front, 
                               initial_scenes[1:n_initial_scenes])
end

function check_feature_extract_implemented(em::EM)

    runlog = create_log(1,1)
    road   = StraightRoadway(1,1.0)

    unsupported = AbstractFeature[]
    for f in get_indicators(em)
        CarEM.clear_unimp()
        get(f, runlog, road, 1.0, 1, 1)
        if CarEM.unimplemented
            push!(unsupported, f)
        end
    end
    unsupported
end
function check_feature_extract_implemented{S<:String}(emstats_arr::Vector{S})
    unsupported = Set{AbstractFeature}()
    for emstats_file in EMSTATS_ARR
        em = load_em(INPUT_EMSTATS_FOLDER * emstats_file)
        for f in check_feature_extract_implemented(em)
            push!(unsupported, f)
        end
    end
    for f in unsupported
        println("get for $(symbol(f)) unsupported!")
    end
    return isempty(unsupported)
end

histobin_params = ParamsHistobin(
                    LinearDiscretizer(linspace(SIM_HORIZON - 1.0,SIM_HORIZON + 1.0,25)),
                    LinearDiscretizer(linspace(-3.5, 3.5,25)))

const TARGET_HISTOBIN_EXTRACT_STATS = OrigHistobinExtractStats(
        # "freeflow",  AbstractFeature[SUBSET_FREE_FLOW, SUBSET_AT_SIXTYFIVE], 0.25, deg2rad(1.0), deg2rad(0.5), Inf, 0.05, Inf, Inf, SIM_HORIZON
        "carfollow", AbstractFeature[],   0.75, Inf, Inf, Inf, Inf, Inf, Inf, SIM_HORIZON
        # "lanechange", AbstractFeature[], 0.5, Inf, Inf, Inf, Inf, Inf, Inf, SIM_HORIZON
    )
const TARGET_HISTOBIN_FILE = "/media/tim/DATAPART1/Data/Bosch/processed/plots/sim/target_histobin_" * TARGET_HISTOBIN_EXTRACT_STATS.name * ".jld"

# orig_histobin_res = extract_original_histobin(TARGET_HISTOBIN_EXTRACT_STATS, histobin_params) 
orig_histobin_res = extract_original_histobin_carfollow(TARGET_HISTOBIN_EXTRACT_STATS, histobin_params)
# orig_histobin_res = extract_original_histobin_lanechange(TARGET_HISTOBIN_EXTRACT_STATS, histobin_params) 
save(TARGET_HISTOBIN_FILE, "histobin", orig_histobin_res)

orig_histobin_res = load(TARGET_HISTOBIN_FILE, "histobin")

histobin_map_dirichlet = Dict{String, Matrix{Float64}}()
histobin_map_categorical = Dict{String, Matrix{Float64}}()

if !check_feature_extract_implemented(EMSTATS_ARR)
    println("Terminating due to unsupported get")
    quit()
end

fout = open(OUTPUT_FOLDER_DIRICHLET * "output.txt", "w")
for emstats_file in EMSTATS_ARR
    println("loading ", emstats_file)

    print(fout, emstats_file, ", ")

    em = load_em(INPUT_EMSTATS_FOLDER * emstats_file)
    println("counts: ", counts(em))

    scenario_selector = ScenarioSelector(EM_FREEFLOW, em, EM_FREEFLOW)

    sim_params, best_KLdiv = coordinate_descent_histobin_fit(ROAD, scenario_selector, orig_histobin_res, histobin_params,
                                SEC_PER_FRAME, N_FRAMES, N_EULER_STEPS, :Dirichlet, verbosity=2)

    histobin_map_dirichlet[emstats_file] = calc_histobin(ROAD, scenario_selector, orig_histobin_res, sim_params, histobin_params, create_log(1, N_FRAMES))

    aggmetrics = aggregate_metrics(orig_histobin_res.initial_scenes, ROAD, scenario_selector, sim_params)
    
    print(fout, best_KLdiv, ", ")
    print_results_csv_readable(fout, sim_params, aggmetrics)
    print(fout, "\n")
end
close(fout)

fout = open(OUTPUT_FOLDER_CATEGORICAL * "output.txt", "w")
for emstats_file in EMSTATS_ARR
    println("loading ", emstats_file)

    print(fout, emstats_file, ", ")

    em = load_em(INPUT_EMSTATS_FOLDER * emstats_file)
    println("counts: ", counts(em))

    scenario_selector = ScenarioSelector(EM_FREEFLOW, em, EM_FREEFLOW)

    sim_params, best_KLdiv = coordinate_descent_histobin_fit(ROAD, scenario_selector, orig_histobin_res, histobin_params,
                                SEC_PER_FRAME, N_FRAMES, N_EULER_STEPS, :Categorical, verbosity=2)
    
    histobin_map_categorical[emstats_file] = calc_histobin(ROAD, scenario_selector, orig_histobin_res, sim_params, histobin_params, create_log(1, N_FRAMES))

    aggmetrics = aggregate_metrics(orig_histobin_res.initial_scenes, ROAD, scenario_selector, sim_params)
    
    print(fout, best_KLdiv, ", ")
    print_results_csv_readable(fout, sim_params, aggmetrics)
    print(fout, "\n")
end
close(fout)

using PGFPlots
function histobin_image(
    filename::String, 
    histobin::Matrix{Float64}, 
    params::ParamsHistobin; 
    log_intensity_scale :: Bool = false,
    div :: Float64 = 4.0
    )

    xrange = (params.discx.binedges[1], params.discx.binedges[end])
    yrange = (params.discy.binedges[1], params.discy.binedges[end])

    histobin .+= 1.0
    if log_intensity_scale
        histobin = Base.log10(histobin)
        histobin ./= div
    else
        histobin = histobin ./ maximum(histobin)
    end

    histobin = rotr90(histobin)
    histobin = 1.0 - histobin 

    println(extrema(histobin))
    
    @assert(maximum(histobin) ≤ 1.0)
    @assert(minimum(histobin) ≥ 0.0)

    p = PGFPlots.Image(histobin, yrange, xrange, scale_mat=false)
    ax = Axis(p, xlabel="Lateral Deviation", ylabel="Longitudinal Deviation")
    PGFPlots.save(filename, ax)
end

max_count = maximum(orig_histobin_res.histobin)
for emstats_file in EMSTATS_ARR
    max_count = max(max_count, maximum(histobin_map_dirichlet[emstats_file]))
end
div = Base.log10(max_count+1)

histobin_image("where_to_export" * TARGET_HISTOBIN_EXTRACT_STATS.name * ".pdf", orig_histobin_res.histobin, histobin_params, log_intensity_scale=true, div=div)

for emstats_file in EMSTATS_ARR

    outfile = OUTPUT_FOLDER_DIRICHLET * emstats_file[1:end-4] * "_histobin.pdf"
    sim_histobin = histobin_map_dirichlet[emstats_file]
    histobin_image(outfile, sim_histobin, histobin_params, log_intensity_scale=true, div=div)

    outfile = OUTPUT_FOLDER_CATEGORICAL * emstats_file[1:end-4] * "_histobin.pdf"
    sim_histobin = histobin_map_categorical[emstats_file]
    histobin_image(outfile, sim_histobin, histobin_params, log_intensity_scale=true)
end

