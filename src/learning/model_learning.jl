
# -------------------------------------------------------

push!(LOAD_PATH, "/home/tim/Documents/wheelerworkspace/Bosch/model/")

# -------------------------------------------------------

using DataFrames
using StreetMap
using Features
using BinMaps
using Smile
using SmileExtra

using BayesNets
using Graphs
using CarEM
using HDF5, JLD

import BayesNets: statistics

# --------

const BASE_OUTPUT_FOLDER = "/media/tim/DATAPART1/Data/Bosch/processed/plots/"

type ModelSet
	name :: String
	filters :: Vector{AbstractFeature}
end
type ModelTargets
	lat :: AbstractFeature
	lon :: AbstractFeature
end
type ParentFeatures
	lat :: Vector{AbstractFeature}
	lon :: Vector{AbstractFeature}
end

type GraphLearningResult
	fileroot     :: String
	target_lat   :: AbstractFeature
	target_lon   :: AbstractFeature     
	parents_lat  :: Vector{AbstractFeature}
	parents_lon  :: Vector{AbstractFeature}
	features     :: Vector{AbstractFeature}
	adj          :: BitMatrix
	net          :: Network
	stats        :: Vector{Matrix{Float64}} # NOTE(tim): this does not include prior counts
	score_raw    :: Float64
	score_K2     :: Float64
	score_BDeu   :: Float64

	function GraphLearningResult(
		basefolder     :: String,
		features       :: Vector{AbstractFeature},
		ind_lat        :: Int,
		ind_lon        :: Int,
		parentinds_lat :: Vector{Int},
		parentinds_lon :: Vector{Int},
		score_raw      :: Float64,
		mat_full       :: Matrix{Int},
		s_r            :: SharedVector{Int},
		s_d            :: SharedMatrix{Int}
		)
		
		parents_lat = features[parentinds_lat]
		parents_lon = features[parentinds_lon]

		f_inds = feature_indeces_in_net(features, ind_lat, ind_lon, parentinds_lat, parentinds_lon)
		net_features = features[f_inds]
		net   = to_net(features, parentinds_lat, parentinds_lon, mat_full, f_inds)
		adj   = SmileExtra.adjacency_matrix(net)
		stats = convert(Vector{Matrix{Float64}}, SmileExtra.statistics(adj, s_r[f_inds], s_d[f_inds,:]))

		score_K2   = net_score(net, s_r, s_d)
		score_BDeu = net_score_BDeu(net, s_r, s_d)

		target_lat = features[ind_lat]
		target_lon = features[ind_lon]

		new(basefolder, target_lat, target_lon, parents_lat, parents_lon, 
			net_features, adj, net, stats, score_raw, score_K2, score_BDeu)
	end
end

const DECIMATION = 1               # NOTE(tim): a negative decimation ends up being a multiplication factor
const RANDOMINIT_N_INIT = 100
const RANDOMINIT_N_PARENTS = 2
const IMPROVEMENT_THRESHOLD = 1e-2 # new_score - cur_score > threshold to add feature
                                   # NOTE(tim): this prevents a feature from being added due to roundoff error

# --------

function create_output_folder(
	outputfolder::String = strftime("graph_feature_selection_%Y%m%d_%H%M/", time())
	)
	
	if !isabspath(outputfolder)
		outputfolder = joinpath(BASE_OUTPUT_FOLDER, outputfolder)
	end

	if !isdir(outputfolder)
		mkdir(outputfolder)
	end

	outputfolder
end
function pull_relevant_features{F<:AbstractFeature}(data::DataFrame, features::Vector{F})
	symbs = map(f->symbol(f), features)
	data[symbs] # pull all of the features we could want
end
function remove_na_target_values(data::DataFrame, target::AbstractFeature)
	n = size(data, 1)
	inds = trues(n)
	for i = 1 : n
		val = data[i, target]
		if isna(val) || val == NA_ALIAS || isnan(val)
			inds[i] = false
		end
	end
	data[inds,:]
end
function remove_na_target_values{F<:AbstractFeature}(data::DataFrame, targets::Vector{F})
	n,m = size(data,1), size(data,2)
	inds = trues(n)
	for i = 1 : n
		for f in targets
			val = data[i, symbol(f)]
			if isna(val) || val == NA_ALIAS || isnan(val)
				inds[i] = false
				break
			end
		end
	end
	data[inds,:]
end
function filter_features{F<:AbstractFeature}(data::DataFrame, filterfeatures::Vector{F})
	if isempty(filterfeatures)
		return data
	end

	n = size(data, 1)
	inds = trues(n)
	for i = 1 : n
		for f in filterfeatures
			val = data[i, symbol(f)]
			if isapprox(val, 0.0) || isna(val) || val == NA_ALIAS || isnan(val)
				inds[i] = false
				break
			end
		end
	end
	data[inds,:]
end


function short_name(F::AbstractFeature)
	if F == FUTURETURNRATE_250MS
		return "turnrate250"
	elseif F == FUTURETURNRATE_500MS
		return "turnrate500"
	elseif F == FUTUREDESIREDANGLE_250MS
		return "desang250"
	elseif F == FUTUREDESIREDANGLE_500MS
		return "desang500"
	elseif F == FUTUREACCELERATION_250MS
		return "acc250"
	elseif F == FUTUREACCELERATION_500MS
		return "acc500"
	elseif F == FUTUREDESIREDSPEED_250MS
		return "desspeed250"
	elseif F == FUTUREDESIREDSPEED_500MS
		return "desspeed500"
	else
		return "unknown"
	end
end
function discretize{T<:Integer}( binmaps, features, data, A::Type{T}=Uint8)
	# returns an M × N matrix where
	#    M = # of valid sample rows
	#    N = # of features

	M = size(data, 1)
	m = 1
	mat = Array(T, M, length(features))
	for i = 1 : M
		dropped = false
		
		for (j,f) in enumerate(features)
			sym = symbol(f)
			var  = data[i, sym]
			dmap = binmaps[sym]
			if (isa(var, NAtype) || isinf(var)) && !supports_encoding(dmap, NA)
				dropped = true
				break
			end
			if isa(var, NAtype) || isinf(var)
				mat[m,j] = convert(T, BinMaps.encode(dmap, NA))
			else
				mat[m,j] = convert(T, BinMaps.encode(dmap, var))
			end
		end

		if !dropped
			m += 1
		end
	end

	mat[1:m-1,:]
end
function writetex(net::Network, filename::String, features)

	latexnames = map(f->lsymbol(f), features)
	all_nodes  = get_all_nodes(net)
	G = simple_graph(length(all_nodes))
	for n in all_nodes
		for p in get_parents(net, n)
			add_edge!(G, p+1, n+1)
		end
	end


	splitext(filename)[2] == ".tex" || error("writetex only supports .tex format")

	fout = open(filename, "w")

	println(fout, "\\documentclass[tikz,border=10pt]{standalone}")
	println(fout, "\\usepackage{amsmath}")
	println(fout, "\\usetikzlibrary{graphdrawing}")
	println(fout, "\\usetikzlibrary{graphs}")
	println(fout, "\\usetikzlibrary{arrows.meta}")
	println(fout, "\\usegdlibrary{layered}")
	println(fout, "\\begin{document}")
	println(fout, "")
	println(fout, "\\begin{tikzpicture}[nodes={text height=.7em, text depth=.2em, draw=black!20, thick, fill=white, font=\\footnotesize}, rounded corners, semithick, >=stealth]")
	println(fout, "\\graph [layered layout, level distance=1cm, sibling sep=.5em, sibling distance=1cm] {")

	for n in topological_sort_by_dfs(G)
		println(fout, "\t\"", latexnames[n], "\" -> { ")
		for p in out_neighbors(n, G)
			println(fout, "\t\t\"", latexnames[p], "\", ")
		end
		println(fout, "\t};")
	end

	println(fout, "\t};")
	println(fout, "\\end{tikzpicture}")
	println(fout, "\\end{document}")
	close(fout)
end

function parent_indeces{F<:AbstractFeature}(parents::ParentFeatures, features::Vector{F})
	parents_lat = find(f->in(f, parents.lat), features)
	parents_lon = find(f->in(f, parents.lon), features)
	(parents_lat, parents_lon)
end

function feature_indeces_in_net(features, find_lat, find_lon, parents_lat, parents_lon)
	retval = [find_lat,find_lon]
	append!(retval, parents_lat)
	append!(retval, parents_lon)
	sort!(unique(retval))
end
function to_net(features, parents_lat, parents_lon, mat, feature_indeces_in_net)

	dset = Smile.matrix_to_dataset(mat[:,feature_indeces_in_net], 
		         map(i->string(symbol(features[i]))::String, feature_indeces_in_net))

	pat = Pattern()
	set_size(pat, length(feature_indeces_in_net))

	for i in parents_lat
		ind = findfirst(feature_indeces_in_net, i)
		set_edge(pat, ind-1, 0, DSL_EDGETYPE_DIRECTED)
	end
	for i in parents_lon
		ind = findfirst(feature_indeces_in_net, i)
		set_edge(pat, ind-1, 1, DSL_EDGETYPE_DIRECTED)
	end

	@assert(is_DAG(pat))

	net = Network()
	to_network(pat, dset, net)
	net
end
function print_parents(features, parents_lat, parents_lon)
	println("Final Structure: ")
	println("\tParents for ", symbol(features[1]), " (", length(parents_lat), ")")
	for i in parents_lat
		println("\t\t", symbol(features[i]))
	end
	println("\tParents for ", symbol(features[2]), " (", length(parents_lon), ")")
	for i in parents_lon
		println("\t\t", symbol(features[i]))
	end
end
function empty_states{S<:Real}(N::Matrix{S}, threshold=1)

	sum(row_sum->row_sum < threshold, sum(N,1))
end

function prior_counts_a{I<:Integer, J<:Integer}(
	i_afut :: Int, 
	i_velFx :: Int,
	parents_a :: AbstractVector{I}, 
	r :: AbstractVector{J},
	binmaps :: Dict{Symbol, AbstractBinMap};
	integration_discretization :: Float64 = 1e-3,
	prior_weight :: Float64 = 10.0
	)

	const v_limit = 29.0 # [m/s]
	const Δt = 0.25 # [s]

	@assert(in(i_velFx, parents_a)) # velFx must be a parent of afut

	pseudocounts = zeros(Float64, r[i_velFx], r[i_afut])

	for bin_velFx = 1 : r[i_velFx]
		velFx_lo, velFx_hi = bin_bounds(binmaps[:velFx], bin_velFx)
		n_counts = 0
		for velFx = velFx_lo : integration_discretization : velFx_hi
			afut = (v_limit - velFx)/Δt
			bin_afut  = encode(binmaps[:f_accel_250ms], afut)
			pseudocounts[bin_velFx, bin_afut] += 1.0
			n_counts += 1
		end
		pseudocounts[bin_velFx,:] .*= prior_weight / n_counts

		# println("bin: ", bin_velFx, "  ", (v_limit - velFx_lo)/Δt, " <-> ", (v_limit - velFx_hi)/Δt, ":  ", pseudocounts[bin_velFx,:])
	end
	
	# run through all parental instantiations and assign alpha
	dims = tuple(r[parents_a]...)
	ind_of_velFx_in_dims = findfirst(parents_a, i_velFx)
	n_parental_instantiations = prod(dims)
	alpha = zeros(Float64, r[i_afut], isempty(parents_a) ? 1 : n_parental_instantiations)
	for q = 1 : size(alpha,2)
		bin_velFx = ind2sub(dims, q)[ind_of_velFx_in_dims]
		alpha[:,q] = pseudocounts[bin_velFx, :]
	end

	# add uniform prior
	alpha .+= 1.0
	return alpha
end
function prior_counts_lon{I<:Integer, J<:Integer}(
	i_pfut::Int, 
	i_yaw::Int,
	i_d_cl::Int,
	parents_lon::AbstractVector{I},
	r::AbstractVector{J},
	binmaps::Dict{Symbol, AbstractBinMap};
	integration_discretization_yaw::Float64 = 1e-3,
	integration_discretization_d_cl::Float64 = 1e-3,
	prior_weight::Float64 = 10.0
	)

	const κ = 0.25 # [-]
	const ρ = 0.01 # [rad/m]
	const Δt = 0.25 # [s]

	@assert(in(i_yaw, parents_lon)) # velFx must be a parent of afut
	@assert(in(i_d_cl, parents_lon)) # velFx must be a parent of afut
	@assert(!in(i_pfut, parents_lon))

	pseudocounts = zeros(Float64, r[i_yaw], r[i_d_cl], r[i_pfut])

	for bin_yaw = 1 : r[i_yaw]
		yaw_lo, yaw_hi = bin_bounds(binmaps[:yaw], bin_yaw)
		for bin_d_cl = 1 : r[i_d_cl]
			d_cl_lo, d_cl_hi = bin_bounds(binmaps[:d_cl], bin_d_cl)

			n_counts = 0
			for yaw = yaw_lo : integration_discretization_yaw : yaw_hi
				for d_cl = d_cl_lo : integration_discretization_d_cl : d_cl_hi
					pfut = κ*(-ρ*d_cl-yaw)/Δt
					bin_pfut  = encode(binmaps[:f_turnrate_250ms], pfut)
					pseudocounts[bin_yaw, bin_d_cl, bin_pfut] += 1.0
					n_counts += 1
				end
			end
			pseudocounts[bin_yaw, bin_d_cl, :] .*= prior_weight / n_counts
		end
	end
	
	# run through all parental instantiations and assign alpha
	dims = tuple(r[parents_lon]...)
	ind_of_yaw_in_dims = findfirst(parents_lon, i_yaw)
	ind_of_d_cl_in_dims = findfirst(parents_lon, i_d_cl)
	n_parental_instantiations = prod(dims)
	alpha = zeros(Float64, r[i_pfut], isempty(parents_lon) ? 1 : n_parental_instantiations)
	for q = 1 : size(alpha,2)
		subs = ind2sub(dims, q)
		bin_yaw = subs[ind_of_yaw_in_dims]
		bin_d_cl = subs[ind_of_d_cl_in_dims]
		alpha[:,q] = pseudocounts[bin_yaw, bin_d_cl, :]
	end

	# add uniform prior
	alpha .+= 1.0
	return alpha
end
function prior_counts_a{I<:Integer, J<:Integer}(
	i_afut::Int, 
	i_velFx::Int,
	parents_a::AbstractVector{I}, 
	r::AbstractVector{J},
	binmaps::Dict{Symbol, AbstractBinMap},
	cache::Dict{Vector{Int}, Matrix{Float64}};
	integration_discretization::Float64 = 1e-3,
	prior_weight::Float64 = 10.0
	)
	
	if haskey(cache, parents_a)
        return cache[parents_a]
    end
    return (cache[parents_a] = prior_counts_a(i_afut, i_velFx, parents_a, r, binmaps, 
    	integration_discretization=integration_discretization, prior_weight=prior_weight))
end
function prior_counts_lon{I<:Integer, J<:Integer}(
	i_pfut::Int, 
	i_yaw::Int,
	i_d_cl::Int,
	parents_lon::AbstractVector{I},
	r::AbstractVector{J},
	binmaps::Dict{Symbol, AbstractBinMap},
	cache::Dict{Vector{Int}, Matrix{Float64}};
	integration_discretization_yaw::Float64 = 1e-3,
	integration_discretization_d_cl::Float64 = 1e-3,
	prior_weight::Float64 = 10.0
	)
	
	if haskey(cache, parents_lon)
        return cache[parents_lon]
    end
    return (cache[parents_lon] = prior_counts_lon(i_pfut, i_yaw, i_d_cl, parents_lon, r, binmaps, 
    	integration_discretization_yaw=integration_discretization_yaw, 
    	integration_discretization_d_cl=integration_discretization_d_cl,
    	prior_weight=prior_weight))
end

function learn_edge_addition{F<:AbstractFeature}(
	features        :: Vector{F},
	n_targets       :: Int,
	n_indicators    :: Int,
	forced          :: ParentFeatures,
	s_r             :: SharedVector{Int},
	s_d             :: SharedMatrix{Int},
	binmaps         :: Dict{Symbol,AbstractBinMap},
	score_cache_lat :: Dict{Vector{Int}, Float64},
	score_cache_lon :: Dict{Vector{Int}, Float64},
	ind_lat         :: Int,
	ind_lon         :: Int;
	verbosity       :: Int = 0
	)

	# edge addition
	#   - start with an empty network and continuously add the next best edge

	parents_lat, parents_lon = parent_indeces(forced, features)
	chosen_lat  = map(i->in(i, parents_lat), [1:n_indicators]+n_targets)
	chosen_lon  = map(i->in(i, parents_lon), [1:n_indicators]+n_targets)
	score_lat   = log_bayes_score_component(ind_lat, parents_lat, s_r, s_d, score_cache_lat)
	score_lon   = log_bayes_score_component(ind_lon, parents_lon, s_r, s_d, score_cache_lon)
	score       = score_lat + score_lon

	if verbosity > 0
		@printf("starting score: %15.3f %15.3f %15.3f\n", score_lat, score_lon, score)
	end

	n_iter = 0
	done = false
	while !done
		n_iter += 1

		if verbosity > 1
			println("iter $n_iter")
		end

		best_diff = IMPROVEMENT_THRESHOLD
		selected_lat = false
		best_p = 0
		best_score_lat = 0.0
		best_score_lon = 0.0
		best_parents_lat = parents_lat
		best_parents_lon = parents_lon

		# check edges for indicators -> lat
		for i = 1 : n_indicators
			# add edge if it does not exist
			if !chosen_lat[i]
				new_parents = sort!(push!(copy(parents_lat), n_targets+i))
				new_score_lat = log_bayes_score_component(ind_lat, new_parents, s_r, s_d, score_cache_lat)

				if new_score_lat - score_lat > best_diff
					selected_lat = true
					best_p = i
					best_score_lat = new_score_lat
					best_parents_lat = new_parents
					best_diff = new_score_lat - score_lat
					if verbosity > 2
						println("\tadding $(symbol(indicators[i])) -> lat ", new_score_lat - score_lat)
					end
				end
			end
		end

		# check edges for indicators -> lon
		for i = 1 : n_indicators
			# add edge if it does not exist
			if !chosen_lon[i]
				new_parents = sort!(push!(copy(parents_lon), n_targets+i))
				new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, score_cache_lon)
				if new_score_lon - score_lon > best_diff
					selected_lat = false
					best_p = i
					best_score_lon = new_score_lon
					best_parents_lon = new_parents
					best_diff = new_score_lon - score_lon
					if verbosity > 2
						println("\tadding $(symbol(indicators[i])) -> lon ", new_score_lon - score_lon)
					end
				end
			end
		end

		# check edge between lat <-> lon
		if !in(ind_lon, parents_lat) && !in(ind_lat, parents_lon)
			# lon -> lat
			new_parents = unshift!(copy(parents_lat), ind_lon)
			new_score_lat = log_bayes_score_component(ind_lat, new_parents, s_r, s_d, score_cache_lat)

			if new_score_lat - score_lat > best_diff
				selected_lat = true
				best_p = -1
				best_score_lat = new_score_lat
				best_parents_lat = new_parents
				best_diff = new_score_lat - score_lat
				if verbosity > 2
					println("\tadding lon -> lat ", new_score_lat - score_lat)
				end
			end

			# lat -> lon
			new_parents = unshift!(copy(parents_lon), ind_lat)
			new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, score_cache_lon)

			if new_score_lon - score_lon > best_diff
				selected_lat = false
				best_p = -1
				best_score_lon = new_score_lon
				best_parents_lon = new_parents
				best_diff = new_score_lon - score_lon
				if verbosity > 2
					println("\tadding lat -> lon ", new_score_lon - score_lon)
				end
			end
		end

		# select best
		if best_diff > IMPROVEMENT_THRESHOLD
			if selected_lat
				parents_lat = best_parents_lat
				score_lat = best_score_lat
				if best_p != -1
					if verbosity > 1
						println("\tChose ", symbol(indicators[best_p]), " -> lat")
					end
					chosen_lat[best_p] = true
				else
					if verbosity > 1
						println("\tChose lon -> lat")
					end
				end
			else
				parents_lon = best_parents_lon
				score_lon = best_score_lon
				if best_p != -1
					if verbosity > 1
						println("\tChose ", symbol(indicators[best_p]), " -> lon")
					end
					chosen_lon[best_p] = true
				else
					if verbosity > 1
						println("\tChose lat -> lon")
					end
				end
			end
			score = score_lat + score_lon

			if verbosity > 0
				@printf("iteration %d:    %15.3f %15.3f %15.3f\n", n_iter, score_lat, score_lon, score)
			end
		else
			done = true
		end
	end

	(parents_lat, parents_lon, score)
end
function learn_graph_traversal{F<:AbstractFeature}(
	features        :: Vector{F},
	n_targets       :: Int,
	n_indicators    :: Int,
	forced          :: ParentFeatures,
	s_r             :: SharedVector{Int},
	s_d             :: SharedMatrix{Int},
	binmaps         :: Dict{Symbol,AbstractBinMap},
	score_cache_lat :: Dict{Vector{Int}, Float64},
	score_cache_lon :: Dict{Vector{Int}, Float64},
	ind_lat         :: Int,
	ind_lon         :: Int;
	verbosity       :: Int = 0
	)

	# graph traversal
	#   - start with an empty network and continuously add or remove the next best edge

	parents_lat, parents_lon = parent_indeces(forced, features)
	chosen_lat  = map(i->in(n_targets+i, parents_lat), [1:n_indicators])
	chosen_lon  = map(i->in(n_targets+i, parents_lon), [1:n_indicators])
	score_lat   = log_bayes_score_component(ind_lat, parents_lat, s_r, s_d, score_cache_lat)
	score_lon   = log_bayes_score_component(ind_lon, parents_lon, s_r, s_d, score_cache_lon)
	score       = score_lat + score_lon

	if verbosity > 0
		@printf("starting score: %15.3f %15.3f %15.3f\n", score_lat, score_lon, score)
	end

	n_iter = 0
	done = false
	while !done
		n_iter += 1

		best_diff = IMPROVEMENT_THRESHOLD
		selected_lat = false
		best_score_lat = 0.0
		best_score_lon = 0.0
		best_parents_lat = parents_lat
		best_parents_lon = parents_lon

		# check edges for indicators -> lat
		for i = 1 : n_indicators
			# add edge if it does not exist
			if !chosen_lat[i]
				new_parents = sort!(push!(copy(parents_lat), n_targets+i))
				new_score_lat = log_bayes_score_component(ind_lat, new_parents, s_r, s_d, score_cache_lat)
				if new_score_lat - score_lat > best_diff
					selected_lat = true
					best_score_lat = new_score_lat
					best_parents_lat = new_parents
					best_diff = new_score_lat - score_lat
				end
			end
		end
		for (idx, i) in enumerate(parents_lat)
			# remove edge if it does exist
			if !in(features[i], forced.lat)
				new_parents = deleteat!(copy(parents_lat), idx)
				new_score_lat = log_bayes_score_component(ind_lat, new_parents, s_r, s_d, score_cache_lat)
				if new_score_lat - score_lat > best_diff
					selected_lat = true
					best_score_lat = new_score_lat
					best_parents_lat = new_parents
					best_diff = new_score_lat - score_lat
				end
			end
		end

		# check edges for indicators -> lon
		for i = 1 : n_indicators
			# add edge if it does not exist
			if !chosen_lon[i]
				new_parents = sort!(push!(copy(parents_lon), n_targets+i))
				new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, score_cache_lon)
				if new_score_lon - score_lon > best_diff
					selected_lat = false
					best_score_lon = new_score_lon
					best_parents_lon = new_parents
					best_diff = new_score_lon - score_lon
				end
			end
		end
		for (idx, i) in enumerate(parents_lon)
			# remove edge if it does exist
			if !in(features[i], forced.lon)
				new_parents = deleteat!(copy(parents_lon), idx)
				new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, score_cache_lon)
				if new_score_lon - score_lon > best_diff
					selected_lat = false
					best_score_lon = new_score_lon
					best_parents_lon = new_parents
					best_diff = new_score_lon - score_lon
				end
			end
		end

		# check edge between lat <-> lon
		if !in(ind_lon, parents_lat) && !in(ind_lat, parents_lon)
			# lon -> lat
			new_parents = unshift!(copy(parents_lat), ind_lon)
			new_score_lat = log_bayes_score_component(ind_lat, new_parents, s_r, s_d, score_cache_lat)
			if new_score_lat - score_lat > best_diff
				selected_lat = true
				best_score_lat = new_score_lat
				best_parents_lat = new_parents
				best_diff = new_score_lat - score_lat
			end

			# lat -> lon
			new_parents = unshift!(copy(parents_lon), ind_lat)
			new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, score_cache_lon)
			if new_score_lon - score_lon > best_diff
				selected_lat = false
				best_score_lon = new_score_lon
				best_parents_lon = new_parents
				best_diff = new_score_lon - score_lon
			end
		elseif in(ind_lon, parents_lat) && !in(features[ind_lon], forced.lat)

			# try edge removal
			new_parents = deleteat!(copy(parents_lat), ind_lat)
			new_score_lat = log_bayes_score_component(ind_lat, new_parents, s_r, s_d, score_cache_lat)
			if new_score_lat - score_lat > best_diff
				selected_lat = true
				best_score_lat = new_score_lat
				best_parents_lat = new_parents
				best_diff = new_score_lat - score_lat
			end

			# try edge reversal (a -> lon)
			new_parents = unshift!(copy(parents_lon), ind_lat)
			# alpha_lon = prior_counts_lon(ind_lon, ind_yaw, ind_d_cl, new_parents, s_r, binmaps, alpha_cache_lon)
			# new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, alpha_lon, score_cache_lon)
			new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, score_cache_lon)
			if new_score_lon - score_lon > best_diff
				selected_lat = false
				best_score_lon = new_score_lon
				best_parents_lon = new_parents
				best_diff = new_score_lon - score_lon
			end
		elseif in(ind_lat, parents_lon)  && !in(features[ind_lat], forced.lon)
			# try edge removal
			new_parents = deleteat!(copy(parents_lon), ind_lat)
			new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, score_cache_lon)
			if new_score_lon - score_lon > best_diff
				selected_lat = false
				best_score_lon = new_score_lon
				best_parents_lon = new_parents
				best_diff = new_score_lon - score_lon
			end

			# try edge reversal (lon -> lat)
			new_parents = unshift!(copy(parents_lat), ind_lon)
			new_score_lat = log_bayes_score_component(ind_lat, new_parents, s_r, s_d, score_cache_lat)
			if new_score_lat - score_lat > best_diff
				selected_lat = true
				best_score_lat = new_score_lat
				best_parents_lat = new_parents
				best_diff = new_score_lat - score_lat
			end
		end

		# select best
		if best_diff > IMPROVEMENT_THRESHOLD
			if selected_lat
				parents_lat = best_parents_lat
				score_lat = best_score_lat
				chosen_lat = map(i->in(n_targets+i, parents_lat), [1:n_indicators])
			else
				parents_lon = best_parents_lon
				score_lon = best_score_lon
				chosen_lon = map(i->in(n_targets+i, parents_lon), [1:n_indicators])
			end
			score = score_lat + score_lon
		else
			done = true
		end
	end

	(parents_lat, parents_lon, score)
end
function learn_random_init{F<:AbstractFeature}(
	features        :: Vector{F},
	n_targets       :: Int,
	n_indicators    :: Int,
	forced          :: ParentFeatures,
	s_r             :: SharedVector{Int},
	s_d             :: SharedMatrix{Int},
	binmaps         :: Dict{Symbol,AbstractBinMap},
	score_cache_lat :: Dict{Vector{Int}, Float64},
	score_cache_lon :: Dict{Vector{Int}, Float64},
	ind_lat         :: Int,
	ind_lon         :: Int;
	verbosity         :: Int = 0,
	n_initializations :: Int = RANDOMINIT_N_INIT,
	n_init_parents    :: Int = RANDOMINIT_N_PARENTS
	)

	# random initialization
	#   - graph traversal with random initialization

	final_parents_lat = Int[]
	final_parents_lon = Int[]
	final_score = -Inf

	while n_initializations > 0
		n_initializations -= 1

		parents_lat, parents_lon = parent_indeces(forced, features)
		parents_lat = sort!(unique(append!(parents_lat, randperm(n_indicators)[1:n_init_parents]+n_targets)))
		parents_lon = sort!(unique(append!(parents_lon, randperm(n_indicators)[1:n_init_parents]+n_targets)))
		chosen_lat  = map(i->in(n_targets+i, parents_lat), [1:n_indicators])
		chosen_lon  = map(i->in(n_targets+i, parents_lon), [1:n_indicators])
		score_lat   = log_bayes_score_component(ind_lat, parents_lat, s_r, s_d, score_cache_lat)
		score_lon   = log_bayes_score_component(ind_lon, parents_lon, s_r, s_d, score_cache_lon)
		score       = score_lat + score_lon

		if verbosity > 1
			@printf("starting score: %15.3f %15.3f %15.3f\n", score_lat, score_lon, score)
		end

		n_iter = 0
		done = false
		while !done
			n_iter += 1

			best_diff = IMPROVEMENT_THRESHOLD
			selected_lat = false
			best_score_lat = 0.0
			best_score_lon = 0.0
			best_parents_lat = parents_lat
			best_parents_lon = parents_lon

			# check edges for indicators -> lat
			for i = 1 : n_indicators
				# add edge if it does not exist
				if !chosen_lat[i]
					new_parents = sort!(push!(copy(parents_lat), n_targets+i))
					new_score_lat = log_bayes_score_component(ind_lat, new_parents, s_r, s_d, score_cache_lat)
					if new_score_lat - score_lat > best_diff
						selected_lat = true
						best_score_lat = new_score_lat
						best_parents_lat = new_parents
						best_diff = new_score_lat - score_lat
					end
				end
			end
			for (idx, i) in enumerate(parents_lat)
				# remove edge if it does exist
				if !in(features[i], forced.lat)
					new_parents = deleteat!(copy(parents_lat), idx)
					new_score_lat = log_bayes_score_component(ind_lat, new_parents, s_r, s_d, score_cache_lat)
					if new_score_lat - score_lat > best_diff
						selected_lat = true
						best_score_lat = new_score_lat
						best_parents_lat = new_parents
						best_diff = new_score_lat - score_lat
					end
				end
			end

			# check edges for indicators -> lon
			for i = 1 : n_indicators
				# add edge if it does not exist
				if !chosen_lon[i]
					new_parents = sort!(push!(copy(parents_lon), n_targets+i))
					new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, score_cache_lon)
					if new_score_lon - score_lon > best_diff
						selected_lat = false
						best_score_lon = new_score_lon
						best_parents_lon = new_parents
						best_diff = new_score_lon - score_lon
					end
				end
			end
			for (idx, i) in enumerate(parents_lon)
				# remove edge if it does exist
				if !in(features[i], forced.lon)
					new_parents = deleteat!(copy(parents_lon), idx)
					new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, score_cache_lon)
					if new_score_lon - score_lon > best_diff
						selected_lat = false
						best_score_lon = new_score_lon
						best_parents_lon = new_parents
						best_diff = new_score_lon - score_lon
					end
				end
			end

			# check edge between lat <-> lon
			if !in(ind_lon, parents_lat) && !in(ind_lat, parents_lon)
				# lon -> lat
				new_parents = unshift!(copy(parents_lat), ind_lon)
				new_score_lat = log_bayes_score_component(ind_lat, new_parents, s_r, s_d, score_cache_lat)
				if new_score_lat - score_lat > best_diff
					selected_lat = true
					best_score_lat = new_score_lat
					best_parents_lat = new_parents
					best_diff = new_score_lat - score_lat
				end

				# lat -> lon
				new_parents = unshift!(copy(parents_lon), ind_lat)
				new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, score_cache_lon)
				if new_score_lon - score_lon > best_diff
					selected_lat = false
					best_score_lon = new_score_lon
					best_parents_lon = new_parents
					best_diff = new_score_lon - score_lon
				end
			elseif in(ind_lon, parents_lat) && !in(features[ind_lon], forced.lat)

				# try edge removal
				new_parents = deleteat!(copy(parents_lat), ind_lat)
				new_score_lat = log_bayes_score_component(ind_lat, new_parents, s_r, s_d, score_cache_lat)
				if new_score_lat - score_lat > best_diff
					selected_lat = true
					best_score_lat = new_score_lat
					best_parents_lat = new_parents
					best_diff = new_score_lat - score_lat
				end

				# try edge reversal (a -> lon)
				new_parents = unshift!(copy(parents_lon), ind_lat)
				# alpha_lon = prior_counts_lon(ind_lon, ind_yaw, ind_d_cl, new_parents, s_r, binmaps, alpha_cache_lon)
				# new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, alpha_lon, score_cache_lon)
				new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, score_cache_lon)
				if new_score_lon - score_lon > best_diff
					selected_lat = false
					best_score_lon = new_score_lon
					best_parents_lon = new_parents
					best_diff = new_score_lon - score_lon
				end
			elseif in(ind_lat, parents_lon)  && !in(features[ind_lat], forced.lon)
				# try edge removal
				new_parents = deleteat!(copy(parents_lon), ind_lat)
				new_score_lon = log_bayes_score_component(ind_lon, new_parents, s_r, s_d, score_cache_lon)
				if new_score_lon - score_lon > best_diff
					selected_lat = false
					best_score_lon = new_score_lon
					best_parents_lon = new_parents
					best_diff = new_score_lon - score_lon
				end

				# try edge reversal (lon -> lat)
				new_parents = unshift!(copy(parents_lat), ind_lon)
				new_score_lat = log_bayes_score_component(ind_lat, new_parents, s_r, s_d, score_cache_lat)
				if new_score_lat - score_lat > best_diff
					selected_lat = true
					best_score_lat = new_score_lat
					best_parents_lat = new_parents
					best_diff = new_score_lat - score_lat
				end
			end

			# select best
			if best_diff > IMPROVEMENT_THRESHOLD
				if selected_lat
					parents_lat = best_parents_lat
					score_lat = best_score_lat
					chosen_lat = map(i->in(n_targets+i, parents_lat), [1:n_indicators])
				else
					parents_lon = best_parents_lon
					score_lon = best_score_lon
					chosen_lon = map(i->in(n_targets+i, parents_lon), [1:n_indicators])
				end
				score = score_lat + score_lon
			else
				done = true
			end
		end

		if verbosity > 0
			println(n_initializations, " score: ", score)
		end

		if score > final_score
			final_score = score
			final_parents_lat = parents_lat
			final_parents_lon = parents_lon
		end
	end

	(final_parents_lat, final_parents_lon, final_score)
end

function graph_feature_selection(
	modelset     :: ModelSet,
	targetset    :: ModelTargets,
	indicatorset :: (String, Vector{AbstractFeature}),
	forced       :: ParentFeatures,
	rawdata      :: DataFrame,
	outputfolder :: String;
	verbosity 				   :: Int  = 0,
	flag_learn_forced          :: Bool = false,
	flag_learn_edge_addition   :: Bool = false,
	flag_learn_graph_traversal :: Bool = false,
	flag_learn_random_init     :: Bool = false,
	pass_forced_to_learn       :: Bool = true
	)

	@assert(flag_learn_forced          || flag_learn_edge_addition ||
		    flag_learn_graph_traversal || flag_learn_random_init)

	target_set_name    = short_name(targetset.lat) * "_" * short_name(targetset.lon)
	indicator_set_name = indicatorset[1]
	indicators         = indicatorset[2]
	set_name           = modelset.name * "_" * target_set_name * "_" * indicator_set_name
	basefolder         = outputfolder * "/" * set_name

	targets  = [targetset.lat, targetset.lon]
	features = [targets, indicators]
	filters  = modelset.filters

	data = deepcopy(rawdata)
	data = filter_features(data, filters)
	data = pull_relevant_features(data, features)
	data = remove_na_target_values(data, targets)

	@assert(size(data, 1) > 0)

	# find target indeces
	ind_lat   = findfirst(f->f==targetset.lat, features)
	ind_lon   = findfirst(f->f==targetset.lon, features)
	@assert(ind_lat  != 0)
	@assert(ind_lon  != 0)

	mat_full  = discretize(binmaps, features, data, Int)

	if DECIMATION > 0
		mat_discr = mat_full[[1:DECIMATION:size(mat_full,1)], :]
	else
		@assert(DECIMATION < 0)
		mult = -DECIMATION
		mat_discr = repmat(mat_full, mult, 1)
		mat_full  = mat_discr
	end

	println("datapoints after decimation: ", size(mat_discr, 1))

	s_d       = SharedArray(Int, (size(mat_discr,2), size(mat_discr,1)))
	s_d[:,:]  = mat_discr'
	s_r       = SharedArray(Int, (size(s_d,1),))
	s_r[:]    = map(f->nlabels(binmaps[symbol(f)]), features)

	# graph search
	##################################################

	n_targets       = length(targets)
	n_indicators    = length(indicators)

	score_cache_lat = Dict{Vector{Int}, Float64}()
	score_cache_lon = Dict{Vector{Int}, Float64}()

	retval = Dict{Symbol, GraphLearningResult}()

	if flag_learn_forced
		parents_lat, parents_lon = parent_indeces(forced, features)
		retval[:forced] = GraphLearningResult(basefolder*"_forced", features, ind_lat, ind_lon, parents_lat,
			                      			  parents_lon, NaN, mat_full, s_r, s_d)
	end

	if !pass_forced_to_learn
		forced = ParentFeatures(AbstractFeature[], AbstractFeature[])
	end

	if flag_learn_edge_addition
		parents_lat, parents_lon, score = learn_edge_addition(features, n_targets, n_indicators, forced,
			                                                  s_r, s_d, binmaps, score_cache_lat, 
			                                                  score_cache_lon, ind_lat, ind_lon,
			                                                  verbosity = verbosity)
		retval[:edge_addition] = GraphLearningResult(basefolder*"_edge_addition", features, 
													 ind_lat, ind_lon, parents_lat,
			                      			  		 parents_lon, score, mat_full, s_r, s_d)
	end

	
	if flag_learn_graph_traversal
		parents_lat, parents_lon, score = learn_graph_traversal(features, n_targets, n_indicators, forced,
			                                                s_r, s_d, binmaps, score_cache_lat, 
			                                                score_cache_lon, ind_lat, ind_lon,
			                                                verbosity = verbosity)
		retval[:graph_traversal] = GraphLearningResult(basefolder*"_graph_traversal", features, 
													   ind_lat, ind_lon, parents_lat,
			                      			  		   parents_lon, score, mat_full, s_r, s_d)
	end

	if flag_learn_random_init
		parents_lat, parents_lon, score = learn_random_init(features, n_targets, n_indicators, forced,
			                                                s_r, s_d, binmaps, score_cache_lat, 
			                                                score_cache_lon, ind_lat, ind_lon,
			                                                verbosity = verbosity)
		retval[:random_init] = GraphLearningResult(basefolder*"_random_init", features, 
												   ind_lat, ind_lon, parents_lat,
			                      			  	   parents_lon, score, mat_full, s_r, s_d)
	end

	retval
end
function display_results(res::GraphLearningResult)

	println("target lat: ", symbol(res.target_lat))
	for f in res.parents_lat
		println("\t\t", symbol(f))
	end
	println("target lon: ", symbol(res.target_lon))
	for f in res.parents_lon
		println("\t\t", symbol(f))
	end

	println("score raw:  ", res.score_raw)
	println("score K2:   ", res.score_K2)
	println("score BDeu: ", res.score_BDeu)
end
function export_results(res::GraphLearningResult, binmaps)

	write_file(res.net, res.fileroot * ".xdsl")
	writetex(res.net, res.fileroot * ".tex", res.features)
	save(res.fileroot * ".jld", 
		"statistics", res.stats, 
		"adjacency",  res.adj,
		"binmaps",    binmaps,
		"targets",    [res.target_lat, res.target_lon],
		"indicators", res.features[3:end])
end

const MODEL_SETS = ModelSet[
	ModelSet("freeflow",  AbstractFeature[SUBSET_FREE_FLOW, SUBSET_AT_SIXTYFIVE]),
	ModelSet("following", AbstractFeature[SUBSET_CAR_FOLLOWING]),
	ModelSet("carfollow", AbstractFeature[]),
	]
const TARGET_SETS = [
	ModelTargets(FUTURETURNRATE_250MS,     FUTUREACCELERATION_250MS),
	ModelTargets(FUTURETURNRATE_500MS,     FUTUREACCELERATION_500MS),
	ModelTargets(FUTURETURNRATE_250MS,     FUTUREDESIREDSPEED_250MS),
	ModelTargets(FUTURETURNRATE_500MS,     FUTUREDESIREDSPEED_500MS),
	ModelTargets(FUTUREDESIREDANGLE_250MS, FUTUREACCELERATION_250MS), # this is what we eneded up using
	ModelTargets(FUTUREDESIREDANGLE_500MS, FUTUREACCELERATION_500MS),
	ModelTargets(FUTUREDESIREDANGLE_250MS, FUTUREDESIREDSPEED_250MS),
	ModelTargets(FUTUREDESIREDANGLE_500MS, FUTUREDESIREDSPEED_500MS),
    ]
const INDICATOR_SETS = [
		("medium", AbstractFeature[
			    YAW, SPEED, VELFX, VELFY, DELTA_SPEED_LIMIT,
		        D_CL, D_ML, D_MR, D_MERGE, D_SPLIT, 
		        TIMETOCROSSING_RIGHT, TIMETOCROSSING_LEFT, TIMESINCELANECROSSING,
		        N_LANE_L, N_LANE_R, HAS_LANE_L, HAS_LANE_R,
		        TURNRATE, TURNRATE_GLOBAL, ACC, ACCFX, ACCFY, A_REQ_STAYINLANE, LANECURVATURE,

		        HAS_FRONT, D_X_FRONT, D_Y_FRONT, V_X_FRONT, V_Y_FRONT, YAW_FRONT, TURNRATE_FRONT,
		        HAS_REAR,  D_X_REAR,  D_Y_REAR,  V_X_REAR,  V_Y_REAR,  YAW_REAR,  TURNRATE_REAR,
		                   D_X_LEFT,  D_Y_LEFT,  V_X_LEFT,  V_Y_LEFT,  YAW_LEFT,  TURNRATE_LEFT,
		                   D_X_RIGHT, D_Y_RIGHT, V_X_RIGHT, V_Y_RIGHT, YAW_RIGHT, TURNRATE_RIGHT,
		        A_REQ_FRONT, TTC_X_FRONT, TIMEGAP_X_FRONT,
		        A_REQ_REAR,  TTC_X_REAR,  TIMEGAP_X_REAR,
		        A_REQ_LEFT,  TTC_X_LEFT,  TIMEGAP_X_LEFT,
		        A_REQ_RIGHT, TTC_X_RIGHT, TIMEGAP_X_RIGHT,

		        SCENEVELFX,

		        TIME_CONSECUTIVE_BRAKE, TIME_CONSECUTIVE_ACCEL,
		             PASTACC250MS,      PASTACC500MS,      PASTACC750MS,      PASTACC1S,
        		PASTTURNRATE250MS, PASTTURNRATE500MS, PASTTURNRATE750MS, PASTTURNRATE1S,
        		   PASTVELFY250MS,    PASTVELFY500MS,    PASTVELFY750MS,    PASTVELFY1S,
        		    PASTD_CL250MS,     PASTD_CL500MS,     PASTD_CL750MS,     PASTD_CL1S,
			]),
		("large", AbstractFeature[
				YAW, SPEED, VELFX, VELFY, DELTA_SPEED_LIMIT,
		        D_CL, D_ML, D_MR, D_MERGE, D_SPLIT, 
		        TIMETOCROSSING_RIGHT, TIMETOCROSSING_LEFT, TIMESINCELANECROSSING,
		        N_LANE_L, N_LANE_R, HAS_LANE_L, HAS_LANE_R,
		        TURNRATE, TURNRATE_GLOBAL, ACC, ACCFX, ACCFY, A_REQ_STAYINLANE, LANECURVATURE,

		        HAS_FRONT, D_X_FRONT, D_Y_FRONT, V_X_FRONT, V_Y_FRONT, YAW_FRONT, TURNRATE_FRONT,
		        HAS_REAR,  D_X_REAR,  D_Y_REAR,  V_X_REAR,  V_Y_REAR,  YAW_REAR,  TURNRATE_REAR,
		                   D_X_LEFT,  D_Y_LEFT,  V_X_LEFT,  V_Y_LEFT,  YAW_LEFT,  TURNRATE_LEFT,
		                   D_X_RIGHT, D_Y_RIGHT, V_X_RIGHT, V_Y_RIGHT, YAW_RIGHT, TURNRATE_RIGHT,
		        A_REQ_FRONT, TTC_X_FRONT, TIMEGAP_X_FRONT,
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
			])
		]
const FORCED = 
		# ParentFeatures(AbstractFeature[D_CL, YAW, VELFY, HAS_LANE_R, HAS_LANE_L], # following
	 #                   AbstractFeature[PASTACC500MS, TIME_CONSECUTIVE_BRAKE, TIME_CONSECUTIVE_ACCEL, D_X_FRONT, V_X_FRONT])
		# ParentFeatures(AbstractFeature[D_CL, YAW, HAS_LANE_R, HAS_LANE_L], # freeflow
	 #                   AbstractFeature[PASTACC500MS, TIME_CONSECUTIVE_BRAKE, TIME_CONSECUTIVE_ACCEL, DELTA_SPEED_LIMIT])
		  # ParentFeatures(AbstractFeature[D_CL], # lanechange
	   #                   AbstractFeature[])
	      # ParentFeatures(AbstractFeature[YAW, VELFY, D_CL, ACCFY, PASTVELFY250MS, TIMESINCELANECROSSING, TIMETOCROSSING_LEFT], # lanechange_forced
	      #                AbstractFeature[ACC, ACCFX, TIME_CONSECUTIVE_BRAKE, TIME_CONSECUTIVE_ACCEL, PASTACC1S])
          ParentFeatures(AbstractFeature[], AbstractFeature[])

outputfolder = create_output_folder("graph_feature_selection_NEW")
reload("path_to/feature_binmaps.jl")
data_orig  = load("path_to_featureset.jld", "data")
println(size(data_orig))

for modelset in MODEL_SETS
	for targetset in TARGET_SETS
		for indicatorset in INDICATOR_SETS

			flagset = "0111"
			res = graph_feature_selection(modelset, targetset, indicatorset, FORCED, data_orig, outputfolder,
										  flag_learn_forced         = (flagset[1] == '1'), 
										  flag_learn_edge_addition  = (flagset[2] == '1'),
										  flag_learn_graph_traversal= (flagset[3] == '1'), 
										  flag_learn_random_init    = (flagset[4] == '1'),
										  pass_forced_to_learn      = true)

			if (flagset[1] == '1')
				println("FORCED")
				display_results(res[:forced])
				export_results(res[:forced], binmaps)
			end

			if (flagset[2] == '1')
				println("EDGE ADDITION")
				display_results(res[:edge_addition])
				export_results(res[:edge_addition], binmaps)
			end

			if (flagset[3] == '1')
				println("GRAPH TRAVERSAL")
				display_results(res[:graph_traversal])
				export_results(res[:graph_traversal], binmaps)
			end

			if (flagset[4] == '1')
				println("RANDOM INIT")
				display_results(res[:random_init])
				export_results(res[:random_init], binmaps)
			end

		end
	end
end