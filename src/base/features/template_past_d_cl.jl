
for ticks in tick_list_short

	counts, unit = ticks_to_time_string(ticks)

	# Past Acceleration
	fname_str = "PastD_CL$counts$unit"
	feature_name = symbol("Feature_" * fname_str)
	str_description  = "the lateral distance covered over past $counts $unit"
	sym_feature = symbol("pastd_cl$counts$unit")
	lstr = latexstring("d_{cl}^\\text{-$counts$unit}")

	create_feature_basics( fname_str, UnextractableFeature, "m", false, false, Inf, -Inf, true, sym_feature, lstr, str_description)

	@eval begin
		function _get(::$feature_name, pdset::PrimaryDataset, sn::StreetNetwork, carind::Int, validfind::Int)

			jump = -$ticks

			carid = carind == CARIND_EGO ? CARID_EGO : carind2id(pdset, carind, validfind)

			jvfind1 = int(jumpframe(pdset, validfind, jump))
			if jvfind1 == 0 # Does not exist
				return NA_ALIAS
			end

			if carind == CARIND_EGO

				posGxA = gete_validfind(pdset, :posGx, validfind)
				posGyA = gete_validfind(pdset, :posGy, validfind)
				posGxB = gete_validfind(pdset, :posGx, jvfind1)
				posGyB = gete_validfind(pdset, :posGy, jvfind1)

				Δs, Δd = frenet_distance_between_points(sn, posGxA, posGyA, posGxB, posGyB)

				if !isnan(Δd)
					return Δd
				else
					return NA_ALIAS
				end

			elseif idinframe(pdset, carid, jvfind1) && idinframe(pdset, carid, jvfind2)
				ind = carid2ind(pdset, carid, jvfind1)

				posGxA = get(pdset, "posGx", carind, validfind)
				posGyA = get(pdset, "posGy", carind, validfind)
				posGxB = get(pdset, "posGx", carind, jvfind1)
				posGyB = get(pdset, "posGy", carind, jvfind1)

				Δs, Δd = frenet_distance_between_points(sn, posGxA, posGyA, posGxB, posGyB)

				if !isnan(Δd)
					return Δd
				else
					return NA_ALIAS
				end
			end

			return NA_ALIAS
		end
	end


end

