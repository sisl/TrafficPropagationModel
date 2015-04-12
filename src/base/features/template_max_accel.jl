
for ticks in tick_list

	counts, unit = ticks_to_time_string(ticks)

	# Max Acceleration F x
	fname_str = "MaxAccFx$counts$unit"
	feature_name = symbol("Feature_" * fname_str)
	str_description  = "the maximum acceleration along the lane for the last $counts $unit"
	sym_feature = symbol("maxaccFx$counts$unit")
	lstr = latexstring("\\hat{a}^F_{x,$counts$unit}")

	create_feature_basics( fname_str, UnextractableFeature, "m/s2", false, false, Inf, -Inf, true, sym_feature, lstr, str_description)


	@eval begin
		function _get(::$feature_name, pdset::PrimaryDataset, carind::Int, validfind::Int)

			maxjump = -$ticks

			val = 0.0
			count = 0

			carid = carind == CARIND_EGO ? CARID_EGO : carind2id(pdset, carind, validfind)

			for jump = -1:-1:maxjump

				jvfind1 = jumpframe(pdset, validfind, jump)
				jvfind2 = jumpframe(pdset, validfind, jump+1)
				if jvfind1 == 0 || jvfind2 == 0 # Does not exist
					continue
				end

				cur, fut = 0.0, 0.0
				if carind == CARIND_EGO
					f1 = validfind2frameind(pdset, jvfind1)
					f2 = validfind2frameind(pdset, jvfind2)
					cur = gete(pdset, :velFx, f1)
					fut = gete(pdset, :velFx, f2)
				elseif idinframe(pdset, carid, jvfind1) && idinframe(pdset, carid, jvfind2)
					ind1 = carid2ind(pdset, carid, jvfind1)
					ind2 = carid2ind(pdset, carid, jvfind2)
					cur = getc(pdset, "velFx", ind1, jvfind1)
					fut = getc(pdset, "velFx", ind2, jvfind2)	
				else
					continue
				end

				count += 1
				acc = (fut - cur)/SEC_PER_FRAME
				val = abs(acc) > abs(val) ? acc : val
			end

			if count == 0
				return NA_ALIAS
			end
			val
		end
	end

	# Max Acceleration F x
	fname_str = "MaxAccFy$counts$unit"
	feature_name = symbol("Feature_" * fname_str)
	str_description  = "the maximum acceleration perpendicular to the lane for the last $counts $unit"
	sym_feature = symbol("maxaccFy$counts$unit")
	lstr = latexstring("\\hat{a}^F_{y,$counts$unit}")

	create_feature_basics( fname_str, UnextractableFeature, "m/s2", false, false, Inf, -Inf, true, sym_feature, lstr, str_description)

	@eval begin
		function _get(::$feature_name, pdset::PrimaryDataset, carind::Int, validfind::Int)

			maxjump = -$ticks

			val = 0.0
			count = 0

			carid = carind == CARIND_EGO ? CARID_EGO : carind2id(pdset, carind, validfind)

			for jump = -1:-1:maxjump

				jvfind1 = jumpframe(pdset, validfind, jump)
				jvfind2 = jumpframe(pdset, validfind, jump+1)
				if jvfind1 == 0 || jvfind2 == 0 # Does not exist
					continue
				end

				cur, fut = 0.0, 0.0
				if carind == CARIND_EGO
					f1 = validfind2frameind(pdset, jvfind1)
					f2 = validfind2frameind(pdset, jvfind2)
					cur = gete(pdset, :velFy, f1)
					fut = gete(pdset, :velFy, f2)
				elseif idinframe(pdset, carid, jvfind1) && idinframe(pdset, carid, jvfind2)
					ind1 = carid2ind(pdset, carid, jvfind1)
					ind2 = carid2ind(pdset, carid, jvfind2)
					cur = getc(pdset, "velFy", ind1, jvfind1)
					fut = getc(pdset, "velFy", ind2, jvfind2)	
				else
					continue
				end

				count += 1
				acc = (fut - cur)/SEC_PER_FRAME
				val = abs(acc) > abs(val) ? acc : val
			end

			if count == 0
				return NA_ALIAS
			end
			val
		end
	end

end

