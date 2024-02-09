# This script was written by Jin Xu and available on Github
# https://github.com/SunnyXu/artificial_random_signaling_network

cd(dirname(@__FILE__))
#include("rr_funcs-Jin.jl")
using Random
using RoadRunner
using StatsBase # random pick reaction by weight
#using Suppressor
RoadRunner.disableLoggingToConsole() # try to disable some warnings like NLEQ

# number of floating species
global nSpecies = 25
# this number excludes the input and output species but includes nSpecies_gene
#Therefore, the total number of species should be (nSpecies+2)

# maximum number of reactions
global nRxns = 35

# number of gene species
global nSpecies_gene = 15
# some restrictions:s
# (nSpecies-nSpecies_gene) >= 4 for the case of BIBI
# nSpecies >= 5 for the case of double catalyzation
# nSpecies >= 13: 6 to initiate input and output reactions + next input and output layer 7

# random values assignment
# random number generation for species concentration [0,10)
rnd_species = 10.
# random number generation for parameters [0,1)
rnd_parameter = 1.
# doubling the concentration at the input species.
concentration_perturb = 2.
# number of sampleNetwork.xml to generate
sampleSize = 1
# i.e. sampleNetwork-1.xml, sampleNetwork-2.xml

RoadRunner.setConfigBool("ROADRUNNER_DISABLE_WARNINGS", 1)
function rv_specs(ids_species, ids_rv)
    for i = 1:size(ids_rv)[1]
        if ids_rv[i] in ids_species
            filter!(e->e≠ids_rv[i], ids_species)
        end
    end
    return ids_species
end

function randomNetwork(nSpecies, nSpecies_gene, nRxns)
    INPUT_RXN_MECH = ["UNICAT", "UNIBI", "BIUNI", "BIBI"]
    INPUT_RXN_MECH_WEIGHT = [0.25, 0.25, 0.25, 0.25]
    RXN_MECH = ["UNICAT", "UNIBI", "BIUNI", "BIBI", "CIRCLE", "DBCIRCLE"]
    RXN_MECH_WEIGHT = [0.2, 0.2, 0.2, 0.2, 0.1, 0.1]
    GN_RXN_MECH = ["UNICAT", "CIRCLE", "DBCIRCLE"]
    GN_RXN_MECH_WEIGHT = [0.25, 0.33, 0.42]

    species = ["S$i" for i = 1:nSpecies]
    gene_species = species[1:nSpecies_gene]
    rct_counter = zeros(nSpecies)
    prd_counter = zeros(nSpecies)

    global num_common_specs = 0
    while num_common_specs == 0
        global rxn_counter = 0
        global rxn_specs = Dict{String, Any}()
        ids_tot = [] #to check if all the species are involved in reactions
        global species_ids = collect(1:nSpecies)
        global non_gene_species_ids = collect((nSpecies_gene+1):nSpecies)

        rxn_counter += 1
        rxn_mechanism = sample(INPUT_RXN_MECH, Weights(INPUT_RXN_MECH_WEIGHT))
        if rxn_mechanism == "UNICAT"
            parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter", "K1_J$rxn_counter", "K2_J$rxn_counter"]
            ids = sample(species_ids, 2, replace = false)
            prds_id = ids[1]
            cats_id = ids[2]
            rv_ids = []
            append!(rv_ids, prds_id)
            append!(rv_ids, cats_id)
            species_ids = rv_specs(species_ids, rv_ids)
            non_gene_species_ids = rv_specs(non_gene_species_ids, rv_ids)
            rcts = ["S_in"]
            prds = ["S$prds_id"]
            cats = ["S$cats_id"]
            rxn_mech = ["S$cats_id * (kf_J$rxn_counter*S_in/K1_J$rxn_counter - kr_J$rxn_counter*S$prds_id/K2_J$rxn_counter) /(1 + S_in/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
            rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
        elseif rxn_mechanism == "UNIBI"
            parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter"]
            ids = sample(non_gene_species_ids, 2, replace = false)
            prds_id  = ids[1]
            prds_id2 = ids[2]
            rv_ids = []
            append!(rv_ids, prds_id)
            append!(rv_ids, prds_id2)
            species_ids = rv_specs(species_ids, rv_ids)
            non_gene_species_ids = rv_specs(non_gene_species_ids, rv_ids)
            rcts = ["S_in"]
            prds = ["S$prds_id", "S$prds_id2"]
            cats = []
            rxn_mech = ["kf_J$rxn_counter * S_in - kr_J$rxn_counter * S$prds_id * S$prds_id2"]
            rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
        elseif rxn_mechanism == "BIUNI"
            parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter"]
            ids = sample(non_gene_species_ids, 2, replace = false)
            rcts_id = ids[1]
            prds_id  = ids[2]
            rv_ids = []
            append!(rv_ids, rcts_id)
            append!(rv_ids, prds_id)
            species_ids = rv_specs(species_ids, rv_ids)
            non_gene_species_ids = rv_specs(non_gene_species_ids, rv_ids)
            rcts = ["S_in", "S$rcts_id"]
            prds = ["S$prds_id"]
            cats = []
            rxn_mech = ["kf_J$rxn_counter * S_in * S$rcts_id - kr_J$rxn_counter * S$prds_id"]
            rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
        else rxn_mechanism == "BIBI"
            parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter"]
            ids = sample(non_gene_species_ids, 3, replace = false)
            rcts_id  = ids[1]
            prds_id  = ids[2]
            prds_id2 = ids[3]
            rv_ids = []
            append!(rv_ids, rcts_id)
            append!(rv_ids, prds_id)
            append!(rv_ids, prds_id2)
            species_ids = rv_specs(species_ids, rv_ids)
            non_gene_species_ids = rv_specs(non_gene_species_ids, rv_ids)
            rcts = ["S_in", "S$rcts_id"]
            prds = ["S$prds_id", "S$prds_id2"]
            cats = []
            rxn_mech = ["kf_J$rxn_counter * S_in * S$rcts_id - kr_J$rxn_counter * S$prds_id * S$prds_id2"]
            rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
        end

        append!(ids_tot, ids)

        rxn_counter += 1

        rxn_mechanism == "CIRCLE"
        parameters1 = ["kf_J$rxn_counter", "kr_J$rxn_counter", "K1_J$rxn_counter", "K2_J$rxn_counter"]
        parameters2 = ["kf_J$(rxn_counter+1)", "kr_J$(rxn_counter+1)", "K1_J$(rxn_counter+1)", "K2_J$(rxn_counter+1)"]
        ids = sample(species_ids, 3, replace = false)
        rcts_id = ids[1]
        cats_id = ids[2]
        cats_id2 = ids[3]
        rv_ids = []
        append!(rv_ids, rcts_id)
        append!(rv_ids, cats_id)
        append!(rv_ids, cats_id2)
        species_ids = rv_specs(species_ids, rv_ids)
        non_gene_species_ids = rv_specs(non_gene_species_ids, rv_ids)
        rcts = ["S$rcts_id"]
        prds = ["S_out"]
        cats1 = ["S$cats_id"]
        cats2 = ["S$cats_id2"]
        rxn_mech1 = ["S$cats_id * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter - kr_J$rxn_counter*S_out/K2_J$rxn_counter) /(1 + S$rcts_id/K1_J$rxn_counter + S_out/K2_J$rxn_counter)"]
        rxn_mech2 = ["S$cats_id2 * (kf_J$(rxn_counter+1)*S_out/K2_J$(rxn_counter+1) - kr_J$(rxn_counter+1)*S$rcts_id/K1_J$(rxn_counter+1)) /(1 + S_out/K2_J$(rxn_counter+1) + S$rcts_id/K1_J$(rxn_counter+1))"]
        rxn_specs["r$rxn_counter"]     = Dict{String, Array{String, 1}}("parameters" => parameters1, "rcts" => rcts, "prds" => prds, "cats" => cats1, "rxn_mech" => rxn_mech1)
        rxn_specs["r$(rxn_counter+1)"] = Dict{String, Array{String, 1}}("parameters" => parameters2, "rcts" => prds, "prds" => rcts, "cats" => cats2, "rxn_mech" => rxn_mech2)
        rxn_counter += 1

        append!(ids_tot, ids)

        specs_input = [] #species connected to input layer
        specs_output = [] #species connected to output layer
        # input layer prds can only be next layer's input or catalyzation
        # for CIRCLE and DBCIRCLE only can be the catalyzation
        # note that the cat was not considered here

        # output layer rcts and cats can be the closest upper layer's prds only
        rcts = rxn_specs["r2"]["rcts"]
        append!(specs_output, rcts)
        cats = rxn_specs["r2"]["cats"]
        append!(specs_output, cats)
        cats = rxn_specs["r3"]["cats"]
        append!(specs_output, cats)
        specs_output_selec = specs_output[rand(1:size(specs_output)[1])]

        species_ids_out = []
        append!(species_ids_out, species_ids)
        non_gene_species_ids_out = []
        append!(non_gene_species_ids_out, non_gene_species_ids)
        specs_output = rv_specs(specs_output, [specs_output_selec])
        if (size(specs_output)[1] > 0)
            for i = 1:nSpecies
                for j = 1:size(specs_output)[1]
                    if "S$i" == specs_output[j]
                        if "S$i" in gene_species
                            append!(species_ids_out, i)
                        else
                            append!(non_gene_species_ids_out, i)
                            append!(species_ids_out, i)
                        end
                    end
                end
            end
        end

        rxn_counter += 1
        if specs_output_selec in gene_species #rxn_mechanism == "UNICAT"
            parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter", "K1_J$rxn_counter", "K2_J$rxn_counter"]
            ids = sample(species_ids, 2, replace = false)
            ids_out = []
            rcts_id = ids[1]
            cats_id = ids[2]
            rcts = ["S$rcts_id"]
            prds = [specs_output_selec]
            cats = ["S$cats_id"]
            rxn_mech = ["S$cats_id * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter - kr_J$rxn_counter*$specs_output_selec/K2_J$rxn_counter) /(1 + S$rcts_id/K1_J$rxn_counter + $specs_output_selec/K2_J$rxn_counter)"]
            rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
        else
            rxn_mechanism = sample(INPUT_RXN_MECH, Weights(INPUT_RXN_MECH_WEIGHT))
            if rxn_mechanism == "UNICAT"
                parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter", "K1_J$rxn_counter", "K2_J$rxn_counter"]
                ids = sample(species_ids, 2, replace = false)
                ids_out = []
                rcts_id = ids[1]
                cats_id = ids[2]
                rcts = ["S$rcts_id"]
                prds = [specs_output_selec]
                cats = ["S$cats_id"]
                rxn_mech = ["S$cats_id * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter - kr_J$rxn_counter*$specs_output_selec/K2_J$rxn_counter) /(1 + S$rcts_id/K1_J$rxn_counter + $specs_output_selec/K2_J$rxn_counter)"]
                rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
            elseif rxn_mechanism == "UNIBI"
                parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter"]
                ids = sample(non_gene_species_ids, 1, replace = false)
                non_gene_species_ids_out = rv_specs(non_gene_species_ids_out, ids)
                ids_out = sample(non_gene_species_ids_out, 1, replace = false)
                append!(non_gene_species_ids_out, ids)
                if rand() < 0.5
                    rcts_id  = ids[1]
                    prds_id2 = ids_out[1]
                    rcts = ["S$rcts_id"]
                    prds = [specs_output_selec, "S$prds_id2"]
                    cats = []
                    rxn_mech = ["kf_J$rxn_counter * S$rcts_id - kr_J$rxn_counter * $specs_output_selec * S$prds_id2"]
                else
                    rcts_id  = ids[1]
                    prds_id  = ids_out[1]
                    rcts = ["S$rcts_id"]
                    prds = ["S$prds_id", specs_output_selec]
                    cats = []
                    rxn_mech = ["kf_J$rxn_counter * S$rcts_id - kr_J$rxn_counter * S$prds_id * $specs_output_selec"]
                end
                rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
            elseif rxn_mechanism == "BIUNI"
                parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter"]
                ids = sample(non_gene_species_ids, 2, replace = false)
                ids_out = []
                rcts_id  = ids[1]
                rcts_id2 = ids[2]
                rcts = ["S$rcts_id", "S$rcts_id2"]
                prds = [specs_output_selec]
                cats = []
                rxn_mech = ["kf_J$rxn_counter * S$rcts_id * S$rcts_id2 - kr_J$rxn_counter * $specs_output_selec"]
                rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
            else rxn_mechanism == "BIBI"
                parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter"]
                ids = sample(non_gene_species_ids, 2, replace = false)
                non_gene_species_ids_out = rv_specs(non_gene_species_ids_out, ids)
                ids_out = sample(non_gene_species_ids_out, 1, replace = false)
                append!(non_gene_species_ids_out, ids)
                if rand() < 0.5
                    rcts_id  = ids[1]
                    rcts_id2 = ids[2]
                    prds_id2 = ids_out[1]
                    rcts = ["S$rcts_id", "S$rcts_id2"]
                    prds = [specs_output_selec, "S$prds_id2"]
                    cats = []
                    rxn_mech = ["kf_J$rxn_counter * S$rcts_id * S$rcts_id2 - kr_J$rxn_counter * $specs_output_selec * S$prds_id2"]
                else
                    rcts_id  = ids[1]
                    rcts_id2 = ids[2]
                    prds_id  = ids_out[1]
                    rcts = ["S$rcts_id", "S$rcts_id2"]
                    prds = ["S$prds_id", specs_output_selec]
                    cats = []
                    rxn_mech = ["kf_J$rxn_counter * S$rcts_id * S$rcts_id2 - kr_J$rxn_counter * S$prds_id * $specs_output_selec"]
                end
                rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
            end
        end
        append!(ids_tot, ids)
        append!(ids_tot, ids_out)

        prds = rxn_specs["r1"]["prds"]
        append!(specs_input, prds)

#from the layer next to the input layout, get started with random connected networks

        all_species_involved = false
        # check if all species are involved
        nSpecies_involved = 0
        for i=1:nSpecies
            if i in ids_tot
                nSpecies_involved += 1
            end
        end
        if nSpecies_involved == nSpecies
            global all_species_involved = true
        end

        valid_nRxns = false
        if rxn_counter == nRxns
            valid_nRxns = true
        end

        while all_species_involved == false || valid_nRxns == false
            # check if all species are involved and fit the nRxns
            nSpecies_involved = 0
            for i=1:nSpecies
                if i in ids_tot
                    nSpecies_involved += 1
                end
            end
            if nSpecies_involved == nSpecies
                global all_species_involved = true
            end

            # at least one of the specs_input should be connected to the next layer
			# random select one species from the prds
			specs_input_selec = specs_input[rand(1:size(specs_input)[1])]
			for i = 1: nSpecies
				if "S$i" == specs_input_selec
					global specs_input_selec_ids = [i]
				end
			end

			species_ids_in = []
			append!(species_ids_in, species_ids)
			non_gene_species_ids_in = []
			append!(non_gene_species_ids_in, non_gene_species_ids)
			specs_input = rv_specs(specs_input, [specs_input_selec])
			if (size(specs_input)[1] > 0)
				for i = 1:nSpecies
					for j = 1:size(specs_input)[1]
						if "S$i" == specs_input[j]
							if "S$i" in gene_species
								append!(species_ids_in, i)
							else
								append!(non_gene_species_ids_in, i)
								append!(species_ids_in, i)
							end
						end
					end
				end
			end

	        rxn_counter_start = rxn_counter
            global rxn_counter += 1

		    if specs_input_selec in gene_species
				rxn_mechanism = sample(GN_RXN_MECH, Weights(GN_RXN_MECH_WEIGHT))
				if rxn_mechanism == "UNICAT"
					parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter", "K1_J$rxn_counter", "K2_J$rxn_counter"]
					try
						ids = sample(species_ids, 1, replace = false)
						ids_in = sample(species_ids_in, 1, replace = false)
						append!(species_ids_in, ids)
					catch e
						temp_ids = collect(1:nSpecies)
						ids = sample(temp_ids, 1, replace = false)
						species_ids_in = rv_specs(temp_ids, ids)
						ids_in = sample(species_ids_in, 1, replace = false)
					end
					if rand() < 0.5
						prds_id = ids[1]
						cats_id = ids_in[1]
						rcts = [specs_input_selec]
						prds = ["S$prds_id"]
						cats = ["S$cats_id"]
						rxn_mech = ["S$cats_id * (kf_J$rxn_counter*$specs_input_selec/K1_J$rxn_counter - kr_J$rxn_counter * S$prds_id/K2_J$rxn_counter)/(1 + $specs_input_selec/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
					else
						rcts_id = ids_in[1]
						prds_id = ids[1]
						rcts = ["S$rcts_id"]
						prds = ["S$prds_id"]
						cats = [specs_input_selec]
						rxn_mech = ["$specs_input_selec * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter - kr_J$rxn_counter*S$prds_id/K2_J$rxn_counter) /(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
					end
					rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
				elseif rxn_mechanism == "CIRCLE"
					parameters1 = ["kf_J$rxn_counter", "kr_J$rxn_counter", "K1_J$rxn_counter", "K2_J$rxn_counter"]
					parameters2 = ["kf_J$(rxn_counter+1)", "kr_J$(rxn_counter+1)", "K1_J$(rxn_counter+1)", "K2_J$(rxn_counter+1)"]
					try
						ids = sample(species_ids, 2, replace = false)
						ids_in = sample(species_ids_in, 1, replace = false)
						append!(species_ids_in, ids)
					catch
						temp_ids = collect(1:nSpecies)
						ids = sample(temp_ids, 2, replace = false)
						ids_in = sample(species_ids_in, 1, replace = false)
					end
					if rand() < 0.5
						rcts_id = ids[1]
						prds_id = ids[2]
						cats_id2 = ids_in[1]
						rcts = ["S$rcts_id"]
						prds = ["S$prds_id"]
						cats1 = [specs_input_selec]
						cats2 = ["S$cats_id2"]
						rxn_mech1 = ["$specs_input_selec * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter - kr_J$rxn_counter*S$prds_id/K2_J$rxn_counter) /(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
						rxn_mech2 = ["S$cats_id2* (kf_J$(rxn_counter+1)*S$prds_id/K2_J$(rxn_counter+1) - kr_J$(rxn_counter+1)*S$rcts_id/K1_J$(rxn_counter+1)) /(1 + S$prds_id/K2_J$(rxn_counter+1) + S$rcts_id/K1_J$(rxn_counter+1))"]
					else
						rcts_id = ids[1]
						prds_id = ids[2]
						cats_id = ids_in[1]
						rcts = ["S$rcts_id"]
						prds = ["S$prds_id"]
						cats1 = ["S$cats_id"]
						cats2 = [specs_input_selec]
						rxn_mech1 = ["S$cats_id * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter - kr_J$rxn_counter*S$prds_id/K2_J$rxn_counter) /(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
						rxn_mech2 = ["$specs_input_selec * (kf_J$(rxn_counter+1)*S$prds_id/K2_J$(rxn_counter+1) - kr_J$(rxn_counter+1)*S$rcts_id/K1_J$(rxn_counter+1)) /(1 + S$prds_id/K2_J$(rxn_counter+1) + S$rcts_id/K1_J$(rxn_counter+1))"]
					end
					rxn_specs["r$rxn_counter"]     = Dict{String, Array{String, 1}}("parameters" => parameters1, "rcts" => rcts, "prds" => prds, "cats" => cats1, "rxn_mech" => rxn_mech1)
					rxn_specs["r$(rxn_counter+1)"] = Dict{String, Array{String, 1}}("parameters" => parameters2, "rcts" => prds, "prds" => rcts, "cats" => cats2, "rxn_mech" => rxn_mech2)
					rxn_counter += 1
				else rxn_mechanism == "DBCIRCLE"
                    # parameters1 = ["kf_J$rxn_counter", "K1_J$rxn_counter", "K2_J$rxn_counter"]
                    # parameters2 = ["kf_J$(rxn_counter+1)"]
                    # parameters3 = ["kf_J$(rxn_counter+2)", "K3_J$(rxn_counter+2)", "K4_J$(rxn_counter+2)", "K5_J$(rxn_counter+2)"]
                    # parameters4 = ["kf_J$(rxn_counter+3)"]
					# try
					# 	ids = sample(species_ids, 3, replace = false)
					# 	ids_in = sample(species_ids_in, 1, replace = false)
					# 	append!(species_ids_in, ids)
					# catch
					# 	temp_ids = collect(1:nSpecies)
					# 	ids = sample(temp_ids, 3, replace = false)
					# 	ids_in = sample(species_ids_in, 1, replace = false)
					# end
					# if rand() < 0.5
					# 	rcts_id  = ids[1]
					# 	prds_id  = ids[2]
					# 	prds_id2 = ids[3]
					# 	cats_id2 = ids_in[1]
					# 	rcts1 = ["S$rcts_id"]
					# 	prds1 = ["S$prds_id"]
					# 	prds2 = ["S$prds_id2"]
					# 	cats1 = [specs_input_selec]
					# 	cats2 = ["S$cats_id2"]
                    #     rxn_mech1 = ["$specs_input_selec  * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter)/(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
                    #     rxn_mech2 = ["$specs_input_selec  * (kf_J$(rxn_counter+1)*S$prds_id/K2_J$rxn_counter)/(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
                    #     rxn_mech3 = ["S$cats_id2 * (kf_J$(rxn_counter+2)*S$prds_id2/K3_J$(rxn_counter+2))/(1 + S$prds_id2/K3_J$(rxn_counter+2) + S$prds_id/K4_J$(rxn_counter+2) + S$rcts_id/K5_J$(rxn_counter+2))"]
                    #     rxn_mech4 = ["S$cats_id2 * (kf_J$(rxn_counter+3)*S$prds_id/K4_J$(rxn_counter+2))/(1 + S$prds_id2/K3_J$(rxn_counter+2) + S$prds_id/K4_J$(rxn_counter+2) + S$rcts_id/K5_J$(rxn_counter+2))"]
					# else
					# 	rcts_id  = ids[1]
					# 	prds_id  = ids[2]
					# 	prds_id2 = ids[3]
					# 	cats_id  = ids_in[1]
					# 	rcts1 = ["S$rcts_id"]
					# 	prds1 = ["S$prds_id"]
					# 	prds2 = ["S$prds_id2"]
					# 	cats1 = ["S$cats_id"]
					# 	cats2 = [specs_input_selec]
                    #     rxn_mech1 = ["S$cats_id * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter)/(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
                    #     rxn_mech2 = ["S$cats_id * (kf_J$(rxn_counter+1)*S$prds_id/K2_J$rxn_counter)/(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
                    #     rxn_mech3 = ["$specs_input_selec * (kf_J$(rxn_counter+2)*S$prds_id2/K3_J$(rxn_counter+2))/(1 + S$prds_id2/K3_J$(rxn_counter+2) + S$prds_id/K4_J$(rxn_counter+2) + S$rcts_id/K5_J$(rxn_counter+2))"]
                    #     rxn_mech4 = ["$specs_input_selec * (kf_J$(rxn_counter+3)*S$prds_id/K4_J$(rxn_counter+2))/(1 + S$prds_id2/K3_J$(rxn_counter+2) + S$prds_id/K4_J$(rxn_counter+2) + S$rcts_id/K5_J$(rxn_counter+2))"]
					# end
                    parameters1 = ["kf_J$rxn_counter", "kr_J$rxn_counter", "K1_J$rxn_counter", "K2_J$rxn_counter"]
                    parameters2 = ["kf_J$(rxn_counter+1)", "kr_J$(rxn_counter+1)", "K1_J$(rxn_counter+1)", "K2_J$(rxn_counter+1)"]
                    parameters3 = ["kf_J$(rxn_counter+2)", "kr_J$(rxn_counter+2)", "K1_J$(rxn_counter+2)", "K2_J$(rxn_counter+2)"]
                    parameters4 = ["kf_J$(rxn_counter+3)", "kr_J$(rxn_counter+3)", "K1_J$(rxn_counter+3)", "K2_J$(rxn_counter+3)"]
                    try
                        ids = sample(species_ids, 3, replace = false)
                        species_ids_in = rv_specs(species_ids_in, ids)
                        ids_in = sample(species_ids_in, 1, replace = false)
                        append!(species_ids_in, ids)
                    catch
                        temp_ids = collect(1:nSpecies)
                        temp_ids = rv_specs(temp_ids, specs_input_selec_ids)
                        ids = sample(temp_ids, 3, replace = false)
                        species_ids_in = rv_specs(temp_ids, ids)
                        ids_in = sample(species_ids_in, 1, replace = false)
                    end
                    if rand() < 0.5
                        rcts_id  = ids[1]
                        prds_id  = ids[2]
                        prds_id2 = ids[3]
                        cats_id2 = ids_in[1]
                        rv_ids = []
                        append!(rv_ids, rcts_id)
                        append!(rv_ids, prds_id)
                        append!(rv_ids, prds_id2)
                        append!(rv_ids, cats_id2)
                        species_ids = rv_specs(species_ids, rv_ids)
                        non_gene_species_ids = rv_specs(non_gene_species_ids, rv_ids)
                        rcts1 = ["S$rcts_id"]
                        prds1 = ["S$prds_id"]
                        prds2 = ["S$prds_id2"]
                        cats1 = [specs_input_selec]
                        cats2 = ["S$cats_id2"]
                        rxn_mech1 = ["$specs_input_selec  * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter - kr_J$rxn_counter*S$prds_id/K2_J$rxn_counter) /(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
                        rxn_mech2 = ["$specs_input_selec  * (kf_J$(rxn_counter+1)*S$prds_id/K1_J$(rxn_counter+1) - kr_J$(rxn_counter+1)*S$prds_id2/K2_J$(rxn_counter+1)) /(1 + S$prds_id/K1_J$(rxn_counter+1) + S$prds_id2/K2_J$(rxn_counter+1))"]
                        rxn_mech3 = ["S$cats_id2*(kf_J$(rxn_counter+2)*S$prds_id2/K1_J$(rxn_counter+2) - kr_J$(rxn_counter+2)*S$prds_id/K2_J$(rxn_counter+2)) /(1 + S$prds_id2/K1_J$(rxn_counter+2) + S$prds_id/K2_J$(rxn_counter+2))"]
                        rxn_mech4 = ["S$cats_id2 * (kf_J$(rxn_counter+3)*S$prds_id/K1_J$(rxn_counter+3) - kr_J$(rxn_counter+3)*S$rcts_id/K2_J$(rxn_counter+3)) /(1 + S$prds_id/K1_J$(rxn_counter+3) + S$rcts_id/K2_J$(rxn_counter+3))"]
                    else
                        rcts_id  = ids[1]
                        prds_id  = ids[2]
                        prds_id2 = ids[3]
                        cats_id  = ids_in[1]
                        rv_ids = []
                        append!(rv_ids, rcts_id)
                        append!(rv_ids, prds_id)
                        append!(rv_ids, prds_id2)
                        append!(rv_ids, cats_id)
                        species_ids = rv_specs(species_ids, rv_ids)
                        non_gene_species_ids = rv_specs(non_gene_species_ids, rv_ids)
                        rcts1 = ["S$rcts_id"]
                        prds1 = ["S$prds_id"]
                        prds2 = ["S$prds_id2"]
                        cats1 = ["S$cats_id"]
                        cats2 = [specs_input_selec]
                        rxn_mech1 = ["S$cats_id  * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter - kr_J$rxn_counter*S$prds_id/K2_J$rxn_counter) /(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
                        rxn_mech2 = ["S$cats_id  * (kf_J$(rxn_counter+1)*S$prds_id/K1_J$(rxn_counter+1) - kr_J$(rxn_counter+1)*S$prds_id2/K2_J$(rxn_counter+1)) /(1 + S$prds_id/K1_J$(rxn_counter+1) + S$prds_id2/K2_J$(rxn_counter+1))"]
                        rxn_mech3 = ["$specs_input_selec * (kf_J$(rxn_counter+2)*S$prds_id2/K1_J$(rxn_counter+2) - kr_J$(rxn_counter+2)*S$prds_id/K2_J$(rxn_counter+2)) /(1 + S$prds_id2/K1_J$(rxn_counter+2) + S$prds_id/K2_J$(rxn_counter+2))"]
                        rxn_mech4 = ["$specs_input_selec * (kf_J$(rxn_counter+3)*S$prds_id/K1_J$(rxn_counter+3) - kr_J$(rxn_counter+3)*S$rcts_id/K2_J$(rxn_counter+3)) /(1 + S$prds_id/K1_J$(rxn_counter+3) + S$rcts_id/K2_J$(rxn_counter+3))"]
                    end
					rxn_specs["r$rxn_counter"]     = Dict{String, Array{String, 1}}("parameters" => parameters1, "rcts" => rcts1, "prds" => prds1, "cats" => cats1, "rxn_mech" => rxn_mech1)
					rxn_specs["r$(rxn_counter+1)"] = Dict{String, Array{String, 1}}("parameters" => parameters2, "rcts" => prds1, "prds" => prds2, "cats" => cats1, "rxn_mech" => rxn_mech2)
					rxn_specs["r$(rxn_counter+2)"] = Dict{String, Array{String, 1}}("parameters" => parameters3, "rcts" => prds2, "prds" => prds1, "cats" => cats2, "rxn_mech" => rxn_mech3)
					rxn_specs["r$(rxn_counter+3)"] = Dict{String, Array{String, 1}}("parameters" => parameters4, "rcts" => prds1, "prds" => rcts1, "cats" => cats2, "rxn_mech" => rxn_mech4)
					rxn_counter += 3
				end

			else # specs_input_selec not in gene_species
				rxn_mechanism = sample(RXN_MECH, Weights(RXN_MECH_WEIGHT))
				if rxn_mechanism == "UNICAT"
					parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter", "K1_J$rxn_counter", "K2_J$rxn_counter"]
					try
						ids = sample(species_ids, 1, replace = false)
						ids_in = sample(species_ids_in, 1, replace = false)
						append!(species_ids_in, ids)
					catch
						temp_ids = collect(1:nSpecies)
						ids = sample(temp_ids, 1, replace = false)
						ids_in = sample(species_ids_in, 1, replace = false)
					end

					if rand() < 0.5
						prds_id = ids[1]
						cats_id = ids_in[1]
						rcts = [specs_input_selec]
						prds = ["S$prds_id"]
						cats = ["S$cats_id"]
						rxn_mech = ["S$cats_id * (kf_J$rxn_counter*$specs_input_selec/K1_J$rxn_counter - kr_J$rxn_counter*S$prds_id/K2_J$rxn_counter) /(1 + $specs_input_selec/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
					else
						rcts_id = ids_in[1]
						prds_id = ids[1]
						rcts = ["S$rcts_id"]
						prds = ["S$prds_id"]
						cats = [specs_input_selec]
						rxn_mech = ["$specs_input_selec * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter - kr_J$rxn_counter*S$prds_id/K2_J$rxn_counter) /(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
					end
					rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
				elseif rxn_mechanism == "UNIBI"
					parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter"]
					try
						ids = sample(non_gene_species_ids, 2, replace = false)
					catch
						temp_ids = collect((nSpecies_gene+1):nSpecies)
						ids = sample(temp_ids, 2, replace = false)
					end
					ids_in = []
					prds_id  = ids[1]
					prds_id2 = ids[2]
					rcts = [specs_input_selec]
					prds = ["S$prds_id", "S$prds_id2"]
					cats = []
					rxn_mech = ["kf_J$rxn_counter * $specs_input_selec - kr_J$rxn_counter * S$prds_id * S$prds_id2"]
					rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
				elseif rxn_mechanism == "BIUNI"
					parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter"]
					try
						ids = sample(non_gene_species_ids, 1, replace = false)
						ids_in = sample(non_gene_species_ids_in, 1, replace = false)
						append!(non_gene_species_ids_in, ids)
					catch
						temp_ids = collect((nSpecies_gene+1):nSpecies)
						ids = sample(temp_ids, 1, replace = false)
						ids_in = sample(species_ids_in, 1, replace = false)
					end
					if rand() < 0.5
						rcts_id2 = ids_in[1]
						prds_id  = ids[1]
						rcts = [specs_input_selec, "S$rcts_id2"]
						prds = ["S$prds_id"]
						cats = []
						rxn_mech = ["kf_J$rxn_counter * $specs_input_selec * S$rcts_id2 - kr_J$rxn_counter * S$prds_id"]
					else
						rcts_id  = ids_in[1]
						prds_id  = ids[1]
						rcts = ["S$rcts_id", specs_input_selec]
						prds = ["S$prds_id"]
						cats = []
						rxn_mech = ["kf_J$rxn_counter * S$rcts_id * $specs_input_selec - kr_J$rxn_counter * S$prds_id"]
					end
					rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
				elseif rxn_mechanism == "BIBI"
					parameters = ["kf_J$rxn_counter", "kr_J$rxn_counter"]
					try
						ids = sample(non_gene_species_ids, 2, replace = false)
						ids_in = sample(non_gene_species_ids_in, 1, replace = false)
						append!(non_gene_species_ids_in, ids)
					catch
						temp_ids = collect((nSpecies_gene+1):nSpecies)
						ids = sample(temp_ids, 2, replace = false)
						ids_in = sample(species_ids_in, 1, replace = false)
					end
					if rand() < 0.5
						rcts_id2 = ids_in[1]
						prds_id  = ids[1]
						prds_id2 = ids[2]
						rcts = [specs_input_selec, "S$rcts_id2"]
						prds = ["S$prds_id", "S$prds_id2"]
						cats = []
						rxn_mech = ["kf_J$rxn_counter * $specs_input_selec * S$rcts_id2 - kr_J$rxn_counter * S$prds_id * S$prds_id2"]
					else
						rcts_id  = ids_in[1]
						prds_id  = ids[1]
						prds_id2 = ids[2]
						rcts = ["S$rcts_id", specs_input_selec]
						prds = ["S$prds_id", "S$prds_id2"]
						cats = []
						rxn_mech = ["kf_J$rxn_counter * S$rcts_id * $specs_input_selec - kr_J$rxn_counter * S$prds_id * S$prds_id2"]
					end
					rxn_specs["r$rxn_counter"] = Dict{String, Array{String, 1}}("parameters" => parameters, "rcts" => rcts, "prds" => prds, "cats" => cats, "rxn_mech" => rxn_mech)
				elseif rxn_mechanism == "CIRCLE"
					parameters1 = ["kf_J$rxn_counter", "kr_J$rxn_counter", "K1_J$rxn_counter", "K2_J$rxn_counter"]
					parameters2 = ["kf_J$(rxn_counter+1)", "kr_J$(rxn_counter+1)", "K1_J$(rxn_counter+1)", "K2_J$(rxn_counter+1)"]
					try
						ids = sample(species_ids, 2, replace = false)
						ids_in = sample(species_ids_in, 1, replace = false)
						append!(species_ids_in, ids)
					catch
						temp_ids = collect(1:nSpecies)
						ids = sample(temp_ids, 2, replace = false)
						ids_in = sample(species_ids_in, 1, replace = false)
					end
					if rand() < 0.5
						rcts_id = ids[1]
						prds_id = ids[2]
						cats_id2 = ids_in[1]
						rcts = ["S$rcts_id"]
						prds = ["S$prds_id"]
						cats1 = [specs_input_selec]
						cats2 = ["S$cats_id2"]
						rxn_mech1 = ["$specs_input_selec * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter - kr_J$rxn_counter*S$prds_id/K2_J$rxn_counter) /(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
						rxn_mech2 = ["S$cats_id2 * (kf_J$(rxn_counter+1)*S$prds_id/K2_J$(rxn_counter+1) - kr_J$(rxn_counter+1)*S$rcts_id/K1_J$(rxn_counter+1)) /(1 + S$prds_id/K2_J$(rxn_counter+1) + S$rcts_id/K1_J$(rxn_counter+1))"]
					else
						rcts_id = ids[1]
						prds_id = ids[2]
						cats_id = ids_in[1]
						rcts = ["S$rcts_id"]
						prds = ["S$prds_id"]
						cats1 = ["S$cats_id"]
						cats2 = [specs_input_selec]
						rxn_mech1 = ["S$cats_id * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter - kr_J$rxn_counter*S$prds_id/K2_J$rxn_counter) /(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
						rxn_mech2 = ["$specs_input_selec * (kf_J$(rxn_counter+1)*S$prds_id/K2_J$(rxn_counter+1) - kr_J$(rxn_counter+1)*S$rcts_id/K1_J$(rxn_counter+1)) /(1 + S$prds_id/K2_J$(rxn_counter+1) + S$rcts_id/K1_J$(rxn_counter+1))"]
					end
					rxn_specs["r$rxn_counter"]     = Dict{String, Array{String, 1}}("parameters" => parameters1, "rcts" => rcts, "prds" => prds, "cats" => cats1, "rxn_mech" => rxn_mech1)
					rxn_specs["r$(rxn_counter+1)"] = Dict{String, Array{String, 1}}("parameters" => parameters2, "rcts" => prds, "prds" => rcts, "cats" => cats2, "rxn_mech" => rxn_mech2)
					rxn_counter += 1
				else rxn_mechanism == "DBCIRCLE"
                    # parameters1 = ["kf_J$rxn_counter", "K1_J$rxn_counter", "K2_J$rxn_counter"]
                    # parameters2 = ["kf_J$(rxn_counter+1)"]
                    # parameters3 = ["kf_J$(rxn_counter+2)", "K3_J$(rxn_counter+2)", "K4_J$(rxn_counter+2)", "K5_J$(rxn_counter+2)"]
                    # parameters4 = ["kf_J$(rxn_counter+3)"]
					# try
					# 	ids = sample(species_ids, 3, replace = false)
					# 	ids_in = sample(species_ids_in, 1, replace = false)
					# 	append!(species_ids_in, ids)
					# catch
					# 	temp_ids = collect(1:nSpecies)
					# 	ids = sample(temp_ids, 3, replace = false)
					# 	ids_in = sample(species_ids_in, 1, replace = false)
					# end
					# if rand() < 0.5
					# 	rcts_id  = ids[1]
					# 	prds_id  = ids[2]
					# 	prds_id2 = ids[3]
					# 	cats_id2 = ids_in[1]
					# 	rcts1 = ["S$rcts_id"]
					# 	prds1 = ["S$prds_id"]
					# 	prds2 = ["S$prds_id2"]
					# 	cats1 = [specs_input_selec]
					# 	cats2 = ["S$cats_id2"]
                    #     rxn_mech1 = ["$specs_input_selec  * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter)/(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
                    #     rxn_mech2 = ["$specs_input_selec  * (kf_J$(rxn_counter+1)*S$prds_id/K2_J$rxn_counter)/(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
                    #     rxn_mech3 = ["S$cats_id2*(kf_J$(rxn_counter+2)*S$prds_id2/K3_J$(rxn_counter+2))/(1 + S$prds_id2/K3_J$(rxn_counter+2) + S$prds_id/K4_J$(rxn_counter+2) + S$rcts_id/K5_J$(rxn_counter+2))"]
                    #     rxn_mech4 = ["S$cats_id2*(kf_J$(rxn_counter+3)*S$prds_id/K4_J$(rxn_counter+2))/(1 + S$prds_id2/K3_J$(rxn_counter+2) + S$prds_id/K4_J$(rxn_counter+2) + S$rcts_id/K5_J$(rxn_counter+2))"]
					# else
					# 	rcts_id  = ids[1]
					# 	prds_id  = ids[2]
					# 	prds_id2 = ids[3]
					# 	cats_id  = ids_in[1]
					# 	rcts1 = ["S$rcts_id"]
					# 	prds1 = ["S$prds_id"]
					# 	prds2 = ["S$prds_id2"]
					# 	cats1 = ["S$cats_id"]
					# 	cats2 = [specs_input_selec]
                    #     rxn_mech1 = ["S$cats_id * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter)/(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
                    #     rxn_mech2 = ["S$cats_id * (kf_J$(rxn_counter+1)*S$prds_id/K2_J$rxn_counter)/(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
                    #     rxn_mech3 = ["$specs_input_selec * (kf_J$(rxn_counter+2)*S$prds_id2/K3_J$(rxn_counter+2))/(1 + S$prds_id2/K3_J$(rxn_counter+2) + S$prds_id/K4_J$(rxn_counter+2) + S$rcts_id/K5_J$(rxn_counter+2))"]
                    #     rxn_mech4 = ["$specs_input_selec * (kf_J$(rxn_counter+3)*S$prds_id/K4_J$(rxn_counter+2))/(1 + S$prds_id2/K3_J$(rxn_counter+2) + S$prds_id/K4_J$(rxn_counter+2) + S$rcts_id/K5_J$(rxn_counter+2))"]
					# end
                    parameters1 = ["kf_J$rxn_counter", "kr_J$rxn_counter", "K1_J$rxn_counter", "K2_J$rxn_counter"]
                    parameters2 = ["kf_J$(rxn_counter+1)", "kr_J$(rxn_counter+1)", "K1_J$(rxn_counter+1)", "K2_J$(rxn_counter+1)"]
                    parameters3 = ["kf_J$(rxn_counter+2)", "kr_J$(rxn_counter+2)", "K1_J$(rxn_counter+2)", "K2_J$(rxn_counter+2)"]
                    parameters4 = ["kf_J$(rxn_counter+3)", "kr_J$(rxn_counter+3)", "K1_J$(rxn_counter+3)", "K2_J$(rxn_counter+3)"]
                    try
                        ids = sample(species_ids, 3, replace = false)
                        species_ids_in = rv_specs(species_ids_in, ids)
                        ids_in = sample(species_ids_in, 1, replace = false)
                        append!(species_ids_in, ids)
                    catch
                        temp_ids = collect(1:nSpecies)
                        temp_ids = rv_specs(temp_ids, specs_input_selec_ids)
                        ids = sample(temp_ids, 3, replace = false)
                        species_ids_in = rv_specs(temp_ids, ids)
                        ids_in = sample(species_ids_in, 1, replace = false)
                    end
                    if rand() < 0.5
                        rcts_id  = ids[1]
                        prds_id  = ids[2]
                        prds_id2 = ids[3]
                        cats_id2 = ids_in[1]
                        rv_ids = []
                        append!(rv_ids, rcts_id)
                        append!(rv_ids, prds_id)
                        append!(rv_ids, prds_id2)
                        append!(rv_ids, cats_id2)
                        species_ids = rv_specs(species_ids, rv_ids)
                        non_gene_species_ids = rv_specs(non_gene_species_ids, rv_ids)
                        rcts1 = ["S$rcts_id"]
                        prds1 = ["S$prds_id"]
                        prds2 = ["S$prds_id2"]
                        cats1 = [specs_input_selec]
                        cats2 = ["S$cats_id2"]
                        rxn_mech1 = ["$specs_input_selec  * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter - kr_J$rxn_counter*S$prds_id/K2_J$rxn_counter) /(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
                        rxn_mech2 = ["$specs_input_selec  * (kf_J$(rxn_counter+1)*S$prds_id/K1_J$(rxn_counter+1) - kr_J$(rxn_counter+1)*S$prds_id2/K2_J$(rxn_counter+1)) /(1 + S$prds_id/K1_J$(rxn_counter+1) + S$prds_id2/K2_J$(rxn_counter+1))"]
                        rxn_mech3 = ["S$cats_id2 * (kf_J$(rxn_counter+2)*S$prds_id2/K1_J$(rxn_counter+2) - kr_J$(rxn_counter+2)*S$prds_id/K2_J$(rxn_counter+2)) /(1 + S$prds_id2/K1_J$(rxn_counter+2) + S$prds_id/K2_J$(rxn_counter+2))"]
                        rxn_mech4 = ["S$cats_id2 * (kf_J$(rxn_counter+3)*S$prds_id/K1_J$(rxn_counter+3) - kr_J$(rxn_counter+3)*S$rcts_id/K2_J$(rxn_counter+3)) /(1 + S$prds_id/K1_J$(rxn_counter+3) + S$rcts_id/K2_J$(rxn_counter+3))"]
                    else
                        rcts_id  = ids[1]
                        prds_id  = ids[2]
                        prds_id2 = ids[3]
                        cats_id  = ids_in[1]
                        rv_ids = []
                        append!(rv_ids, rcts_id)
                        append!(rv_ids, prds_id)
                        append!(rv_ids, prds_id2)
                        append!(rv_ids, cats_id)
                        species_ids = rv_specs(species_ids, rv_ids)
                        non_gene_species_ids = rv_specs(non_gene_species_ids, rv_ids)
                        rcts1 = ["S$rcts_id"]
                        prds1 = ["S$prds_id"]
                        prds2 = ["S$prds_id2"]
                        cats1 = ["S$cats_id"]
                        cats2 = [specs_input_selec]
                        rxn_mech1 = ["S$cats_id  * (kf_J$rxn_counter*S$rcts_id/K1_J$rxn_counter - kr_J$rxn_counter*S$prds_id/K2_J$rxn_counter) /(1 + S$rcts_id/K1_J$rxn_counter + S$prds_id/K2_J$rxn_counter)"]
                        rxn_mech2 = ["S$cats_id  * (kf_J$(rxn_counter+1)*S$prds_id/K1_J$(rxn_counter+1) - kr_J$(rxn_counter+1)*S$prds_id2/K2_J$(rxn_counter+1)) /(1 + S$prds_id/K1_J$(rxn_counter+1) + S$prds_id2/K2_J$(rxn_counter+1))"]
                        rxn_mech3 = ["$specs_input_selec * (kf_J$(rxn_counter+2)*S$prds_id2/K1_J$(rxn_counter+2) - kr_J$(rxn_counter+2)*S$prds_id/K2_J$(rxn_counter+2)) /(1 + S$prds_id2/K1_J$(rxn_counter+2) + S$prds_id/K2_J$(rxn_counter+2))"]
                        rxn_mech4 = ["$specs_input_selec * (kf_J$(rxn_counter+3) * S$prds_id/K1_J$(rxn_counter+3) - kr_J$(rxn_counter+3) * S$rcts_id/K2_J$(rxn_counter+3)) /(1 + S$prds_id/K1_J$(rxn_counter+3) + S$rcts_id/K2_J$(rxn_counter+3))"]
                    end
					rxn_specs["r$rxn_counter"]     = Dict{String, Array{String, 1}}("parameters" => parameters1, "rcts" => rcts1, "prds" => prds1, "cats" => cats1, "rxn_mech" => rxn_mech1)
					rxn_specs["r$(rxn_counter+1)"] = Dict{String, Array{String, 1}}("parameters" => parameters2, "rcts" => prds1, "prds" => prds2, "cats" => cats1, "rxn_mech" => rxn_mech2)
					rxn_specs["r$(rxn_counter+2)"] = Dict{String, Array{String, 1}}("parameters" => parameters3, "rcts" => prds2, "prds" => prds1, "cats" => cats2, "rxn_mech" => rxn_mech3)
					rxn_specs["r$(rxn_counter+3)"] = Dict{String, Array{String, 1}}("parameters" => parameters4, "rcts" => prds1, "prds" => rcts1, "cats" => cats2, "rxn_mech" => rxn_mech4)
					rxn_counter += 3
				end
			end

			append!(ids_tot, ids)
			append!(ids_tot, ids_in)

			rxn_counter_end = rxn_counter
			specs_input = []
			for i = (rxn_counter_start+1):rxn_counter_end
				prds = rxn_specs["r$i"]["prds"]
				append!(specs_input, prds)
				if size(specs_input)[1] == 4 # if select "DBCIRCLE", remove the first species
					deleteat!(specs_input, 1)
				end
			end

			if rxn_counter >= nRxns
				valid_nRxns = true
			end

        end

        if rxn_counter == nRxns
            # last prds to connect to next output layer with its rcts or cats
            specs_r4 = []
            rcts = rxn_specs["r4"]["rcts"]
            cats = rxn_specs["r4"]["cats"]
            append!(specs_r4, rcts)
            append!(specs_r4, cats)


            specs_last = []
            rcts = rxn_specs["r$rxn_counter"]["rcts"]
            prds = rxn_specs["r$rxn_counter"]["prds"]
            append!(specs_last, rcts)
            append!(specs_last, cats)

            rand_specs_last = rand(1: size(specs_last)[1])
            rand_specs_r4 = rand(1: size(specs_r4)[1])
            specs_r4[rand_specs_r4] = specs_last[rand_specs_last]


            # count the number of effective steps
            #S_in is in reaction r1, S_out is in reaction r2 and r3
            #step 1: does r1 and r2+r3 have common species?
            specs_in = []
            rcts = rxn_specs["r1"]["rcts"]
            prds = rxn_specs["r1"]["prds"]
            cats = rxn_specs["r1"]["cats"]
            specs_in = append!(specs_in, rcts)
            specs_in = append!(specs_in, prds)
            specs_in = append!(specs_in, cats)

            specs_out = []
            rcts = rxn_specs["r2"]["rcts"]
            prds = rxn_specs["r2"]["prds"]
            cats = rxn_specs["r2"]["cats"]
            cats2= rxn_specs["r3"]["cats"]
            specs_out = append!(specs_out, rcts)
            specs_out = append!(specs_out, prds)
            specs_out = append!(specs_out, cats)
            specs_out = append!(specs_out, cats2)

            common_specs = intersect(specs_in, specs_out)
            global num_common_specs = size(common_specs)[1]

            eff_rxs = Array{Any,1}
            eff_rxs_num = 3 # 3 = r1(S_in)+r2(S_out)+r3(S_out)
            eff_rxs = ["r1", "r2", "r3"]

            rxs_in  = []
            rxs_out = []
            global trial = 1

            while num_common_specs == 0 && trial <= nRxns
                trial += 1
                #global specs_in, specs_out
                rxs_in  = []
                rxs_out = []
                for i = 4:rxn_counter # r1 includes S_in, r2 and r3 include S_out
                    specs = []
                    rcts = rxn_specs["r$i"]["rcts"]
                    prds = rxn_specs["r$i"]["prds"]
                    cats = rxn_specs["r$i"]["cats"]
                    specs = append!(specs, rcts)
                    specs = append!(specs, prds)
                    specs = append!(specs, cats)
                    if size(intersect(specs_out, specs))[1] > 0
                        append!(rxs_out, ["r$i"])
                    end
                end

                specs_out = []
                for j = 1:size(rxs_out)[1]
                    rcts = rxn_specs[rxs_out[j]]["rcts"]
                    prds = rxn_specs[rxs_out[j]]["prds"]
                    cats = rxn_specs[rxs_out[j]]["cats"]
                    specs_out = append!(specs_out, rcts)
                    specs_out = append!(specs_out, prds)
                    specs_out = append!(specs_out, cats)
                end

                common_specs = intersect(specs_in, specs_out)
                global num_common_specs = size(common_specs)[1]

                #println("num_common_specs:", num_common_specs)
            end
        else
            global num_common_specs = 0
        end
    end

    global eff_steps = trial

    return (rxn_specs, eff_steps)
end

function negativeConcentration(rr)
    C_S = RoadRunner.getFloatingSpeciesConcentrations(rr)
    nSpecies_boundary = RoadRunner.getNumberOfBoundarySpecies(rr)
    nSpecies_floating = RoadRunner.getNumberOfFloatingSpecies(rr)
    neg_c = 0
    #for i = 0:(nSpecies_floating-1)
	for i = 1:nSpecies_floating
        #C_element = getVectorElement(C_S, i)
		C_element = C_S[i]
        if C_element < 0
            neg_c = 1
        end
    end
    #freeVector(C_S)
    C_B = RoadRunner.getBoundarySpeciesConcentrations(rr)
    #for i = 0:(nSpecies_boundary-1)
	for i = 1:nSpecies_boundary
        #C_element = getVectorElement(C_B, i)
		C_element = C_B[i]
        if C_element < 0
            neg_c = 1
        end
    end
    #freeVector(C_B)
    return neg_c
end


# Main code starts from here
f = open("timer.txt", "w+");
tim=@elapsed begin # @elapsed returns time in second
    goodSample = 0
    while goodSample < sampleSize
        global species = ["S$i" for i = 1:nSpecies]
        global gene_species = species[1:nSpecies_gene]
        global rr_real = RoadRunner.createRRInstance()

        RoadRunner.addCompartment(rr_real, "compartment", 1.0, false)

        for s in species
            RoadRunner.addSpecies(rr_real, s, "compartment", 0.1, "", false)
            RoadRunner.setBoundary(rr_real, s, false, false)
        end
        RoadRunner.addSpecies(rr_real, "S_in", "compartment", 0.1, "", false)

        RoadRunner.setBoundary(rr_real, "S_in", true, false)
        #default is boundary species

        RoadRunner.addSpecies(rr_real, "S_out", "compartment", 0.1, "", false)

        RoadRunner.setBoundary(rr_real, "S_out", false, false)

        try
            (rxn_specs, eff_steps) = randomNetwork(nSpecies, nSpecies_gene, nRxns)

            st = zeros((nSpecies, nRxns))
            for i = 1: rxn_counter
                rct_list = rxn_specs["r$i"]["rcts"]
                prd_list = rxn_specs["r$i"]["prds"]
                rct_num = size(rct_list)[1]
                prd_num = size(prd_list)[1]
                for j = 1: rct_num
                    for k = 1 : nSpecies
                        if rct_list[j] == "S$k"
                            st[k,i] = -1
                        end
                    end
                end
                for j = 1: prd_num
                    for k = 1: nSpecies
                        if prd_list[j] == "S$k"
                            st[k,i] = 1
                        end
                    end
                end
            end

            for i = 1: nSpecies
                x = st[[i],:]
                if all(==(0), x)

                else
                    if all(>=(0), x)
                        RoadRunner.setBoundary(rr_real, "S$i", true, false)
                    elseif all(<=(0), x)
                        RoadRunner.setBoundary(rr_real, "S$i", true, false)
                    end
                end
            end

            for i = 1: nRxns
                parameters = rxn_specs["r$i"]["parameters"]
                parameters_len = size(parameters)[1]
                for j = 1:parameters_len
                    RoadRunner.addParameter(rr_real, parameters[j], 0.1, false)
                end
            end

            for i = 1: nRxns
                rcts = rxn_specs["r$i"]["rcts"]
                prds = rxn_specs["r$i"]["prds"]
                rxn_mech = rxn_specs["r$i"]["rxn_mech"]
                i == nRxns ? regen = true : regen = false
                RoadRunner.addReaction(rr_real, "r$i", rcts, prds, rxn_mech[1], regen)
            end

            P_num = RoadRunner.getNumberOfGlobalParameters(rr_real)
            nSpecies_boundary = RoadRunner.getNumberOfBoundarySpecies(rr_real)
            nSpecies_floating = RoadRunner.getNumberOfFloatingSpecies(rr_real)
            #assign random values to spieces and parameters
            for i = 0:(nSpecies_floating-1)
                RoadRunner.setFloatingSpeciesInitialConcentrationByIndex(rr_real, i, rnd_species*rand()) # random number [0,1)
            end
            for i = 0:(nSpecies_boundary-1)
                RoadRunner.setBoundarySpeciesByIndex(rr_real, i, rnd_species*rand())
            end
            for i = 0:(P_num-1)
                RoadRunner.setGlobalParameterByIndex(rr_real, i, rnd_parameter*rand())
            end

			neg_c = negativeConcentration(rr_real)
			if neg_c == 0
				##save sbml files before steadyState
			    sbml_sample0 = RoadRunner.getCurrentSBML(rr_real)
				f_sample0 = open("sampleNetwork0-$(goodSample+1).xml", "w+");
				write(f_sample0, sbml_sample0)
				close(f_sample0)

				RoadRunner.setConfigInt("LOADSBMLOPTIONS_CONSERVED_MOIETIES", 1)
				# Check if the ground truth model has a steady state
				try
					RoadRunner.steadyState(rr_real)
					#check if it is a valid steadystated
                    C_B = RoadRunner.getBoundarySpeciesConcentrations(rr_real)
                    #C_Sin = getVectorElement(C_B, nSpecies_boundary-1)
                    C_Sin = C_B[nSpecies_boundary]
                    C_S = RoadRunner.getFloatingSpeciesConcentrations(rr_real)
                    #C_Sout = getVectorElement(C_S, nSpecies_floating-1)
                    C_Sout = C_S[nSpecies_floating]
                    #freeVector(C_B)
                    #freeVector(C_S)
                    RoadRunner.setBoundarySpeciesByIndex(rr_real, nSpecies_boundary-1, C_Sin*concentration_perturb)

                    RoadRunner.setConfigInt("LOADSBMLOPTIONS_CONSERVED_MOIETIES", 1)
                    try
                        RoadRunner.steadyState(rr_real)
                        C_S = RoadRunner.getFloatingSpeciesConcentrations(rr_real)
                        #C_Sout_after = getVectorElement(C_S, nSpecies_floating-1)
                        C_Sout_after = C_S[nSpecies_floating]
                        C_Sout_change = C_Sout_after - C_Sout
                        #freeVector(C_S)
                        if (abs.(C_Sout_change) > 0 && eff_steps >= 3)
                            global goodSample += 1
                            sbml_sample = RoadRunner.getCurrentSBML(rr_real)
                            f_sample = open("sampleNetwork-$goodSample.xml", "w+");
                            write(f_sample, sbml_sample)
                            close(f_sample)
                        end
                    catch e
                        continue
                    end

				catch e
					continue
				end
		    else
                println("Error message: negative concentrations")
            end
        catch e
            println(e)
        finally
            RoadRunner.freeRRInstance(rr_real)
        end
    end
end
write(f, "$tim\n");
close(f);