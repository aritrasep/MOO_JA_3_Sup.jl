###############################################################################
#                                                                             #
#  This file is part of the julia module for Multi Objective Optimization     #
#  (c) Copyright 2017 by Aritra Pal                                           #
#                                                                             #
#  Permission is granted for academic research use.  For other uses,          #
#  contact the authors for licensing options.                                 #
#                                                                             #
#  Use at your own risk. I make no guarantees about the correctness or        #          
#  usefulness of this code.                                                   #
#                                                                             #
###############################################################################

import Modof, Modolib, Modoalgos, Match, DataFrames, DataFramesMeta
@everywhere using Modof, Modolib, Modoalgos, Match, DataFrames, DataFramesMeta

params = Dict()
params[:algorithm] = "BBM"
params[:type_of_operation] = "lex_min"
params[:mip_emphasis] = true
params[:type_of_boost] = "NB"
params[:turn_off_cplex_heur] = true
params[:timelimit] = Inf
params[:total_threads] = 1

@everywhere function solve_bopinstance_exactly(params)
	if params[:read_instance] != :read_boscp_xavier
		instance, true_non_dom_sols = eval(params[:read_instance])(params[:read_instance_input])
		println("Solving Instance = $(params[:read_instance_input])")
	else
		instance, true_non_dom_sols = eval(params[:read_instance])(params[:read_instance_input]...)
		println("Solving Instance = $(params[:read_instance_input]...)")
	end
	
	if params[:total_threads] == 1
		println("Solving Instance on a Single Thread")
		master0::Float64 = time()
		non_dom_sols, rects, ip_solved = compute_exact_frontier(instance, params)
		master1::Float64 = time()
		println("Finished")
		non_dom_sols = select_non_dom_sols(non_dom_sols)
		slave0::Float64, slave1::Float64 = 0.0, 0.0
	else
		println("Decomposing Instance into $(params[:total_threads]) Subproblems")
		master0 = time()
		non_dom_sols_, rects_, ip_solved_ = Modoalgos.decompose_a_problem_into_smaller_problem(instance, params)
		println("Finished Decomposition")
		master1 = time()
		
		algorithm = " "
		if params[:algorithm] == "AUTO"
			algorithm = "AUTO"
			density::Float64 = length(findn(vec(instance.A)))/(size(instance.A)[1]*size(instance.A)[2])
			if density >= 0.42
				params[:algorithm] = "EBBM"
				params[:type_of_boost] = "O_SH"
				params[:turn_off_cplex_heur] = true
				params[:type_of_operation] = "lex_min"
				params[:mip_emphasis] = true
			else
				params[:algorithm] = "UDECM"
				params[:type_of_boost] = "NB"
				params[:turn_off_cplex_heur] = false
				params[:type_of_operation] = "dynamic"
				params[:mip_emphasis] = false
				if floor(abs(log10(size(instance.A)[1]/size(instance.A)[2]))) > 1.0
					params[:algorithm] = "EBBM"
				end
			end
		end
		
		println("Started Solving Subproblem")
		slave0 = time()
		non_dom_sols, rects, ip_solved = @match params[:algorithm] begin
			"UDECM" || "BDECM" => Modoalgos.ϵ_constraint_method(rects_[params[:thread_id]], instance, non_dom_sols_[params[:thread_id]], params)
			_ => Modoalgos.balanced_box_method(rects_[params[:thread_id]], instance, non_dom_sols_[params[:thread_id]], params)
		end
		slave1 = time()
		ip_solved += ip_solved_
		println("Finished Solving Subproblem")
		non_dom_sols = select_non_dom_sols(non_dom_sols)
		if algorithm == "AUTO"
			params[:algorithm] = "AUTO"
			params[:type_of_boost] = "NB"
			params[:turn_off_cplex_heur] = false
			params[:type_of_operation] = "lex_min"
			params[:mip_emphasis] = false
		end
	end
	
	filename::String = "../results/Experimental_Results_$(myid()).csv"
	cor_::Float64 = cor(instance.c1, instance.c2)
	sparsity::Float64 = 1-(length(findn(vec(instance.A)))/(size(instance.A)[1]*size(instance.A)[2]))
	
	if isfile(filename) == false
		results::DataFrame = DataFrame(Problem_Type=params[:problem_type], Instance_Type=params[:instance_type], Instance_ID=params[:instance_id], Objectives=2, Correlation=cor_, Variables=size(instance.A)[2], Constraints=size(instance.A)[1], Sparsity=sparsity, Algorithm=params[:algorithm], Type_of_Operation=params[:type_of_operation], MIP_Emphasis=params[:mip_emphasis], Type_of_Boost=params[:type_of_boost], Cplex_Heur=!params[:turn_off_cplex_heur], Parallelism=params[:parallelism], Total_Threads=params[:total_threads], Thread_ID=params[:thread_id], Master_Time=master1-master0, Slave_Time=slave1-slave0, IP_Solved=ip_solved, Num_Non_Dom_Pts=size(true_non_dom_sols)[1], Num_Non_Dom_Pts_Found=length(non_dom_sols))
	else
		results = readtable(filename)
		results[:Slave_Time] = float(results[:Slave_Time])
		push!(results, [params[:problem_type], params[:instance_type], params[:instance_id], 2, cor_, size(instance.A)[2], size(instance.A)[1], sparsity, params[:algorithm], params[:type_of_operation], params[:mip_emphasis], params[:type_of_boost], !params[:turn_off_cplex_heur], params[:parallelism], params[:total_threads], params[:thread_id], master1-master0, slave1-slave0, ip_solved, size(true_non_dom_sols)[1], length(non_dom_sols)])
	end
	writetable(filename, results)
end

@everywhere instance, true_non_dom_sols = read_bokp_xavier1("2KP50-1A")

###############################################################################
# Testing all Methods                                                         #
###############################################################################

@everywhere function count_true_non_dom_sols(true_::Array{Float64, 2}, non_dom_sols::Vector{BOPSolution})
	count = 0
	inds = Int64[]
	for i in 1:size(true_)[1]
		for j in 1:length(non_dom_sols)
			if true_[i, :] == [non_dom_sols[j].obj_val1, non_dom_sols[j].obj_val2]
				push!(inds, j)
				count += 1
				break
			end
		end
	end
	count, non_dom_sols[inds]
end

@everywhere function testing_all_methods(instance::BOBPInstance, true_non_dom_sols::Array{Float64, 2})
	params = Dict()
	params[:algorithm] = "UDECM"
	params[:type_of_operation] = "lex_min"
	params[:mip_emphasis] = false
	params[:type_of_boost] = "NB"
	params[:turn_off_cplex_heur] = false
	params[:timelimit] = Inf
	params[:total_threads] = 1
	
	println("Testing 1 Thread UDECM_LM_NB_HON")
	non_dom_sols, rects = compute_exact_frontier(instance, params)
	non_dom_sols = select_non_dom_sols(non_dom_sols)
	count, non_dom_sols = count_true_non_dom_sols(true_non_dom_sols, non_dom_sols)
	if length(non_dom_sols) == size(true_non_dom_sols)[1]
		println("Test Successful")
	else
		println("Test Not Successful - Difference = $(size(true_non_dom_sols)[1]-count)")
	end
	
	println("Testing 1 Thread BBM_LM_NB_HON")
	params[:algorithm] = "BBM"
	non_dom_sols, rects = compute_exact_frontier(instance, params)
	non_dom_sols = select_non_dom_sols(non_dom_sols)
	count, non_dom_sols = count_true_non_dom_sols(true_non_dom_sols, non_dom_sols)
	if length(non_dom_sols) == size(true_non_dom_sols)[1]
		println("Test Successful")
	else
		println("Test Not Successful - Difference = $(size(true_non_dom_sols)[1]-count)")
	end
	
	println("Testing 1 Thread EBBM_LM_NB_HON")
	params[:algorithm] = "EBBM"
	non_dom_sols, rects = compute_exact_frontier(instance, params)
	non_dom_sols = select_non_dom_sols(non_dom_sols)
	count, non_dom_sols = count_true_non_dom_sols(true_non_dom_sols, non_dom_sols)
	if length(non_dom_sols) == size(true_non_dom_sols)[1]
		println("Test Successful")
	else
		println("Test Not Successful - Difference = $(size(true_non_dom_sols)[1]-count)")
	end
	
	println("Testing 1 Thread UDECM_AO_NB_HON")
	params[:algorithm] = "UDECM"
	params[:type_of_operation] = "augmented_operation"
	non_dom_sols, rects = compute_exact_frontier(instance, params)
	non_dom_sols = select_non_dom_sols(non_dom_sols)
	count, non_dom_sols = count_true_non_dom_sols(true_non_dom_sols, non_dom_sols)
	if length(non_dom_sols) == size(true_non_dom_sols)[1]
		println("Test Successful")
	else
		println("Test Not Successful - Difference = $(size(true_non_dom_sols)[1]-count)")
	end
	
	println("Testing 1 Thread EBBM_AO_NB_HON")
	params[:algorithm] = "EBBM"
	non_dom_sols, rects = compute_exact_frontier(instance, params)
	non_dom_sols = select_non_dom_sols(non_dom_sols)
	count, non_dom_sols = count_true_non_dom_sols(true_non_dom_sols, non_dom_sols)
	if length(non_dom_sols) == size(true_non_dom_sols)[1]
		println("Test Successful")
	else
		println("Test Not Successful - Difference = $(size(true_non_dom_sols)[1]-count)")
	end
	
	println("Testing 1 Thread UDECM_D_NB_HON")
	params[:algorithm] = "UDECM"
	params[:type_of_operation] = "dynamic"
	non_dom_sols, rects = compute_exact_frontier(instance, params)
	non_dom_sols = select_non_dom_sols(non_dom_sols)
	count, non_dom_sols = count_true_non_dom_sols(true_non_dom_sols, non_dom_sols)
	if length(non_dom_sols) == size(true_non_dom_sols)[1]
		println("Test Successful")
	else
		println("Test Not Successful - Difference = $(size(true_non_dom_sols)[1]-count)")
	end
	
	println("Testing 1 Thread EBBM_D_NB_HON")
	params[:algorithm] = "EBBM"
	non_dom_sols, rects = compute_exact_frontier(instance, params)
	non_dom_sols = select_non_dom_sols(non_dom_sols)
	count, non_dom_sols = count_true_non_dom_sols(true_non_dom_sols, non_dom_sols)
	if length(non_dom_sols) == size(true_non_dom_sols)[1]
		println("Test Successful")
	else
		println("Test Not Successful - Difference = $(size(true_non_dom_sols)[1]-count)")
	end
	
	println("Testing 1 Thread EBBM_LM2_NB_HON")
	params[:algorithm] = "EBBM"
	params[:type_of_operation] = "lex_min"
	params[:mip_emphasis] = true
	params[:type_of_boost] = "NB"
	params[:turn_off_cplex_heur] = false
	non_dom_sols, rects = compute_exact_frontier(instance, params)
	non_dom_sols = select_non_dom_sols(non_dom_sols)
	count, non_dom_sols = count_true_non_dom_sols(true_non_dom_sols, non_dom_sols)
	if length(non_dom_sols) == size(true_non_dom_sols)[1]
		println("Test Successful")
	else
		println("Test Not Successful - Difference = $(size(true_non_dom_sols)[1]-count)")
	end
	
	println("Testing 1 Thread UDECM_D_B_HON")
	params[:algorithm] = "UDECM"
	params[:type_of_operation] = "dynamic"
	params[:mip_emphasis] = false
	params[:type_of_boost] = "B"
	params[:turn_off_cplex_heur] = true
	non_dom_sols, rects = compute_exact_frontier(instance, params)
	non_dom_sols = select_non_dom_sols(non_dom_sols)
	count, non_dom_sols = count_true_non_dom_sols(true_non_dom_sols, non_dom_sols)
	if length(non_dom_sols) == size(true_non_dom_sols)[1]
		println("Test Successful")
	else
		println("Test Not Successful - Difference = $(size(true_non_dom_sols)[1]-count)")
	end
	
	println("Testing 1 Thread UDECM_D_O_SH_HON")
	params[:type_of_boost] = "O_SH"
	non_dom_sols, rects = compute_exact_frontier(instance, params)
	non_dom_sols = select_non_dom_sols(non_dom_sols)
	count, non_dom_sols = count_true_non_dom_sols(true_non_dom_sols, non_dom_sols)
	if length(non_dom_sols) == size(true_non_dom_sols)[1]
		println("Test Successful")
	else
		println("Test Not Successful - Difference = $(size(true_non_dom_sols)[1]-count)")
	end
	
	println("Testing 1 Thread EBBM_D_B_HON")
	params[:algorithm] = "EBBM"
	params[:type_of_boost] = "B"
	params[:turn_off_cplex_heur] = true
	non_dom_sols, rects = compute_exact_frontier(instance, params)
	non_dom_sols = select_non_dom_sols(non_dom_sols)
	count, non_dom_sols = count_true_non_dom_sols(true_non_dom_sols, non_dom_sols)
	if length(non_dom_sols) == size(true_non_dom_sols)[1]
		println("Test Successful")
	else
		println("Test Not Successful - Difference = $(size(true_non_dom_sols)[1]-count)")
	end
	
	println("Testing 1 Thread EBBM_D_O_SH_HON")
	params[:type_of_boost] = "O_SH"
	non_dom_sols, rects = compute_exact_frontier(instance, params)
	non_dom_sols = select_non_dom_sols(non_dom_sols)
	count, non_dom_sols = count_true_non_dom_sols(true_non_dom_sols, non_dom_sols)
	if length(non_dom_sols) == size(true_non_dom_sols)[1]
		println("Test Successful")
	else
		println("Test Not Successful - Difference = $(size(true_non_dom_sols)[1]-count)")
	end
	
	println("Testing 1 Thread AUTO")
	params[:algorithm] = "AUTO"
	non_dom_sols, rects = compute_exact_frontier(instance, params)
	non_dom_sols = select_non_dom_sols(non_dom_sols)
	count, non_dom_sols = count_true_non_dom_sols(true_non_dom_sols, non_dom_sols)
	if length(non_dom_sols) == size(true_non_dom_sols)[1]
		println("Test Successful")
	else
		println("Test Not Successful - Difference = $(size(true_non_dom_sols)[1]-count)")
	end
end

###############################################################################
# Experiments for All Methods                                                 #
###############################################################################

###############################################################################
## Single Thread                                                             ##
###############################################################################

@everywhere function generate_queue_for_detailed_experiments_on_all_instances()
	queue_ = Dict{Any, Any}[]	
	
	solver_params = Dict()
	solver_params[:algorithm] = ["UDECM", "BBM", "EBBM", "UDECM", "UDECM", "EBBM", "EBBM", "EBBM", "UDECM", "UDECM", "EBBM", "EBBM", "EBBM", "EBBM", "EBBM", "EBBM", "AUTO", "AUTO", "AUTO"]
	solver_params[:type_of_operation] = ["lex_min", "lex_min", "lex_min", "augmented_operation", "dynamic", "augmented_operation", "dynamic", "lex_min", "dynamic", "dynamic", "dynamic", "dynamic", "lex_min", "lex_min", "lex_min", "lex_min", "lex_min", "lex_min", "lex_min"]
	solver_params[:mip_emphasis] = [false, false, false, false, false, false, false, true, false, false, false, false, true, true, true, true, false, false, false]
	solver_params[:type_of_boost] = ["NB", "NB", "NB", "NB", "NB", "NB", "NB", "NB", "B", "O_SH", "B", "O_SH", "B", "O_SH", "B", "O_SH", "NB", "NB", "NB"]
	solver_params[:turn_off_cplex_heur] = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, false, false, false]
	solver_params[:total_threads] = [ones(16)..., 2, 3, 4]
	
	for i in [[100, 10], [200, 40], [600, 60]], j in ["a", "b", "c", "d"], k in 1:length(solver_params[:algorithm])
		for t in 1:solver_params[:total_threads][k]
			params = Dict()
			params[:read_instance] = :read_boscp_xavier
			params[:read_instance_input] = [i..., j]
			params[:problem_type] = "BOSCP"
			params[:instance_type] = "Xavier"
			params[:instance_id] = "SCP$(i[2])x$(i[1])_$(j)"
			params[:algorithm] = solver_params[:algorithm][k]
			params[:type_of_operation] = solver_params[:type_of_operation][k]
			params[:mip_emphasis] = solver_params[:mip_emphasis][k]
			params[:type_of_boost] = solver_params[:type_of_boost][k]
			params[:turn_off_cplex_heur] = solver_params[:turn_off_cplex_heur][k]
			params[:timelimit] = Inf
			params[:total_threads] = solver_params[:total_threads][k]
			params[:thread_id] = t		
			push!(queue_, params)
		end
	end
	
	instances = readdir(Pkg.dir("Modolib/instances/BOKP/Xavier_Instances/non dominated points/"))
	ind = Int64[]
	for i in 1:length(instances)
		instances[i]= split(instances[i],".")[1]
		if contains(instances[i], "2KP")
			push!(ind, i)
		end
    end
    for i in 1:length(instances), k in 1:length(solver_params[:algorithm])
    	for t in 1:solver_params[:total_threads][k]
			params = Dict()
    		if i in ind
				params[:read_instance] = :read_bokp_xavier1
			else
				params[:read_instance] = :read_bokp_xavier2
			end
			params[:read_instance_input] = instances[i]
			params[:problem_type] = "BOKP"
			params[:instance_type] = "Xavier"
			params[:instance_id] = instances[i]
			params[:algorithm] = solver_params[:algorithm][k]
			params[:type_of_operation] = solver_params[:type_of_operation][k]
			params[:mip_emphasis] = solver_params[:mip_emphasis][k]
			params[:type_of_boost] = solver_params[:type_of_boost][k]
			params[:turn_off_cplex_heur] = solver_params[:turn_off_cplex_heur][k]
			params[:timelimit] = Inf
			params[:total_threads] = solver_params[:total_threads][k]
			params[:thread_id] = t		
			push!(queue_, params)
		end
	end
	
	instances = readdir(Pkg.dir("Modolib/instances/BOSPP/Xavier_Instances/non dominated points/"))
	for i in 1:length(instances)
		instances[i]= split(instances[i],".")[1]
    end
    sort!(instances)
	for i in 1:length(instances), k in 1:length(solver_params[:algorithm])
		for t in 1:solver_params[:total_threads][k]
			params = Dict()
			params[:read_instance] = :read_bospp_xavier
			params[:read_instance_input] = instances[i]
			params[:problem_type] = "BOSPP"
			params[:instance_type] = "Xavier"
			params[:instance_id] = instances[i]
			params[:algorithm] = solver_params[:algorithm][k]
			params[:type_of_operation] = solver_params[:type_of_operation][k]
			params[:mip_emphasis] = solver_params[:mip_emphasis][k]
			params[:type_of_boost] = solver_params[:type_of_boost][k]
			params[:turn_off_cplex_heur] = solver_params[:turn_off_cplex_heur][k]
			params[:timelimit] = Inf
			params[:total_threads] = solver_params[:total_threads][k]
			params[:thread_id] = t		
			push!(queue_, params)
		end
	end
	
	for j in 1:20, i in 1:2, k in 1:length(solver_params[:algorithm])
		for t in 1:solver_params[:total_threads][k]
			params = Dict()
			if i == 1
				params[:read_instance] = :read_bokp_hadi
				params[:problem_type] = "BOKP"
			else
				params[:read_instance] = :read_boap_hadi
				params[:problem_type] = "BOAP"
			end
			params[:read_instance_input] = j
			params[:instance_type] = "Hadi"
			params[:instance_id] = "$(j)dat.txt"
			params[:algorithm] = solver_params[:algorithm][k]
			params[:type_of_operation] = solver_params[:type_of_operation][k]
			params[:mip_emphasis] = solver_params[:mip_emphasis][k]
			params[:type_of_boost] = solver_params[:type_of_boost][k]
			params[:turn_off_cplex_heur] = solver_params[:turn_off_cplex_heur][k]
			params[:timelimit] = Inf
			params[:total_threads] = solver_params[:total_threads][k]
			params[:thread_id] = t		
			push!(queue_, params)
		end
	end
	
	try
		files = readdir("../results/")
		files = files[[contains(files[i], "Experimental_Results") for i in 1:length(files)]]
		data = readtable("../results/$(files[1])")
		for i in 2:length(files)
    		data = vcat(data, readtable("../results/$(files[i])"))
		end
		data = data[names(data)[1:21]]
		unique!(data)
		sort!(data)
		[rm("../results/$(files[i])") for i in 1:length(files)]
		writetable("../results/Experimental_Results.csv", data)
		ind = Int64[]
		for i in 1:length(queue_)
			tmp = data[data[:Problem_Type] .== queue_[i][:problem_type],:]
    		tmp = tmp[tmp[:Instance_Type] .== queue_[i][:instance_type],:]
    		tmp = tmp[tmp[:Instance_ID] .== queue_[i][:instance_id],:]
    		tmp = tmp[tmp[:Algorithm] .== queue_[i][:algorithm],:]
    		tmp = tmp[tmp[:Type_of_Operation] .== queue_[i][:type_of_operation],:]    		
    		tmp = tmp[tmp[:MIP_Emphasis] .== queue_[i][:mip_emphasis],:]
    		tmp = tmp[tmp[:Type_of_Boost] .== queue_[i][:type_of_boost],:]
    		tmp = tmp[tmp[:Cplex_Heur] .!= queue_[i][:turn_off_cplex_heur],:]
 			tmp = tmp[tmp[:Parallelism] .== queue_[i][:parallelism],:]
    		tmp = tmp[tmp[:Total_Threads] .== queue_[i][:total_threads],:]
    		tmp = tmp[tmp[:Thread_ID] .== queue_[i][:thread_id],:]
    		if nrow(tmp) == 0
        		push!(ind, i)
        		continue
    		end
		end
		return unique(queue_[ind])
	catch
		return unique(queue_)
	end
end

@everywhere function solve_bopinstance_exactly(params, inds::Vector{Int64})
	for i in inds
		try
			solve_bopinstance_exactly(params[i])
		catch
		end
	end
end

@everywhere function generate_parallel_queues(inds::Vector{Int64})
	procs_::Vector{Int64} = procs()
	inds_ = Dict()
	i::Int64 = 1
	count::Int64 = 1
    while i <= length(inds)
    	p = procs_[count]
    	if p != myid() || length(procs_) == 1
        	if p in keys(inds_)
        		push!(inds_[p], inds[i])
        	else
        		inds_[p] = [inds[i]]
        	end
        	i += 1
       	end
       	count += 1
       	if count > length(procs_)
       		count = 1
       	end
    end
	inds_
end

@everywhere function run_parallel_experiments_on_all_instances()
	procs_::Vector{Int64} = procs()
	@sync begin
        for p in procs_
            @async begin
            	remotecall_fetch(testing_all_methods, p, instance, true_non_dom_sols)
            end
        end
    end
    queue = generate_queue_for_detailed_experiments_on_all_instances()
    inds = generate_parallel_queues([1:length(queue)...])
	@sync begin
        for p in procs_
        	if p != myid() || length(procs_) == 1
            	@async begin
            		remotecall_fetch(solve_bopinstance_exactly, p, queue, inds[p])
            	end
            end
        end
    end
end

#####################################################################
# Execution                                                         #
#####################################################################

#run_parallel_experiments_on_all_instances()
