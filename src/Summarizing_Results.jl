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

using Modof, DataFrames, DataFramesMeta, Match, PyCall, PyPlot

@pyimport seaborn as sns
sns.set_style("white")
sns.set_context("paper")
sns.set(font_scale=2.0)

plot_location = "../plots/"

#####################################################################
# Collecting the Experimental Results in one place                  #
#####################################################################

files = readdir("../results/")
files = files[[contains(files[i], "Experimental_Results") for i in 1:length(files)]]
data = readtable("../results/$(files[1])")
for i in 2:length(files)
	data = vcat(data, readtable("../results/$(files[i])"))
end
data = data[names(data)[1:21]]
data = data[completecases(data),:]
unique!(data)
sort!(data)
[rm("../results/$(files[i])") for i in 1:length(files)]
writetable("../results/Experimental_Results.csv", data)

#####################################################################
# Computing the results for AUTO for a Single Threads               #
#####################################################################

instances = unique(data[[:Problem_Type, :Instance_Type, :Instance_ID, :Correlation, :Objectives, :Variables, :Constraints, :Sparsity]])
algorithms = unique(data[[:Algorithm, :Type_of_Operation, :MIP_Emphasis, :Type_of_Boost, :Cplex_Heur, :Total_Threads]])
instance_inds = Int64[]
algorithm_inds = Int64[]
for i in 1:nrow(instances)
    push!(instance_inds, i)
    if (1-instances[:Sparsity][i]) >= 0.42
        push!(algorithm_inds, 10)
    else
        if floor(abs(log10(instances[:Constraints][i]/instances[:Variables][i]))) > 1.0
            push!(algorithm_inds, 4)
        else
            push!(algorithm_inds, 14)
        end
    end
end
results = hcat(instances[instance_inds,:], algorithms[algorithm_inds, :])
data_ = zeros(nrow(results), 5)
count = 1
for i in 1:nrow(results)
	tmp = data[data[:Problem_Type] .== results[:Problem_Type][i],:]
    tmp = tmp[tmp[:Instance_Type] .== results[:Instance_Type][i],:]
    tmp = tmp[tmp[:Instance_ID] .== results[:Instance_ID][i],:]
    tmp = tmp[tmp[:Algorithm] .== results[:Algorithm][i],:]
    tmp = tmp[tmp[:Type_of_Operation] .== results[:Type_of_Operation][i],:]
    tmp = tmp[tmp[:MIP_Emphasis] .== results[:MIP_Emphasis][i],:]
    tmp = tmp[tmp[:Type_of_Boost] .== results[:Type_of_Boost][i],:]
    tmp = tmp[tmp[:Cplex_Heur] .== results[:Cplex_Heur][i],:]
    tmp = tmp[tmp[:Total_Threads] .== results[:Total_Threads][i],:]
	if nrow(tmp) == 1
		data_[count, 1] = tmp[:Master_Time][1]
        data_[count, 2] = tmp[:Slave_Time][1]
        data_[count, 3] = tmp[:IP_Solved][1]
        data_[count, 4] = tmp[:Num_Non_Dom_Pts][1]
        data_[count, 5] = tmp[:Num_Non_Dom_Pts_Found][1]
	end
	count += 1
end
results[:Algorithm] = "AUTO"
results[:Type_of_Operation] = "lex_min"
results[:MIP_Emphasis] = false
results[:Type_of_Boost] = "NB"
results[:Cplex_Heur] = false
results[:Parallelism] = false
results[:Thread_ID] = 1
results[:Master_Time] = data_[:, 1]
results[:Slave_Time] = data_[:, 2]
results[:IP_Solved] = data_[:, 3]
results[:Num_Non_Dom_Pts] = data_[:, 4]
results[:Num_Non_Dom_Pts_Found] = data_[:, 5]
results = results[names(results)[[1:13..., 15, 14, 16:21...]]]
writetable("../results/Experimental_Results_0.csv", results)

#####################################################################
# Recollecting the Experimental Results in one place                #
#####################################################################

files = readdir("../results/")
files = files[[contains(files[i], "Experimental_Results") for i in 1:length(files)]]
data = readtable("../results/$(files[1])")
for i in 2:length(files)
	data = vcat(data, readtable("../results/$(files[i])"))
end
data = data[names(data)[1:21]]
data = data[completecases(data),:]
unique!(data)
sort!(data)
[rm("../results/$(files[i])") for i in 1:length(files)]
writetable("../results/Experimental_Results.csv", data)

#####################################################################
# Collecting the Experimental Results in one place                  #
#####################################################################

instances = unique(data[[:Problem_Type, :Instance_Type, :Instance_ID, :Correlation, :Objectives, :Variables, :Constraints, :Sparsity]])
algorithms = unique(data[[:Algorithm, :Type_of_Operation, :MIP_Emphasis, :Type_of_Boost, :Cplex_Heur, :Total_Threads]])
instance_inds = Int64[]
algorithm_inds = Int64[]
data_ = zeros(nrow(instances)*nrow(algorithms), 8)
count = 1
for i in 1:nrow(instances)
    tmp = data[data[:Problem_Type] .== instances[:Problem_Type][i],:]
    tmp = tmp[tmp[:Instance_Type] .== instances[:Instance_Type][i],:]
    tmp = tmp[tmp[:Instance_ID] .== instances[:Instance_ID][i],:]
    for j in 1:nrow(algorithms)
    	k = algorithms[:Total_Threads][j]
        tmp2 = tmp[tmp[:Algorithm] .== algorithms[:Algorithm][j],:]
        tmp2 = tmp2[tmp2[:Type_of_Operation] .== algorithms[:Type_of_Operation][j],:]
        tmp2 = tmp2[tmp2[:MIP_Emphasis] .== algorithms[:MIP_Emphasis][j],:]
        tmp2 = tmp2[tmp2[:Type_of_Boost] .== algorithms[:Type_of_Boost][j],:]
        tmp2 = tmp2[tmp2[:Cplex_Heur] .== algorithms[:Cplex_Heur][j],:]
        tmp2 = tmp2[tmp2[:Total_Threads] .== algorithms[:Total_Threads][j],:]        
        if nrow(tmp2) == k
            push!(instance_inds, i)
            push!(algorithm_inds, j)
            data_[count, 1] = maximum(tmp2[:Master_Time]) + maximum(tmp2[:Slave_Time])
            data_[count, 2:(2+k-1)] = [tmp2[:Slave_Time]...]
            data_[count, 2] += maximum(tmp2[:Master_Time])
            data_[count, 6] = sum(tmp2[:IP_Solved])
            data_[count, 7] = mean(tmp2[:Num_Non_Dom_Pts])
            data_[count, 8] = sum(tmp2[:Num_Non_Dom_Pts_Found]) - (k-1)
            count += 1
        end
    end
end
data_ = data_[1:length(instance_inds),:]
results = DataFrame()
results[:Problem_Type] = instances[:Problem_Type][instance_inds]
results[:Instance_Type] = instances[:Instance_Type][instance_inds]
results[:Instance_ID] = instances[:Instance_ID][instance_inds]
results[:Correlation] = instances[:Correlation][instance_inds]
results[:Objectives] = instances[:Objectives][instance_inds]
results[:Variables] = instances[:Variables][instance_inds]
results[:Constraints] = instances[:Constraints][instance_inds]
results[:Sparsity] = instances[:Sparsity][instance_inds]
results[:Algorithm] = algorithms[:Algorithm][algorithm_inds]
results[:Type_of_Operation] = algorithms[:Type_of_Operation][algorithm_inds]
results[:MIP_Emphasis] = algorithms[:MIP_Emphasis][algorithm_inds]
results[:Type_of_Boost] = algorithms[:Type_of_Boost][algorithm_inds]
results[:Cplex_Heur] = algorithms[:Cplex_Heur][algorithm_inds]
results[:Total_Threads] = algorithms[:Total_Threads][algorithm_inds]
results[:Total_Time] = round(data_[:,1],4)
results[:Thread_1_Slave_Time] = round(data_[:,2],4)
results[:Thread_2_Slave_Time] = round(data_[:,3],4)
results[:Thread_3_Slave_Time] = round(data_[:,4],4)
results[:Thread_4_Slave_Time] = round(data_[:,5],4)
results[:IP_Solved] = data_[:,6]
results[:Percent_of_Non_Dom_Pts_Found] = round((data_[:,8] ./ data_[:,7] ) * 100, 4)
sort!(results)
writetable("../results/UER1.csv", results)

#####################################################################
# Compressing Algorithms and Number of Threads                      #
# into a Single Column                                              #
#####################################################################

algorithm_ = String[]
for i in 1:nrow(results)
	push!(algorithm_, results[:Algorithm][i])
	if algorithm_[end] != "AUTO"
		algorithm_[end] = @match results[:Type_of_Operation][i] begin
			"lex_min" => "$(algorithm_[end])_LM"
			"augmented_operation" => "$(algorithm_[end])_AO"
			"dynamic" => "$(algorithm_[end])_D"
		end
		algorithm_[end] = @match results[:MIP_Emphasis][i] begin
			true => "$(algorithm_[end])2"
			_ => algorithm_[end]
		end
		algorithm_[end] = @match results[:Type_of_Boost][i] begin
			"NB" => "$(algorithm_[end])_NB"
			"B" => "$(algorithm_[end])_B"
			"O_SH" => "$(algorithm_[end])_OSH"
			"M_SH" => "$(algorithm_[end])_MSH"
		end
		algorithm_[end] = @match results[:Cplex_Heur][i] begin
			true => "$(algorithm_[end])_HON"
			_ => "$(algorithm_[end])_HOFF"
		end
	end
	algorithm_[end] = "$(algorithm_[end])_$(results[:Total_Threads][i])"
end
results[:Algorithm] = algorithm_
delete!(results, [:Type_of_Operation, :MIP_Emphasis, :Type_of_Boost, :Cplex_Heur, :Total_Threads])
sort!(results)
writetable("../results/UER2.csv", results)

#####################################################################
# Computing the Rank of Algorithms for a Single Thread              #
#####################################################################

minimum_time = Float64[]
ratio_to_minimum_time = Float64[]
rank = Int64[]
for i in 1:nrow(results)
	tmp = results[results[:Problem_Type] .== results[:Problem_Type][i],:]
    tmp = tmp[tmp[:Instance_Type] .== results[:Instance_Type][i],:]
    tmp = tmp[tmp[:Instance_ID] .== results[:Instance_ID][i],:]
	thread = parse(Int64, split(results[:Algorithm][i], "_")[end])
	tmp = tmp[ [ contains(tmp[:Algorithm][j], "_$(thread)") for j in 1:nrow(tmp) ],:]
	push!(minimum_time, minimum(tmp[:Total_Time]))
	push!(ratio_to_minimum_time, results[:Total_Time][i]/minimum_time[end])
	tmp2 = sort(tmp[:Total_Time])
	push!(rank, findin(tmp2, results[:Total_Time][i])[1])
end
results[:Minimum_Time] = minimum_time
results[:Ratio_to_Minimum_Time] = ratio_to_minimum_time
results[:Rank_of_the_Algorithm] = rank
results = results[names(results)[[1:9..., 18, 19, 10:17...]]]
sort!(results)
writetable("../results/UER3.csv", results)

#####################################################################
# Decomposing Results for Single Threads and Multiple Threads       #
#####################################################################

single_thread_results = results[[ contains(results[:Algorithm][i], "_1") for i in 1:nrow(results) ], :]
for i in 1:nrow(single_thread_results)
	single_thread_results[:Algorithm][i] = @match single_thread_results[:Algorithm][i] begin
		"AUTO_1" => "AUTO"
		"UDECM_LM_NB_HON_1" => "V1"
 		"BBM_LM_NB_HON_1" => "V2"
 		"EBBM_LM_NB_HON_1" => "V3" 
 		"UDECM_AO_NB_HON_1" => "V4"
 		"UDECM_D_NB_HON_1" => "V5"
 		"EBBM_LM2_NB_HON_1" => "V6"
 		"EBBM_AO_NB_HON_1" => "V7"
 		"EBBM_D_NB_HON_1" => "V8"
 		"UDECM_D_B_HON_1" => "V9"    
 		"UDECM_D_OSH_HON_1" => "V10"  
 		"EBBM_LM2_B_HON_1" => "V11"   
 		"EBBM_LM2_OSH_HON_1" => "V12" 
 		"EBBM_LM2_B_HOFF_1" => "V13"  
 		"EBBM_LM2_OSH_HOFF_1" => "V14"
 		"EBBM_D_B_HON_1" => "V15"     
 		"EBBM_D_OSH_HON_1" => "V16"  
		_ => single_thread_results[:Algorithm][i]
	end
end
sort!(single_thread_results)
writetable("../results/UER4.csv", single_thread_results)

multiple_thread_results = results[[ contains(results[:Algorithm][i], "AUTO") for i in 1:nrow(results) ], :]
multiple_thread_results[:Algorithm] = [ split(multiple_thread_results[:Algorithm][i], "_")[end] for i in 1:nrow(multiple_thread_results)]
sort!(multiple_thread_results)
writetable("../results/UER5.csv", multiple_thread_results)

#####################################################################
# Decomposing Results for Small and Large Scale Instances           #
#####################################################################

large_instance_minimum_time = 600.0

single_thread_results_of_small_instances = @where(single_thread_results, :Minimum_Time .< large_instance_minimum_time)
sort!(single_thread_results_of_small_instances)
writetable("../results/UER4_SI.csv", single_thread_results_of_small_instances)
multiple_thread_results_of_small_instances = @where(multiple_thread_results, :Minimum_Time .< large_instance_minimum_time)
sort!(multiple_thread_results_of_small_instances)
writetable("../results/UER5_SI.csv", multiple_thread_results_of_small_instances)

single_thread_results_of_large_instances = @where(single_thread_results, :Minimum_Time .>= large_instance_minimum_time)
sort!(single_thread_results_of_large_instances)
writetable("../results/UER4_LI.csv", single_thread_results_of_large_instances)
multiple_thread_results_of_large_instances = @where(multiple_thread_results, :Minimum_Time .>= large_instance_minimum_time)
sort!(multiple_thread_results_of_large_instances)
writetable("../results/UER5_LI.csv", multiple_thread_results_of_large_instances)

#####################################################################
# Generating the DataFrame for additional insights                  #
#####################################################################

data = zeros(17, 15)
algorithm = String[]
for i in 1:17
	if i <= 16
		push!(algorithm, "V$i")
		tmp = @where(single_thread_results, :Algorithm .== "V$i")
	else
		push!(algorithm, "AUTO")
		tmp = @where(single_thread_results, :Algorithm .== "AUTO")
	end
	data[i, 1] = mean(tmp[:Total_Time])
	data[i, 2] = mean(tmp[:Ratio_to_Minimum_Time])
	data[i, 3] = mean(tmp[:Rank_of_the_Algorithm])
	data[i, 4] = mean(tmp[:IP_Solved])
	data[i, 5] = mean(tmp[:Percent_of_Non_Dom_Pts_Found])
	if i <= 16
		tmp = @where(single_thread_results_of_small_instances, :Algorithm .== "V$i")
	else
		tmp = @where(single_thread_results_of_small_instances, :Algorithm .== "AUTO")
	end
	data[i, 6] = mean(tmp[:Total_Time])
	data[i, 7] = mean(tmp[:Ratio_to_Minimum_Time])
	data[i, 8] = mean(tmp[:Rank_of_the_Algorithm])
	data[i, 9] = mean(tmp[:IP_Solved])
	data[i, 10] = mean(tmp[:Percent_of_Non_Dom_Pts_Found])
	if i <= 16
		tmp = @where(single_thread_results_of_large_instances, :Algorithm .== "V$i")
	else
		tmp = @where(single_thread_results_of_large_instances, :Algorithm .== "AUTO")
	end
	data[i, 11] = mean(tmp[:Total_Time])
	data[i, 12] = mean(tmp[:Ratio_to_Minimum_Time])
	data[i, 13] = mean(tmp[:Rank_of_the_Algorithm])
	data[i, 14] = mean(tmp[:IP_Solved])
	data[i, 15] = mean(tmp[:Percent_of_Non_Dom_Pts_Found])
end
results = DataFrame()
results[:Algorithm] = algorithm
results[:All_Instances_Total_Time] = data[:, 1]
results[:All_Instances_Ratio] = data[:, 2]
results[:All_Instances_Rank] = data[:, 3]
results[:All_Instances_IP_Solved] = data[:, 4]
results[:All_Instances_Percent] = data[:, 5]
results[:Small_Instances_Total_Time] = data[:, 6]
results[:Small_Instances_Ratio] = data[:, 7]
results[:Small_Instances_Rank] = data[:, 8]
results[:Small_Instances_IP_Solved] = data[:, 9]
results[:Small_Instances_Percent] = data[:, 10]
results[:Large_Instances_Total_Time] = data[:, 11]
results[:Large_Instances_Ratio] = data[:, 12]
results[:Large_Instances_Rank] = data[:, 13]
results[:Large_Instances_IP_Solved] = data[:, 14]
results[:Large_Instances_Percent] = data[:, 15]
writetable("../results/UER6.csv", results)

#####################################################################
# Generating the Plots                                              #
#####################################################################

#####################################################################
## Plots for the 3 classes of Instances                            ##
#####################################################################

#####################################################################
### Plots for All Instances                                       ###
#####################################################################

plt_performance_profiles(single_thread_results, [:Problem_Type, :Instance_Type, :Instance_ID], ["AUTO", "V1", "V2", "V3"], "../plots/1_Thread_All_Instances_Performance_Profiles.eps")

single_thread_results = single_thread_results[[single_thread_results[:Algorithm][i] in ["AUTO", "V4", "V5", "V6", "V7", "V8", "V10", "V15", "V16"] for i in 1:nrow(single_thread_results)],:]

f = plt[:figure](figsize=(20, 10))
sns.boxplot(x=single_thread_results[:Algorithm], y=single_thread_results[:Total_Time])
ax = gca()
ax[:set_yscale]("log")
xlabel("Algorithms")
ylabel("Computing Time (secs)")
title("All Instances")
f[:savefig]("../plots/1_Thread_All_Instances_Computing_Time.eps")

f = plt[:figure](figsize=(20, 10))
sns.boxplot(x=single_thread_results[:Algorithm], y=single_thread_results[:Ratio_to_Minimum_Time])
ax = gca()
ax[:set_yscale]("log")
xlabel("Algorithms")
ylabel("Ratio to Minimum Time")
title("All Instances")
f[:savefig]("../plots/1_Thread_All_Instances_Ratio_to_Minimum_Time.eps")

f = plt[:figure](figsize=(20, 10))
sns.boxplot(x=single_thread_results[:Algorithm], y=single_thread_results[:Rank_of_the_Algorithm])
xlabel("Algorithms")
ylabel("Rank of the Algorithm")
title("All Instances")
f[:savefig]("../plots/1_Thread_All_Instances_Rank_of_the_Algorithm.eps")

f = plt[:figure](figsize=(20, 10))
sns.boxplot(x=single_thread_results[:Algorithm], y=single_thread_results[:IP_Solved])
xlabel("Algorithms")
ylabel("IPs Solved")
title("All Instances")
f[:savefig]("../plots/1_Thread_All_Instances_IP_Solved.eps")

f = plt[:figure](figsize=(20, 10))
sns.boxplot(x=single_thread_results[:Algorithm], y=single_thread_results[:Percent_of_Non_Dom_Pts_Found])
xlabel("Algorithms")
ylabel("Percentage of Nondominated Points Found")
title("All Instances")
f[:savefig]("../plots/1_Thread_All_Instances_Percent_of_Non_Dom_Pts_Found.eps")

#####################################################################
### Plots for Small Instances                                     ###
#####################################################################

plt_performance_profiles(single_thread_results_of_small_instances, [:Problem_Type, :Instance_Type, :Instance_ID], ["AUTO", "V1", "V2", "V3"], "../plots/1_Thread_Small_Instances_Performance_Profiles.eps")

single_thread_results_of_small_instances = single_thread_results_of_small_instances[[single_thread_results_of_small_instances[:Algorithm][i] in ["AUTO", "V3", "V4", "V5", "V6", "V7", "V8", "V12", "V15", "V16"] for i in 1:nrow(single_thread_results_of_small_instances)],:]

f = plt[:figure](figsize=(20, 10))
ax = sns.boxplot(x=single_thread_results_of_small_instances[:Algorithm], y=single_thread_results_of_small_instances[:Total_Time])
ax[:set_yscale]("log")
xlabel("Algorithms")
ylabel("Computing Time (secs)")
title("Small Instances")
f[:savefig]("../plots/1_Thread_Small_Instances_Computing_Time.eps")

f = plt[:figure](figsize=(20, 10))
ax = sns.boxplot(x=single_thread_results_of_small_instances[:Algorithm], y=single_thread_results_of_small_instances[:Ratio_to_Minimum_Time])
ax[:set_yscale]("log")
xlabel("Algorithms")
ylabel("Ratio to Minimum Time")
title("Small Instances")
f[:savefig]("../plots/1_Thread_Small_Instances_Ratio_to_Minimum_Time.eps")

f = plt[:figure](figsize=(20, 10))
sns.boxplot(x=single_thread_results_of_small_instances[:Algorithm], y=single_thread_results_of_small_instances[:Rank_of_the_Algorithm])
xlabel("Algorithms")
ylabel("Rank of the Algorithm")
title("Small Instances")
f[:savefig]("../plots/1_Thread_Small_Instances_Rank_of_the_Algorithm.eps")

f = plt[:figure](figsize=(20, 10))
sns.boxplot(x=single_thread_results_of_small_instances[:Algorithm], y=single_thread_results_of_small_instances[:IP_Solved])
xlabel("Algorithms")
ylabel("IPs Solved")
title("Small Instances")
f[:savefig]("../plots/1_Thread_Small_Instances_IP_Solved.eps")

f = plt[:figure](figsize=(20, 10))
sns.boxplot(x=single_thread_results_of_small_instances[:Algorithm], y=single_thread_results_of_small_instances[:Percent_of_Non_Dom_Pts_Found])
xlabel("Algorithms")
ylabel("Percentage of Nondominated Points Found")
title("Small Instances")
f[:savefig]("../plots/1_Thread_Small_Instances_Percent_of_Non_Dom_Pts_Found.eps")

#####################################################################
### Plots for Large Instances                                     ###
#####################################################################

plt_performance_profiles(single_thread_results_of_large_instances, [:Problem_Type, :Instance_Type, :Instance_ID], ["AUTO", "V1", "V2", "V3"], "../plots/1_Thread_Large_Instances_Performance_Profiles.eps")

single_thread_results_of_large_instances = single_thread_results_of_large_instances[[single_thread_results_of_large_instances[:Algorithm][i] in ["AUTO", "V5", "V6", "V7", "V8", "V9", "V10", "V13", "V15", "V16"] for i in 1:nrow(single_thread_results_of_large_instances)],:]

f = plt[:figure](figsize=(20, 10))
ax = sns.boxplot(x=single_thread_results_of_large_instances[:Algorithm], y=single_thread_results_of_large_instances[:Total_Time])
ax[:set_yscale]("log")
xlabel("Algorithms")
ylabel("Computing Time (secs)")
title("Large Instances")
f[:savefig]("../plots/1_Thread_Large_Instances_Computing_Time.eps")

f = plt[:figure](figsize=(20, 10))
ax = sns.boxplot(x=single_thread_results_of_large_instances[:Algorithm], y=single_thread_results_of_large_instances[:Ratio_to_Minimum_Time])
ax[:set_yscale]("log")
xlabel("Algorithms")
ylabel("Ratio to Minimum Time")
title("Large Instances")
f[:savefig]("../plots/1_Thread_Large_Instances_Ratio_to_Minimum_Time.eps")

f = plt[:figure](figsize=(20, 10))
sns.boxplot(x=single_thread_results_of_large_instances[:Algorithm], y=single_thread_results_of_large_instances[:Rank_of_the_Algorithm])
xlabel("Algorithms")
ylabel("Rank of the Algorithm")
title("Large Instances")
f[:savefig]("../plots/1_Thread_Large_Instances_Rank_of_the_Algorithm.eps")

f = plt[:figure](figsize=(20, 10))
sns.boxplot(x=single_thread_results_of_large_instances[:Algorithm], y=single_thread_results_of_large_instances[:IP_Solved])
xlabel("Algorithms")
ylabel("IPs Solved")
title("Large Instances")
f[:savefig]("../plots/1_Thread_Large_Instances_IP_Solved.eps")

f = plt[:figure](figsize=(20, 10))
sns.boxplot(x=single_thread_results_of_large_instances[:Algorithm], y=single_thread_results_of_large_instances[:Percent_of_Non_Dom_Pts_Found])
xlabel("Algorithms")
ylabel("Percentage of Nondominated Points Found")
title("Large Instances")
f[:savefig]("../plots/1_Thread_Large_Instances_Percent_of_Non_Dom_Pts_Found.eps")
