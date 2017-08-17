# Supplemental Materials and Experimental Results for "Computational Study of Biobjective Pure Integer Linear Programs" #

## Julia sourcecode: ##

1. [Experimental_Run.jl](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/src/Experimental_Run.jl) - Script to run the whole experiment. The proper way to run the Experiments will be `julia -p 4 Experimental_Run.jl`. In this case, experiments related to MDLS, V1, V2 and V3 will run 4 instances in parallel. However, all experiments related to NSGA-II and V3 with 1, 2, 3 and 4 on the large random instances will be run serially.
2. [Summarizing_Results.jl](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/src/Summarizing_Results.jl) - Script to generate all the plots from the [Experimental Results](https://github.com/aritrasep/MOO_JA_2_Sup.jl/blob/master/results/Experimental_Results.csv).

## Experimental Results: ##

1. [Detailed Experimental Results](https://github.com/aritrasep/MOO_JA_2_Sup.jl/blob/master/results/Experimental_Results.csv)
2. Nondominated Frontiers:
	1. [NSGA-II](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/NSGA-II/)
	2. [MDLS](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/MDLS/)
	3. [V1](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/V1/)
	4. [V2](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/V2/)
	5. [V3](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/V3/)
	6. Random Instances:
		1. [Combined Frontier](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/BOKP_Aritra_Instances/Combined_Frontier/)
		2. [V3 - 1 Thread](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/BOKP_Aritra_Instances/Thread_1/)
		3. [V3 - 2 Threads](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/BOKP_Aritra_Instances/Thread_2/)
		4. [V3 - 3 Threads](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/BOKP_Aritra_Instances/Thread_3/)
		5. [V3 - 4 Threads](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/BOKP_Aritra_Instances/Thread_4/)
3. [Plots](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/plots/)
