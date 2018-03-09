source("/xchip/cga_home/adunford/DLBCL/Consensus_Clustering/local_run/GDAC_selectBestcluster/select_best_cluster_chc.R")
result <- main("-uDLBCL_no_21q22.3.expclu.gct","-mPearson","-cDLBCL_no_21q22.3.cophenetic.coefficient.txt","-wDLBCL_no_21q22.3.membership.txt","-vDLBCL_no_21q22.3","-p/xchip/cga_home/adunford/DLBCL/Gene_Sample_Matrix/DLBCL_mutation_scna_sv_matrix.no_21q22.3.21_July_2017.txt")
