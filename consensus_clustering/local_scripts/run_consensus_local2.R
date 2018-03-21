source("src/GDAC_selectBestcluster/select_best_cluster_chc.R")
result <- main("-uDLBCL.expclu.gct","-mPearson","-cDLBCL.cophenetic.coefficient.txt","-wDLBCL.membership.txt","-vDLBCL","-pinput_data/DLBCL_significant_event_matrix.txt")
