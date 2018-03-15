source("src/GDAC_selectBestcluster/select_best_cluster_chc.R")
result <- main("-uoutput/DLBCL.expclu.gct","-mPearson","-coutput/DLBCL.cophenetic.coefficient.txt","-woutput/DLBCL.membership.txt","-voutput/DLBCL","-pinput_data/DLBCL_significant_event_matrix.txt")
