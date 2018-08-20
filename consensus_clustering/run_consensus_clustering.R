source('src/GDAC_TopgenesforCluster/Topgenes_v1.R')
result <- main("-s./src/GDAC_TopgenesforCluster/","-minput_data/DLBCL_significant_event_matrix.txt","-uALL","-oDLBCL")


source("src/GDAC_NmfConsensusClustering/GDAC_CNMF.R", echo = TRUE)
result <- main("-ssrc/GDAC_NmfConsensusClustering/","-mDLBCL.expclu.gct","-oDLBCL",'-u4','-v10') 
source("src/GDAC_selectBestcluster/select_best_cluster_chc.R")
result <- main("-uDLBCL.expclu.gct","-mPearson","-cDLBCL.cophenetic.coefficient.txt","-wDLBCL.membership.txt","-voutput_dir/DLBCL","-pinput_data/DLBCL_significant_event_matrix.txt")
