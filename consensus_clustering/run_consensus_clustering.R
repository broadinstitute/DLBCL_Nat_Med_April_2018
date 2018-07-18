source('src/GDAC_TopgenesforCluster/Topgenes_v1.R')
result <- main("-s./src/GDAC_TopgenesforCluster/","-minput_data/DLBCL_significant_event_matrix.txt","-uALL","-ooutput_dir/DLBCL")


source("src/GDAC_NmfConsensusClustering/GDAC_CNMF.R", echo = TRUE)
result <- main("-ssrc/GDAC_NmfConsensusClustering/","-moutput_dir/DLBCL.expclu.gct","-ooutput_dir/DLBCL",'-u4','-v10') 
source("src/GDAC_selectBestcluster/select_best_cluster_chc.R")
result <- main("-uoutput_dir/DLBCL.expclu.gct","-mPearson","-coutput_dir/DLBCL.cophenetic.coefficient.txt","-woutput_dir/DLBCL.membership.txt","-voutput_dir/DLBCL","-pinput_data/DLBCL_significant_event_matrix.txt")
