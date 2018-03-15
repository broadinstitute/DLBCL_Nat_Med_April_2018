source('src/GDAC_TopgenesforCluster/Topgenes_v1.R')
result <- main("-s./src/GDAC_TopgenesforCluster/","-minput_data/DLBCL_significant_event_matrix.txt","-uALL","-oDLBCL")


#require(jit)
source("src/GDAC_NmfConsensusClustering/GDAC_CNMF.R", echo = TRUE)
result <- main("-ssrc/GDAC_NmfConsensusClustering/","-mDLBCL.expclu.gct","-oDLBCL",'-u4','-v10') 
