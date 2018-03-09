source('/xchip/cga_home/adunford/DLBCL/Consensus_Clustering/local_run/GDAC_TopgenesforCluster/Topgenes_v1.R')
result <- main("-s/xchip/cga_home/adunford/DLBCL/Consensus_Clustering/local_run/GDAC_TopgenesforCluster/","-m/xchip/cga_home/adunford/DLBCL/Gene_Sample_Matrix/DLBCL_mutation_scna_sv_matrix.no_21q22.3.21_July_2017.txt","-uALL","-oDLBCL_no_21q22.3")


#require(jit)
source("/xchip/cga_home/adunford/DLBCL/Consensus_Clustering/local_run/GDAC_NmfConsensusClustering/GDAC_CNMF.R", echo = TRUE)
result <- main("-s/xchip/cga_home/adunford/DLBCL/Consensus_Clustering/local_run/GDAC_NmfConsensusClustering/","-mDLBCL_no_21q22.3.expclu.gct","-oDLBCL_no_21q22.3",'-u4','-v10') 
