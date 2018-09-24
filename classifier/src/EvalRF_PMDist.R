generateRF_Dist = function(nIter){
  for(k in 1:nIter){
    print(k)
    set.seed(k)
    if((k %% (nIter)/10) == 0){
      print(k)
    }
    source("src/RF_TrainNonExpand_Reduced.R")
    fn = "WrittenFiles/RF_Dist.tmp"
    if(!useSV && !useCNA){
      fn = "WrittenFiles/RF_noSV_noCNA_Dist.tmp"
    } else if(!useSV){
      fn = "WrittenFiles/RF_noSV_Dist.tmp"
    } else if(!useCNA){
      fn = "WrittenFiles/RF_noCNA_Dist.tmp"
    }
    accuracy = top70thConfMat[[3]][[1]]
    
    if(k != 1){
      currTable = read.csv(fn, header = FALSE, sep = ',')
      currTable = rbind(currTable, c(weightedPearson, accuracy, weightedPearson*accuracy))
      write.table(currTable, fn, row.names = FALSE, col.names = FALSE, sep = ',')
    } else {
      currTable = t(c(weightedPearson, accuracy, weightedPearson*accuracy))
      write.table(currTable, fn, row.names = FALSE, col.names = FALSE, sep = ',')
    }
  }
}


generateRF_Dist(200)
fn = "WrittenFiles/RF_Dist.tmp"
if(!useSV && !useCNA){
  fn = "WrittenFiles/RF_noSV_noCNA_Dist.tmp"
} else if(!useSV){
  fn = "WrittenFiles/RF_noSV_Dist.tmp"
} else if(!useCNA){
  fn = "WrittenFiles/RF_noCNA_Dist.tmp"
}
distributionTable = read.csv(fn, header=FALSE, sep=',')
colnames(distributionTable) = c("weightedPearson", "accuracy", "performance")
p = hist(distributionTable$performance, breaks = 20, xlim = c(0,1))
fn = "Plots/RF_dist.jpeg"
if(!useSV && !useCNA){
  fn = "Plots/RF_noSV_noCNA_Dist.jpeg"
} else if(!useSV){
  fn = "Plots/RF_noSV_Dist.jpeg"
} else if(!useCNA){
  fn = "Plots/RF_noCNA_Dist.jpeg"
}
# jpeg(filename = fn, width = 1080, height = 720)
# hist(distributionTable$performance, breaks = 20, xlim = c(0,1), axes=FALSE)
# axis(2)
# axis(1, at=seq(0,1,0.05))
# dev.off()
