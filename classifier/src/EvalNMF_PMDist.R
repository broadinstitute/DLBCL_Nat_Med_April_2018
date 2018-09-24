#A truly awful script that calls EvalNMF.Rmd many times.......
source("src/LoadLibraries.R")

generateNMF_Dist = function(nIter){
  
  for(k in 1:nIter){
    print(k)
    set.seed(k*10)
    source("src/EvalNMF.R")
    fn = "WrittenFiles/NMF_Dist.txt"
    if(!useSV && !useCNA){
      fn = "WrittenFiles/NMF_noSV_noCNA_Dist.txt"
    } else if(!useSV){
      fn = "WrittenFiles/NMF_noSV_Dist.txt"
    } else if(!useCNA){
      fn = "WrittenFiles/NMF_noCNA_Dist.txt"
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
  return(list(useSV, useCNA))
}


vals = generateNMF_Dist(200)
useSV = vals[[1]]
useCNA = vals[[2]]
fn = "WrittenFiles/NMF_Dist.txt"
if(!useSV && !useCNA){
  fn = "WrittenFiles/NMF_noSV_noCNA_Dist.txt"
} else if(!useSV){
  fn = "WrittenFiles/NMF_noSV_Dist.txt"
} else if(!useCNA){
  fn = "WrittenFiles/NMF_noCNA_Dist.txt"
}
distributionTable = read.csv(fn, header=FALSE, sep=',')
colnames(distributionTable) = c("weightedPearson", "accuracy", "performance")
p = hist(distributionTable$performance, breaks = 20, xlim = c(0,1))
fn = "Plots/NMF_dist.jpeg"
if(!useSV && !useCNA){
  fn = "Plots/NMF_noSV_noCNA_Dist.jpeg"
} else if(!useSV){
  fn = "Plots/NMF_noSV_Dist.jpeg"
} else if(!useCNA){
  fn = "Plots/NMF_noCNA_Dist.jpeg"
}
jpeg(filename = fn, width = 1080, height = 720)
hist(distributionTable$performance, breaks = 20, xlim = c(0,1), axes=FALSE)
axis(2)
axis(1, at=seq(0,1,0.05))
dev.off()
