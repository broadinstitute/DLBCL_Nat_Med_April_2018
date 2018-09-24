#Need to sync the two arguments in the system() call with the two variables in the EvalNN.R file. useSV is first argument,
#useCNA is second argument.

source("src/LoadLibraries.R")

generateNN_Dist = function(nIter){
  
  for(k in 1:nIter){
    set.seed(k)
    print(k)
    #this line needs to point to the interpreter that torch was installed for
    setwd("ANN")
    system("/Users/twood/anaconda3/bin/python testFoldsNNs.py True False")
    setwd("..")
    source("src/EvalNN.R")
    fn = "WrittenFiles/NN_Dist.txt"
    if(!useSV && !useCNA){
      fn = "WrittenFiles/NN_noSV_noCNA_Dist.txt"
    } else if(!useSV){
      fn = "WrittenFiles/NN_noSV_Dist.txt"
    } else if(!useCNA){
      fn = "WrittenFiles/NN_noCNA_Dist.txt"
    }

    accuracy = top70thConfMat[[3]][[1]]
    
    if(k != 1){
      currTable = read.csv(fn, header = FALSE, sep = ',')
      currTable = rbind(currTable, c(weightedPearson, accuracy, weightedPearson*accuracy))
      #print(currTable)
      write.table(currTable, fn, row.names = FALSE, col.names = FALSE, sep = ',')
    } else {
      currTable = t(c(weightedPearson, accuracy, weightedPearson*accuracy))
      write.table(currTable, fn, row.names = FALSE, col.names = FALSE, sep = ',')
    }
    
  }
}


generateNN_Dist(200)
fn = "WrittenFiles/NN_Dist.txt"
if(!useSV && !useCNA){
  fn = "WrittenFiles/NN_noSV_noCNA_Dist.txt"
} else if(!useSV){
  fn = "WrittenFiles/NN_noSV_Dist.txt"
} else if(!useCNA){
  fn = "WrittenFiles/NN_noCNA_Dist.txt"
}
distributionTable = read.csv(fn, header=FALSE, sep=',')
colnames(distributionTable) = c("weightedPearson", "accuracy", "performance")
hist(distributionTable$performance)
p = hist(distributionTable$performance, breaks = 20, xlim = c(0,1))
fn = "Plots/NN_dist.jpeg"
if(!useSV && !useCNA){
  fn = "Plots/NN_noSV_noCNA_Dist.jpeg"
} else if(!useSV){
  fn = "Plots/NN_noSV_Dist.jpeg"
} else if(!useCNA){
  fn = "Plots/NN_noCNA_Dist.jpeg"
}
# jpeg(fn, width = 1080, height = 720)
# hist(distributionTable$performance, breaks = 20, xlim = c(0,1), axes=FALSE)
# axis(2)
# axis(1, at=seq(0,1,0.05))
# dev.off()
