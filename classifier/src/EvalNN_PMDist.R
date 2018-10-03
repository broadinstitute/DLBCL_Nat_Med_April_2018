#Need to sync the two arguments in the system() call with the two variables in the EvalNN.R file. useSV is first argument,
#useCNA is second argument, genome doubling the third, use probability labels 4th, fixed validation set 5th,
#161 features 6th argument

source("src/LoadLibraries.R")

generateNN_Dist = function(nIter){
  for(k in 1:nIter){
    set.seed(k)
    print(k)
    #this line needs to point to the interpreter that torch was installed for
    setwd("ANN")
    system("/Users/twood/anaconda3/bin/python testFoldsNNs.py True True False True False True")
    setwd("..")
    source("src/EvalNN.R")
    fn = "WrittenFiles/NN_Dist"
    if(!useSV){
      fn = paste(fn, "_noSV", sep="")
    }
    if(!useCNA){
      fn = paste(fn, "_noCNA", sep="")
    }
    if(GD){
      fn = paste(fn, "_GD", sep="")
    }
    if(!useProbs){
      fn = paste(fn, "_BinaryLabels", sep="")
    }
    if(fixedValidationSets){
      fn = paste(fn, "_FixedSets", sep="")
    }
    if(fullFeatures){
      fn = paste(fn,"_161F", sep="")
    }
    fn = paste(fn, ".txt", sep="")
    
    fn2 = "WrittenFiles/ANN/Table.NN_Dist"
    if(!useSV){
      fn2 = paste(fn2, "_noSV", sep="")
    }
    if(!useCNA){
      fn2 = paste(fn2, "_noCNA", sep="")
    }
    if(GD){
      fn2 = paste(fn2, "_GD", sep="")
    }
    if(!useProbs){
      fn2 = paste(fn2, "_BinaryLabels", sep="")
    }
    if(fixedValidationSets){
      fn2 = paste(fn2, "_FixedSets", sep="")
    }
    if(fullFeatures){
      fn2 = paste(fn2,"_161F", sep="")
    }
    fn2 = paste(fn2,"_iteration",k, ".txt", sep="")
    
    accuracy = top70thConfMat[[3]][[1]]
    
    if(k != 1){
      currTable = read.csv(fn, header = FALSE, sep = ',')
      currTable = rbind(currTable, c(weightedPearson, accuracy, weightedPearson*accuracy))
      #print(currTable)
      write.table(currTable, fn, row.names = FALSE, col.names = FALSE, sep = ',')
      write.table(nnDF, fn2, row.names = TRUE, col.names = TRUE, sep=',')
    } else {
      currTable = t(c(weightedPearson, accuracy, weightedPearson*accuracy))
      write.table(currTable, fn, row.names = FALSE, col.names = FALSE, sep = ',')
      write.table(nnDF, fn2, row.names = TRUE, col.names = TRUE, sep=',')
    }
  }
}


generateNN_Dist(100)
fn = "WrittenFiles/NN_Dist"
if(!useSV){
  fn = paste(fn, "_noSV", sep="")
} 
if(!useCNA){
  fn = paste(fn, "_noCNA", sep="")
}
if(GD){
  fn = paste(fn, "_GD", sep="")
}
if(!useProbs){
  fn = paste(fn, "_BinaryLabels", sep="")
}
if(fixedValidationSets){
  fn = paste(fn, "_FixedSets", sep="")
}
if(fullFeatures){
  fn = paste(fn,"_161F", sep="")
}
fn = paste(fn, ".txt", sep="")
distributionTable = read.csv(fn, header=FALSE, sep=',')
colnames(distributionTable) = c("weightedPearson", "accuracy", "performance")
hist(distributionTable$performance)
p = hist(distributionTable$performance, breaks = 20, xlim = c(0,1))
fn = "Plots/NN_dist"
if(!useSV){
  fn = paste(fn, "_noSV", sep="")
} 
if(!useCNA){
  fn = paste(fn, "_noCNA", sep="")
}
if(GD){
  fn = paste(fn, "_GD", sep="")
}
if(!useProbs){
  fn = paste(fn, "_BinaryLabels", sep="")
}
if(fixedValidationSets){
  fn = paste(fn, "_FixedSets", sep="")
}
if(fullFeatures){
  fn = paste(fn,"_161F", sep="")
}
fn = paste(fn, ".jpeg", sep="")
# jpeg(fn, width = 1080, height = 720)
# hist(distributionTable$performance, breaks = 20, xlim = c(0,1), axes=FALSE)
# axis(2)
# axis(1, at=seq(0,1,0.05))
# dev.off()
