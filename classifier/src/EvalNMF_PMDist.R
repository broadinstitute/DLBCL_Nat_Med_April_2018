#A truly awful script that calls EvalNMF.Rmd many times.......
source("src/LoadLibraries.R")

generateNMF_Dist = function(nIter){
  
  for(k in 1:nIter){
    print(k)
    set.seed(k*10)
    source("src/EvalNMF.R")
    
    fn = "WrittenFiles/NMF_Dist"
    if(!useSV){
      fn = paste(fn, "_noSV", sep="")
    } 
    if(!useCNA){
      fn = paste(fn, "_noCNA", sep="")
    }
    if(GD){
      fn = paste(fn, "_GD", sep="")
    }
    if(fixedValidationSets){
      fn = paste(fn, "_FixedSets", sep="")
    }
    if(fullFeatures){
      fn = paste(fn,"_161F", sep="")
    }
    fn = paste(fn, ".txt", sep="")
    
    accuracy = top70thConfMat[[3]][[1]]
    
    fn2 = "WrittenFiles/NMF/Table.NMF_Dist"
    if(!useSV){
      fn2 = paste(fn2, "_noSV", sep="")
    }
    if(!useCNA){
      fn2 = paste(fn2, "_noCNA", sep="")
    }
    if(GD){
      fn2 = paste(fn2, "_GD", sep="")
    }
    if(fixedValidationSets){
      fn2 = paste(fn2, "_FixedSets", sep="")
    }
    if(fullFeatures){
      fn2 = paste(fn2,"_161F", sep="")
    }
    fn2 = paste(fn2,"_iteration",k, ".txt", sep="")
    
    nmfDF = data.frame(cbind(rownames(Hnorms), Hnorms$maxVal, Hnorms$correct, Hnorms$cluster, Hnorms$trueCluster))
    colnames(nmfDF) = c("Sample","maxVal","correct","predictedCluster","trueCluster")
    
    if(k != 1){
      currTable = read.csv(fn, header = FALSE, sep = ',')
      currTable = rbind(currTable, c(weightedPearson, accuracy, weightedPearson*accuracy))
      write.table(currTable, fn, row.names = FALSE, col.names = FALSE, sep = ',')
      write.table(nmfDF, fn2, row.names = FALSE, col.names = TRUE, sep=",", quote=FALSE)
    } else {
      currTable = t(c(weightedPearson, accuracy, weightedPearson*accuracy))
      write.table(currTable, fn, row.names = FALSE, col.names = FALSE, sep = ',')
      write.table(nmfDF, fn2, row.names = FALSE, col.names = TRUE, sep=",", quote=FALSE)
    }
    
  }
  return(list(useSV, useCNA, GD))
}


vals = generateNMF_Dist(100)
useSV = vals[[1]]
useCNA = vals[[2]]
GD = vals[[3]]

fn = "WrittenFiles/NMF_Dist"
if(!useSV){
  fn = paste(fn, "_noSV", sep="")
} 
if(!useCNA){
  fn = paste(fn, "_noCNA", sep="")
}
if(GD){
  fn = paste(fn, "_GD", sep="")
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
p = hist(distributionTable$performance, breaks = 20, xlim = c(0,1))
fn = "Plots/NMF_dist"

if(!useSV){
  fn = paste(fn, "_noSV", sep="")
} 
if(!useCNA){
  fn = paste(fn, "_noCNA", sep="")
}
if(GD){
  fn = paste(fn, "_GD", sep="")
}
if(fixedValidationSets){
  fn = paste(fn, "_FixedSets", sep="")
}
if(fullFeatures){
  fn = paste(fn,"_161F", sep="")
}
fn = paste(fn, ".jpeg", sep="")
jpeg(filename = fn, width = 1080, height = 720)
hist(distributionTable$performance, breaks = 20, xlim = c(0,1), axes=FALSE)
axis(2)
axis(1, at=seq(0,1,0.05))
dev.off()
