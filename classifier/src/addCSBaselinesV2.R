
featureSum = function(featureVec, fsize){
  if(length(featureVec) == 0){
    return(0)
  }
  retval = 0
  for(i in 1:length(featureVec)){
    toAdd = (fsize-featureVec[i]+1)
    retval = retval+toAdd
  }
  return(retval)
}

fixInconsistencies = function(probDF){
  retval = probDF
  for(i in 1:nrow(retval)){
    maxidx = which(retval[i,2:6] == max(retval[i,2:6]))
    correctCluster = retval[i,'cluster']
    if(correctCluster != maxidx){
      retval[i,(correctCluster+1)] = max(retval[i,2:6])+.01
      retval[i,2:6] = t(apply(retval[i,2:6],1,function(x) x/sum(x)))
    }
  }
  return(retval)
}

computeRanks = function(fullDF, probDF, labels, summedW){
  clus1Order = c("SV.BCL6", "BCL10", "TNFAIP3", "UBE2A", "CD70", "B2M", "NOTCH2", "TMEM30A", "FAS",
                 "X5p.AMP", "SV.TP63", "ZEB2", "HLA.B", "SPEN", "SV.CD274.PDCD1LG2")
  clus3Order = c("TP53", "X17p.DEL", "X17p11.2.DEL", "X21q.AMP", "X9p21.3.DEL", "X9q21.13.DEL",
                 "X4q35.1.DEL","X1p31.1.DEL", "X1p36.11.DEL", "X1p13.1.DEL", "X4q21.22.DEL", "X14q32.31.DEL",
                 "X3p21.31.DEL", "X2p16.1.AMP", "X16q12.1.DEL", "X1p36.32.DEL", "X3q28.DEL", "X1q23.3.AMP" ,
                 "X18q23.DEL", "X8q24.22.AMP", "X17q24.3.AMP", "X13q14.2.DEL", "X19p13.3.DEL", "X5q.AMP",
                 "X11q.AMP", "X13q31.3.AMP", "X6p.AMP", "X2q22.2.DEL", "X12p13.2.DEL", "X6q.DEL",
                 "X3q28.AMP", "X11q23.3.AMP", "X1q42.12.DEL", "X8q12.1.DEL", "X19q13.32.DEL", "X10q23.31.DEL")
  
  clus2Order = c("BCL2", "SV.BCL2", "CREBBP", "EZH2", "KMT2D", "TNFRSF14", "HVCN1", "IRF8", "GNA13", "MEF2B",
                 "PTEN")
  
  clus5Order = c("SGK1", "HIST1H1E", "NFKBIE","BRAF", "CD83", "NFKBIA", "CD58", "HIST1H2BC", "STAT3",
                 "HIST1H1C", "ZFP36L1", "KLHL6", "HIST1H1D", "HIST1H1B", "ETS1", "TOX", "HIST1H2AM",
                 "HIST1H2BK", "RHOA", "ACTB", "LTB", "SF3B1", "CARD11", "HIST1H2AC")
  
  clus4Order = c("X18q.AMP", "X13q.AMP", "CD79B", "X3p.AMP", "MYD88", "ETV6", "X18p.AMP", "PIM1", 
                 "X17q25.1.DEL", "TBL1XR1", "X19q13.42.AMP", "GRHPR", "ZC3H12A", "X19p13.2.DEL",
                 "X19q.AMP", "HLA.A", "PRDM1", "BTG1", "X18q21.33.BCL2..AMP", "SV.MYC")
  
  clus1Order = toupper(clus1Order)
  clus2Order = toupper(clus2Order)
  clus3Order = toupper(clus3Order)
  clus4Order = toupper(clus4Order)
  clus5Order = toupper(clus5Order)
  
  c1DF = fullDF[,clus1Order]
  c2DF = fullDF[,clus2Order]
  c3DF = fullDF[,clus3Order]
  c4DF = fullDF[,clus4Order]
  c5DF = fullDF[,clus5Order]
  
  subsetLabels = labels[rownames(labels) %in% rownames(fullDF),]
  
  rankedDF1 = data.frame()
  rankedDF2 = data.frame()
  rankedDF3 = data.frame()
  rankedDF4 = data.frame()
  rankedDF5 = data.frame()
  n1s = c()
  
  lenC1 = 1:ncol(c1DF)
  lenC2 = 1:ncol(c2DF)
  lenC3 = 1:ncol(c3DF)
  lenC4 = 1:ncol(c4DF)
  lenC5 = 1:ncol(c5DF)
  mainPower = 1
  alpha = .90
  
  numerator = (floor((ncol(fullDF)/5)/10)*10)
  denomCoeff = 1
  denomPower = 1.2
  constantScalar = 1/2
  for(i in 1:nrow(subsetLabels)){
    currCluster = subsetLabels[i,'cluster']
    sample = rownames(subsetLabels)[i]
    if (currCluster == 1){
      fullNonZeros = colnames(fullDF)[which(fullDF[rownames(fullDF) == sample,] != 0)]
      n = numerator/(denomCoeff*sum(match(fullNonZeros, rownames(summedW))**denomPower))
      
      #compute the ratio of ranks for the assigned cluster
      nonZerosMain = colnames(c1DF)[which(c1DF[rownames(c1DF) == sample,] != 0)]
      
      mainRanks = match(nonZerosMain, colnames(c1DF))
      mainRanks = ncol(c1DF)-mainRanks+1
      
      
      mainRanks = mainRanks**mainPower
      term1 = sum(mainRanks)/(sum(lenC1**mainPower))
      
      
      #compute the ratio of ranks for each unassigned cluster
      ranksA = colnames(c2DF)[which(c2DF[rownames(c2DF) == sample,] != 0)]
      ranksA = match(ranksA, colnames(c2DF))
      ranksA = sum(ranksA)
      partA = sum(lenC2)
      
      ranksB = colnames(c3DF)[which(c3DF[rownames(c3DF) == sample,] != 0)]
      ranksB = match(ranksB, colnames(c3DF))
      ranksB = sum(ranksB)
      partB = sum(lenC3)
      
      ranksC = colnames(c4DF)[which(c4DF[rownames(c4DF) == sample,] != 0)]
      ranksC = match(ranksC, colnames(c4DF))
      ranksC = sum(ranksC)
      partC = sum(lenC4)
      
      ranksD = colnames(c5DF)[which(c5DF[rownames(c5DF) == sample,] != 0)]
      ranksD = match(ranksD, colnames(c5DF))
      ranksD = sum(ranksD)
      partD = sum(lenC5)
      
      term2 = (sum(ranksA,ranksB,ranksC,ranksD)/sum(partA,partB,partC,partB))
      
      #if you have no marker features for your assigned cluster, then you autolose
      epsilon = term2/term1
      if(term1 == 0){
        epsilon = 10
      }
      
      toadd = constantScalar*(alpha*n+(1-alpha)*epsilon)
      rankedDF1 = rbind(rankedDF1, toadd)
      rownames(rankedDF1)[nrow(rankedDF1)] = sample
      
    } else if (currCluster == 2){
      fullNonZeros = colnames(fullDF)[which(fullDF[rownames(fullDF) == sample,] != 0)]
      n = numerator/(denomCoeff*sum(match(fullNonZeros, rownames(summedW))**denomPower))
      
      #compute the ratio of ranks for the assigned cluster
      nonZerosMain = colnames(c2DF)[which(c2DF[rownames(c2DF) == sample,] != 0)]
      
      mainRanks = match(nonZerosMain, colnames(c2DF))
      mainRanks = ncol(c2DF)-mainRanks+1
      
      
      mainRanks = mainRanks**mainPower
      term1 = sum(mainRanks)/(sum(lenC2**mainPower))
      
      
      #compute the ratio of ranks for each unassigned cluster
      ranksA = colnames(c1DF)[which(c1DF[rownames(c1DF) == sample,] != 0)]
      ranksA = match(ranksA, colnames(c1DF))
      ranksA = sum(ranksA)
      partA = sum(lenC1)
      
      ranksB = colnames(c3DF)[which(c3DF[rownames(c3DF) == sample,] != 0)]
      ranksB = match(ranksB, colnames(c3DF))
      ranksB = sum(ranksB)
      partB = sum(lenC3)
      
      ranksC = colnames(c4DF)[which(c4DF[rownames(c4DF) == sample,] != 0)]
      ranksC = match(ranksC, colnames(c4DF))
      ranksC = sum(ranksC)
      partC = sum(lenC4)
      
      ranksD = colnames(c5DF)[which(c5DF[rownames(c5DF) == sample,] != 0)]
      ranksD = match(ranksD, colnames(c5DF))
      ranksD = sum(ranksD)
      partD = sum(lenC5)
      
      term2 = (sum(ranksA,ranksB,ranksC,ranksD)/sum(partA,partB,partC,partB))
      
      #if you have no marker features for your assigned cluster, then you autolose
      epsilon = term2/term1
      if(term1 == 0){
        epsilon = 10
      }
      
      toadd = constantScalar*(alpha*n+(1-alpha)*epsilon)
      rankedDF2 = rbind(rankedDF2, toadd)
      rownames(rankedDF2)[nrow(rankedDF2)] = sample
      
      
    } else if (currCluster == 3){
      fullNonZeros = colnames(fullDF)[which(fullDF[rownames(fullDF) == sample,] != 0)]
      n = numerator/(denomCoeff*sum(match(fullNonZeros, rownames(summedW))**denomPower))
      
      #compute the ratio of ranks for the assigned cluster
      nonZerosMain = colnames(c3DF)[which(c3DF[rownames(c3DF) == sample,] != 0)]
      
      mainRanks = match(nonZerosMain, colnames(c3DF))
      mainRanks = ncol(c3DF)-mainRanks+1
      
      
      mainRanks = mainRanks**mainPower
      term1 = sum(mainRanks)/(sum(lenC3**mainPower))
      
      
      #compute the ratio of ranks for each unassigned cluster
      ranksA = colnames(c1DF)[which(c1DF[rownames(c1DF) == sample,] != 0)]
      ranksA = match(ranksA, colnames(c1DF))
      ranksA = sum(ranksA)
      partA = sum(lenC1)
      
      ranksB = colnames(c2DF)[which(c2DF[rownames(c2DF) == sample,] != 0)]
      ranksB = match(ranksB, colnames(c2DF))
      ranksB = sum(ranksB)
      partB = sum(lenC2)
      
      ranksC = colnames(c4DF)[which(c4DF[rownames(c4DF) == sample,] != 0)]
      ranksC = match(ranksC, colnames(c4DF))
      ranksC = sum(ranksC)
      partC = sum(lenC4)
      
      ranksD = colnames(c5DF)[which(c5DF[rownames(c5DF) == sample,] != 0)]
      ranksD = match(ranksD, colnames(c5DF))
      ranksD = sum(ranksD)
      partD = sum(lenC5)
      
      term2 = (sum(ranksA,ranksB,ranksC,ranksD)/sum(partA,partB,partC,partB))
      
      #if you have no marker features for your assigned cluster, then you autolose
      epsilon = term2/term1
      if(term1 == 0){
        epsilon = 10
      }
      
      toadd = constantScalar*(alpha*n+(1-alpha)*epsilon)
      rankedDF3 = rbind(rankedDF3, toadd)
      rownames(rankedDF3)[nrow(rankedDF3)] = sample
      
      
    } else if (currCluster == 4){
      fullNonZeros = colnames(fullDF)[which(fullDF[rownames(fullDF) == sample,] != 0)]
      n = numerator/(denomCoeff*sum(match(fullNonZeros, rownames(summedW))**denomPower))
      
      #compute the ratio of ranks for the assigned cluster
      nonZerosMain = colnames(c4DF)[which(c4DF[rownames(c4DF) == sample,] != 0)]
      
      mainRanks = match(nonZerosMain, colnames(c4DF))
      mainRanks = ncol(c4DF)-mainRanks+1
      
      
      mainRanks = mainRanks**mainPower
      term1 = sum(mainRanks)/(sum(lenC4**mainPower))
      
      
      #compute the ratio of ranks for each unassigned cluster
      ranksA = colnames(c1DF)[which(c1DF[rownames(c1DF) == sample,] != 0)]
      ranksA = match(ranksA, colnames(c1DF))
      ranksA = sum(ranksA)
      partA = sum(lenC1)
      
      ranksB = colnames(c3DF)[which(c3DF[rownames(c3DF) == sample,] != 0)]
      ranksB = match(ranksB, colnames(c3DF))
      ranksB = sum(ranksB)
      partB = sum(lenC3)
      
      ranksC = colnames(c2DF)[which(c2DF[rownames(c2DF) == sample,] != 0)]
      ranksC = match(ranksC, colnames(c2DF))
      ranksC = sum(ranksC)
      partC = sum(lenC2)
      
      ranksD = colnames(c5DF)[which(c5DF[rownames(c5DF) == sample,] != 0)]
      ranksD = match(ranksD, colnames(c5DF))
      ranksD = sum(ranksD)
      partD = sum(lenC5)
      
      term2 = (sum(ranksA,ranksB,ranksC,ranksD)/sum(partA,partB,partC,partB))
      
      #if you have no marker features for your assigned cluster, then you autolose
      epsilon = term2/term1
      if(term1 == 0){
        epsilon = 10
      }
      
      toadd = constantScalar*(alpha*n+(1-alpha)*epsilon)
      rankedDF4 = rbind(rankedDF4, toadd)
      rownames(rankedDF4)[nrow(rankedDF4)] = sample
      
      
    } else if (currCluster == 5){
      fullNonZeros = colnames(fullDF)[which(fullDF[rownames(fullDF) == sample,] != 0)]
      n = numerator/(denomCoeff*sum(match(fullNonZeros, rownames(summedW))**denomPower))
      
      #compute the ratio of ranks for the assigned cluster
      nonZerosMain = colnames(c5DF)[which(c5DF[rownames(c5DF) == sample,] != 0)]
      
      mainRanks = match(nonZerosMain, colnames(c5DF))
      mainRanks = ncol(c5DF)-mainRanks+1
      
      
      mainRanks = mainRanks**mainPower
      term1 = sum(mainRanks)/(sum(lenC5**mainPower))
      
      
      #compute the ratio of ranks for each unassigned cluster
      ranksA = colnames(c1DF)[which(c1DF[rownames(c1DF) == sample,] != 0)]
      ranksA = match(ranksA, colnames(c1DF))
      ranksA = sum(ranksA)
      partA = sum(lenC1)
      
      ranksB = colnames(c3DF)[which(c3DF[rownames(c3DF) == sample,] != 0)]
      ranksB = match(ranksB, colnames(c3DF))
      ranksB = sum(ranksB)
      partB = sum(lenC3)
      
      ranksC = colnames(c4DF)[which(c4DF[rownames(c4DF) == sample,] != 0)]
      ranksC = match(ranksC, colnames(c4DF))
      ranksC = sum(ranksC)
      partC = sum(lenC4)
      
      ranksD = colnames(c2DF)[which(c2DF[rownames(c2DF) == sample,] != 0)]
      ranksD = match(ranksD, colnames(c2DF))
      ranksD = sum(ranksD)
      partD = sum(lenC2)
      
      term2 = (sum(ranksA,ranksB,ranksC,ranksD)/sum(partA,partB,partC,partB))
      
      #if you have no marker features for your assigned cluster, then you autolose
      epsilon = term2/term1
      if(term1 == 0){
        epsilon = 10
      }
      
      toadd = constantScalar*(alpha*n+(1-alpha)*epsilon)
      rankedDF5 = rbind(rankedDF5, toadd)
      rownames(rankedDF5)[nrow(rankedDF5)] = sample
    }
    n1s = c(n1s, alpha*n)
  }
  colnames(rankedDF1)[1] = 'val'
  colnames(rankedDF2)[1] = 'val'
  colnames(rankedDF3)[1] = 'val'
  colnames(rankedDF4)[1] = 'val'
  colnames(rankedDF5)[1] = 'val'
  return(list(rankedDF1,rankedDF2,rankedDF3,rankedDF4,rankedDF5, n1s))
}


csBaselines = function(fullDF, connectivityDF, labels, summedW){
  probPower = 1
  probDF = fixInconsistencies(connectivityDF)
  rankedDFs = computeRanks(fullDF, probDF, labels, summedW)
  
  c1DF = rankedDFs[[1]]
  c2DF = rankedDFs[[2]]
  c3DF = rankedDFs[[3]]
  c4DF = rankedDFs[[4]]
  c5DF = rankedDFs[[5]]
  n1s = rankedDFs[[6]]
  
  c1Ret = data.frame()
  c2Ret = data.frame()
  c3Ret = data.frame()
  c4Ret = data.frame()
  c5Ret = data.frame()
  
  for(i in 1:nrow(probDF)){
    currProbs = probDF[i,2:6]
    currClus = probDF[i,1]
    sample = rownames(probDF)[i]
    if(currClus == 1){
      val = c1DF[rownames(c1DF) == sample,'val']
      currProbs = (currProbs+val)**probPower
      currProbs = t(apply(currProbs,1,function(x) x/sum(x)))
      c1Ret = rbind(c1Ret,c(currProbs,currClus))
      rownames(c1Ret)[nrow(c1Ret)] = sample
    } else if(currClus == 2){
      val = c2DF[rownames(c2DF) == sample,'val']
      currProbs = (currProbs+val)**probPower
      currProbs = t(apply(currProbs,1,function(x) x/sum(x)))
      c2Ret = rbind(c2Ret,c(currProbs,currClus))
      rownames(c2Ret)[nrow(c2Ret)] = sample
    } else if(currClus == 3){
      val = c3DF[rownames(c3DF) == sample,'val']
      currProbs = (currProbs+val)**probPower
      currProbs = t(apply(currProbs,1,function(x) x/sum(x)))
      c3Ret = rbind(c3Ret,c(currProbs,currClus))
      rownames(c3Ret)[nrow(c3Ret)] = sample
    } else if(currClus == 4){
      val = c4DF[rownames(c4DF) == sample,'val']
      currProbs = (currProbs+val)**probPower
      currProbs = t(apply(currProbs,1,function(x) x/sum(x)))
      c4Ret = rbind(c4Ret,c(currProbs,currClus))
      rownames(c4Ret)[nrow(c4Ret)] = sample
    } else if(currClus == 5){
      val = c5DF[rownames(c5DF) == sample,'val']
      currProbs = (currProbs+val)**probPower
      currProbs = t(apply(currProbs,1,function(x) x/sum(x)))
      c5Ret = rbind(c5Ret,c(currProbs,currClus))
      rownames(c5Ret)[nrow(c5Ret)] = sample
    }
  }
  colnames(c1Ret) = c("P1","P2","P3","P4","P5","cluster")
  colnames(c2Ret) = c("P1","P2","P3","P4","P5","cluster")
  colnames(c3Ret) = c("P1","P2","P3","P4","P5","cluster")
  colnames(c4Ret) = c("P1","P2","P3","P4","P5","cluster")
  colnames(c5Ret) = c("P1","P2","P3","P4","P5","cluster")
  
  c1Ret = c1Ret[order(c1Ret$P1), , drop=FALSE]
  c2Ret = c2Ret[order(c2Ret$P2), , drop=FALSE]
  c3Ret = c3Ret[order(c3Ret$P3), , drop=FALSE]
  c4Ret = c4Ret[order(c4Ret$P4), , drop=FALSE]
  c5Ret = c5Ret[order(c5Ret$P5), , drop=FALSE]
  
  retval = rbind(c1Ret,c2Ret,c3Ret,c4Ret,c5Ret)
  return(cbind(n1s,retval))
}

#csBaselines(fullDF, connectivityProbs1.0, labels, summedW)