nBaseline = function(fullDF, probDF, rankedW){
  retval = data.frame()
  for (i in 1:nrow(probDF)){
    newRow = probDF[i,]
    sample = rownames(probDF)[i]
    nonZeros = colnames(fullDF)[which(fullDF[rownames(fullDF) == sample,] != 0)]
    numerator = (floor((ncol(fullDF)/5)/10)*10)
    n = numerator/(sum(match(nonZeros, rownames(summedW))))
    clus = newRow[1]
    newRow = newRow[2:6]
    newRow = newRow+n
    newRow = t(apply(newRow, 1, function(x) x/sum(x)))
    maxval = max(newRow)
    newRow = cbind(newRow,maxval)
    newRow = cbind(clus, newRow)
    retval = rbind(retval, newRow)
  }
  #fix inconsistencies
  maxvals = c()
  for(i in 1:nrow(retval)){
    pred = which(retval[i,2:6] == max(retval[i,2:6]))
    actual = retval[i,'cluster']
    if(pred != actual){
      print(paste("Sample",rownames(retval)[i])) 
      retval[i,(actual+1)] = max(retval[i,2:6]) + .01
      retval[i,2:6] = apply(retval[i,2:6],1,function(x) x/sum(x))
    }
    maxvals = c(maxvals, max(retval[i,2:6]))
  }
  retval$maxval = maxvals
  return(data.frame(retval))
}

nBaselineReduced = function(reducedDF, probDF){
  retval = data.frame()
  probDF = probDF[rownames(probDF) %in% rownames(reducedDF),]
  for (i in 1:nrow(probDF)){
    newRow = probDF[i,]
    sample = rownames(probDF)[i]
    currRow = reducedDF[rownames(reducedDF) == sample,2:(ncol(reducedDF)-1)]
    n = sum(reducedDF[rownames(reducedDF) == sample,2:(ncol(reducedDF)-1)] > 0)
    clus = newRow[1]
    newRow = newRow + (1/(2*n**2))
    newRow = t(apply(newRow[2:6], 1, function(x) x/sum(x)))
    maxval = max(newRow)
    newRow = cbind(newRow,maxval)
    newRow = cbind(clus, newRow)
    retval = rbind(retval, newRow)
  }
  return(data.frame(retval))
}

#tmpDF = nBaseline(fullDF, connectivityProbs1.0, summedW)
