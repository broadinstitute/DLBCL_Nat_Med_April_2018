standardize = function(reducedDF){
  rown = rownames(reducedDF)
  coln = colnames(reducedDF)
  retval = matrix(0L, nrow = nrow(reducedDF), ncol = ncol(reducedDF))
  
  for(i in 1:(ncol(reducedDF))){
    currCol = as.numeric(as.character(reducedDF[,i]))
    newCol = (currCol-min(currCol))/(max(currCol) - min(currCol))
    retval[,i] = newCol
  }
  colnames(retval) = coln
  rownames(retval) = rown
  return(data.frame(retval))
}

