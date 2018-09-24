probabilityTable = function(predFrame){
  retval = c()
  for (i in 1:nrow(predFrame)){
    retval = c(retval,max(table(predFrame[i,]))/length(predFrame[1,]))
  }
  return(retval)
}