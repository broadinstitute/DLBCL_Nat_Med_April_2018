source("src/ComputeErrors.R")

ensembleErrors = function(predictedDF, actual, verb){
  predictions = c()
  for(i in 1:nrow(predictedDF)){
    predictions = c(predictions, as.integer(names(which(table(predictedDF[i,]) == max(table(predictedDF[i,]))))[[1]]))
  }
  if(verb){
    print(predictions)
    print(actual)
  }
  confmat = computeErrors(predictions, actual)
  errors = computeClassErrors(confmat,verbose=verb)
  return(list(confmat,errors, predictions))
}