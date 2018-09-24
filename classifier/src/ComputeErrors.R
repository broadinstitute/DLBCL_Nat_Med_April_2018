computeErrors = function(predicted,actual){
  if(length(predicted) != length(actual)){
    print("fn computeErrors: length wrong")
    return(NA)
  }
  confmat = matrix(0L, nrow=5, ncol=5)
  for(i in 1:length(predicted)){
    x = predicted[i]
    y = actual[i]
    confmat[x,y] = confmat[x,y] + 1
  }
  return(confmat)
}

computeClassErrors = function(mat,verbose=FALSE){
  c1Error = 1-mat[1,1]/sum(mat[,1])
  c2Error = 1-mat[2,2]/sum(mat[,2])
  c3Error = 1-mat[3,3]/sum(mat[,3])
  c4Error = 1-mat[4,4]/sum(mat[,4])
  c5Error = 1-mat[5,5]/sum(mat[,5])
  overallError = 1-(mat[1,1]+mat[2,2]+mat[3,3]+mat[4,4]+mat[5,5])/sum(mat)
  errors = list(c1Error,c2Error,c3Error,c4Error,c5Error,overallError)
  if(verbose){
    print(paste("Class 1 error: ",errors[[1]]))
    print(paste("Class 2 error: ",errors[[2]]))
    print(paste("Class 3 error: ",errors[[3]]))
    print(paste("Class 4 error: ",errors[[4]]))
    print(paste("Class 5 error: ",errors[[5]]))
    print(paste("Class average error: ",errors[[6]]))
  }
  return(errors)
}