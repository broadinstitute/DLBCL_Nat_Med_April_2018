pearson_custom = function(x,y,bins){
  if(length(x) != length(y)){
    print("Different lengths x y")
    return()
  }
  weights = c()
  for(i in 1:length(bins)){
    if(nrow(bins[[i]]) == 0){
      next
    }
    numCorrect = sum(bins[[i]]$correctness)
    trials = nrow(bins[[i]])
    lowerConf = prop.test(c(numCorrect),trials,conf.level=0.84, correct=FALSE)[[6]][1]
    upperConf = prop.test(c(numCorrect),trials,conf.level=0.84, correct=FALSE)[[6]][2]
    errorRange = upperConf-lowerConf
    weights = c(weights,1/(errorRange**2))
  }
  
  numerator = 0
  denominatorT1 = 0
  denominatorT2 = 0
  mx = mean(x)
  my = mean(y)
  for(i in 1:length(x)){
    term1 = (x[[i]]-mx)
    term2 = y[[i]]-my
    toAdd = weights[i]*term1*term2
    numerator = numerator+toAdd
    
    term1 = weights[i]*(term1**2)
    term2 = weights[i]*(term2**2)
    denominatorT1 = denominatorT1 + term1
    denominatorT2 = denominatorT2 + term2
  }
  
  denominatorT1 = sqrt(denominatorT1)
  denominatorT2 = sqrt(denominatorT2)
  denominator = denominatorT1*denominatorT2
  retval = numerator/denominator
  return(retval)
}

#prop.test(c(9),9,conf.level=0.95, correct=FALSE)
#prop.test(c(9),9,conf.level=0.95, correct=FALSE)[[6]][2]