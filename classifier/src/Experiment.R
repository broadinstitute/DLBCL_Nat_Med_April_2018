rm(list = ls())
source("src/GenerateFullDF2.R")
fullDF = ifelse(fullDF > 0, 1, 0)
#the probabilities of each gene occurring (value 1)
numGenes = ncol(fullDF)
marginalGeneProbabilities = vector("list", length=numGenes)
for(i in 1:numGenes){
  currProb = as.vector((table(fullDF[,i])/nrow(fullDF)))
  marginalGeneProbabilities[[i]] = currProb
}

sampleSize = 5000
#the mutational load distribution
mutationalLoadDistribution = c()
for(i in 1:sampleSize){
  numMuts = 0
  for(k in 1:numGenes){
    currProbs = marginalGeneProbabilities[[k]]
    call = sample(c(0,1), 1, prob = currProbs)
    numMuts = numMuts+call
  }
  mutationalLoadDistribution = c(mutationalLoadDistribution, numMuts)
  #mutationalLoadDistribution = c(mutationalLoadDistribution, 45)
}
distributionHistogram = hist(mutationalLoadDistribution, breaks=20)

constructedDF = data.frame()
#make the data frame
for(i in 1:sampleSize){
  if((i %% (sampleSize/10)) == 0){
    print(i)
  }
  #initialize new sample to all zeros
  currSample = matrix(0, nrow=1, ncol=numGenes)
  currSample = data.frame(currSample)
  colnames(currSample) = 1:numGenes
  
  #get the number of mutations for the sample from the mutationalLoadDistribution list
  numMutations = mutationalLoadDistribution[[i]]
  #current mutations start at 0
  currMutations = 0
  
  while(currMutations < numMutations){
    randomIndex = sample(1:numGenes,1)
    currProbs = marginalGeneProbabilities[[randomIndex]]
    currentCall = sample(c(0,1), 1, prob = currProbs)
    if(currentCall == 1){
      if(currSample[,randomIndex] == 0){
        currMutations = currMutations+1
        currSample[,randomIndex] = currentCall
      }
    }
  }
  constructedDF = rbind(constructedDF, currSample)
}

PValsAll = vector("list", length = ncol(constructedDF))
for(i in 1:(ncol(constructedDF))){
  print(i)
  pvals = c()
  for(j in 1:ncol(constructedDF)){
    if(i == j){
      next
    } else {
      p = fisher.test(x=factor(constructedDF[,i]), y=factor(constructedDF[,j]))[[1]]
      pvals = c(pvals, p)
    }
  }
  PValsAll[[i]] = pvals
}

for(i in 1:length(PValsAll)){
  feature = colnames(constructedDF)[[i]]
  title = paste("Pvals for ",feature,sep="")
  df = data.frame(PValsAll[[i]])
  colnames(df)[1] = "Pvals"
  p = ggplot(df, aes(x=Pvals)) +
    geom_histogram(fill="black", alpha=1, position="identity", bins=50) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5),
          title = element_text(size=20))+
    theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="grey20",size=15,angle=0,vjust=0,face="plain"),  
          axis.title.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
  
  fn = paste("Plots/PValuesNullModel/Experiment_PValues_",feature,".jpeg",sep="")
  jpeg(fn, width = 1080, height = 720)
  print(p)
  dev.off()
}

