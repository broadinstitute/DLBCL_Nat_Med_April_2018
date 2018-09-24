
rm(list = ls())
knitr::opts_chunk$set(echo = TRUE)
library(randomForest)
library(reprtree)
library(reshape2)
library(ggplot2)


source("src/GeneratedReducedCOO.R")
source("src/GenerateTestTrainSets.R")
source("src/ComputeErrors.R")
source("src/ComputeEnsembleRF.R")
source("src/ProbabilityTable.R")
source("src/addBaselines.R")
source("src/addCSBaselinesV2.R")
source("src/GenerateProbabilityTable.R")
source("src/SumRankWmats.R")
source("src/StandardizeReducedFeatures.R")
source("src/clusOrderingsFull.R")



baselineProbTable = csBaselines(fullDF, connectivityProbs1.0, labels, summedW)
n1s = baselineProbTable$n1s
baselineProbTable[,-1]

c1Subset = baselineProbTable[baselineProbTable$cluster == 1,]
c2Subset = baselineProbTable[baselineProbTable$cluster == 2,]
c3Subset = baselineProbTable[baselineProbTable$cluster == 3,]
c4Subset = baselineProbTable[baselineProbTable$cluster == 4,]
c5Subset = baselineProbTable[baselineProbTable$cluster == 5,]
#compute max vals for histogram
maxvals = c()
for(i in 1:nrow(baselineProbTable)){
  val = max(baselineProbTable[i,1:5])
  maxvals = c(maxvals,val)
}
baselineProbTable = cbind(baselineProbTable, maxvals)



orderedLabels = labels[rownames(reducedDFCOO),]
reducedDFCOO = cbind(factor(orderedLabels$cluster),reducedDFCOO)
colnames(reducedDFCOO)[1] = "cluster"



set.seed(444)
trainTable = reducedDFCOO[rownames(reducedDFCOO) %in% trainingSet,]
subset1 = trainTable[trainTable$cluster == 1,]
subset2 = trainTable[trainTable$cluster == 2,]
subset3 = trainTable[trainTable$cluster == 3,]
subset4 = trainTable[trainTable$cluster == 4,]
subset5 = trainTable[trainTable$cluster == 5,]
v1 = subset1[sample(nrow(subset1),nrow(subset1)*.2),]
v2 = subset2[sample(nrow(subset2),nrow(subset2)*.2),]
v3 = subset3[sample(nrow(subset3),nrow(subset3)*.2),]
v4 = subset4[sample(nrow(subset4),nrow(subset4)*.2),]
v5 = subset5[sample(nrow(subset5),nrow(subset5)*.2),]
validationTable = rbind(v1,v2,v3,v4,v5)
trainTable = trainTable[!(rownames(trainTable) %in% rownames(validationTable)),]



#generate clusters from probability distribution
n = 100
nforests = 50
allClusterAssignments = vector("list", length = nforests)
for(c in 1:1){
  expandedClusters <- vector("list", length = nrow(trainTable)*n)
  for(i in 1:nrow(trainTable)){
    currSample = rownames(trainTable)[i]
    currProbs = baselineProbTable[rownames(baselineProbTable) == currSample, 1:5]
    clusters = sample(c(1,2,3,4,5), n, replace = TRUE, prob=currProbs)
    expandedClusters[[i]] = clusters
  }
  expandedClusters = unlist(expandedClusters)
  allClusterAssignments[[c]] = expandedClusters
}



trainTable = trainTable[rep(seq_len(nrow(trainTable)), each=n),]



rfList = vector("list", length = nforests)
trainTable$cluster = allClusterAssignments[[1]]
for(i in 1:nforests){
  print(i)
  #trainTable$cluster = allClusterAssignments[[i]]
  rf = randomForest(factor(cluster) ~ ., trainTable, keep.forest = TRUE, norm.votes = FALSE)
  rfList[[i]] = rf
}



allPredictions = vector("list", length = nforests)
for(i in 1:nforests){
  currRF = rfList[[i]]
  predictions = predict(currRF, validationTable, predict.all = TRUE)
  currRows = names(predictions$aggregate)
  preds = as.integer(as.character(predictions$aggregate))
  #allPredictions[[i]] = predictions$aggregate
  names(preds) = currRows
  allPredictions[[i]] = preds
}
allPredictions = as.data.frame(allPredictions)
colnames(allPredictions) = 1:nforests

validationProbFrame = baselineProbTable[rownames(baselineProbTable) %in% rownames(validationTable),]
validationProbFrame = validationProbFrame[rownames(allPredictions),]



consensusPredictions = data.frame()
for(i in 1:nrow(allPredictions)){
  probVec = vector("integer",5)
  for(j in 1:length(allPredictions[i,])){
    val = allPredictions[i,j]
    if(val == "1"){
      probVec[1] = probVec[1] + 1
    } else if (val == "2"){
      probVec[2] = probVec[2] + 1
    } else if (val == "3"){
      probVec[3] = probVec[3] + 1
    } else if (val == "4"){
      probVec[4] = probVec[4] + 1
    } else if (val == "5"){
      probVec[5] = probVec[5] + 1
    }
  }
  probVec = probVec/sum(probVec)
  maxVal = max(probVec)
  assignedCluster = which(max(probVec) == probVec)
  consensusPredictions = rbind(consensusPredictions, c(probVec, maxVal, assignedCluster))
}
rownames(consensusPredictions) = rownames(allPredictions)
colnames(consensusPredictions) = c("P1","P2","P3","P4","P5","maxval","Predicted")



#wrong samples
wrongSamples = rownames(validationProbFrame)[which(!(consensusPredictions == validationProbFrame$cluster))]
validationProbFrame[rownames(validationProbFrame) %in% wrongSamples,]



#top 70%
sortedProbVec = consensusPredictions[order(consensusPredictions$maxval, decreasing = TRUE),]
top70th = sortedProbVec[1:(nrow(sortedProbVec)*.70),]
bot70th = sortedProbVec[!(rownames(sortedProbVec) %in% rownames(top70th)),]



consensusPredictions = consensusPredictions[order(consensusPredictions$maxval),]
validationProbFrame = validationProbFrame[rownames(consensusPredictions),]



plotDF = cbind(consensusPredictions$maxval, validationProbFrame$maxvals)
plotDF = data.frame(plotDF)
colnames(plotDF) = c("maxPred","maxActual")
ggplot(plotDF, aes(x=maxPred, y=maxActual)) + geom_point() + geom_smooth(method=lm) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1))



