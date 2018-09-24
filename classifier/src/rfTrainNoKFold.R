rm(list = ls())
knitr::opts_chunk$set(echo = TRUE)
library(randomForest)
library(reprtree)
library(reshape2)
library(ggplot2)
library(caret)


source("src/GeneratedReducedDF.R")
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

#set.seed(111)
set.seed(1234567)
baselineProbTable = csBaselines(fullDF, connectivityProbs1.0, labels, summedW)
c1Subset = baselineProbTable[baselineProbTable$cluster == 1,]
c2Subset = baselineProbTable[baselineProbTable$cluster == 2,]
c3Subset = baselineProbTable[baselineProbTable$cluster == 3,]
c4Subset = baselineProbTable[baselineProbTable$cluster == 4,]
c5Subset = baselineProbTable[baselineProbTable$cluster == 5,]
#compute max vals for histogram
maxvals = c()
for(i in 1:nrow(baselineProbTable)){
  val = max(baselineProbTable[i,2:6])
  maxvals = c(maxvals,val)
}
baselineProbTable = cbind(baselineProbTable, maxvals)
print(paste("Number > 90%:",sum(baselineProbTable$maxvals > .90)))


orderedLabels = labels[rownames(reducedDF),]
reducedDF = cbind(factor(orderedLabels$cluster),reducedDF)
colnames(reducedDF)[1] = "cluster"

validationSize = .25
trainTable = reducedDF[rownames(reducedDF) %in% trainingSet,]
#trainTable = trainTable[sample(nrow(trainTable), nrow(trainTable)),]
subset1 = trainTable[trainTable$cluster == 1,]
subset2 = trainTable[trainTable$cluster == 2,]
subset3 = trainTable[trainTable$cluster == 3,]
subset4 = trainTable[trainTable$cluster == 4,]
subset5 = trainTable[trainTable$cluster == 5,]
v1 = subset1[sample(nrow(subset1),nrow(subset1)*validationSize),]
v2 = subset2[sample(nrow(subset2),nrow(subset2)*validationSize),]
v3 = subset3[sample(nrow(subset3),nrow(subset3)*validationSize),]
v4 = subset4[sample(nrow(subset4),nrow(subset4)*validationSize),]
v5 = subset5[sample(nrow(subset5),nrow(subset5)*validationSize),]
validationTable = rbind(v1,v2,v3,v4,v5)
trainTable = trainTable[!(rownames(trainTable) %in% rownames(validationTable)),]


#generate clusters from probability distribution
n = 1
nforests = 100
allClusterAssignments = vector("list", length = nforests)
for(c in 1:nforests){
  expandedClusters <- vector("list", length = nrow(trainTable)*n)
  for(i in 1:nrow(trainTable)){
    currSample = rownames(trainTable)[i]
    currProbs = baselineProbTable[rownames(baselineProbTable) == currSample, 2:6]
    clusters = sample(c(1,2,3,4,5), n, replace = TRUE, prob=currProbs)
    expandedClusters[[i]] = clusters
  }
  expandedClusters = unlist(expandedClusters)
  allClusterAssignments[[c]] = expandedClusters
}


#train nforests
rfList = vector("list", length = nforests)

for(i in 1:nforests){
  if((i %% floor(nforests/10)) == 0){
    print(i)
  }
  trainTable$cluster = allClusterAssignments[[i]]
  rf = randomForest(factor(cluster) ~ ., trainTable, keep.forest = TRUE, norm.votes = FALSE)
  rfList[[i]] = rf
}


allPredictions = vector("list", length = nforests)
prevnames = NA
for(i in 1:nforests){
  currRF = rfList[[i]]
  predictions = predict(currRF, validationTable, predict.all = TRUE)
  currRows = names(predictions$aggregate)
  if(i > 1){
    if(!identical(currRows,prevnames)){
      print("Order of predicted different")
    }
  }
  preds = as.integer(as.character(predictions$aggregate))
  #allPredictions[[i]] = predictions$aggregate
  names(preds) = currRows
  allPredictions[[i]] = preds
  prevnames = currRows
}
allPredictions = as.data.frame(allPredictions)
colnames(allPredictions) = 1:nforests

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
  assignedCluster = which(max(probVec) == probVec)[1]
  trueCluster = labels[rownames(labels) == rownames(allPredictions)[i],"cluster"]
  consensusPredictions = rbind(consensusPredictions, c(probVec, maxVal, assignedCluster, trueCluster))
}
rownames(consensusPredictions) = rownames(allPredictions)
colnames(consensusPredictions) = c("P1","P2","P3","P4","P5","maxval","Predicted","TrueCluster")



validationProbFrame = baselineProbTable[rownames(baselineProbTable) %in% rownames(validationTable),]
validationProbFrame = validationProbFrame[rownames(consensusPredictions),]

plotDF = data.frame()
correctness = consensusPredictions$Predicted == consensusPredictions$TrueCluster
plotDF = cbind(consensusPredictions$maxval, validationProbFrame$maxvals, correctness)
plotDF = data.frame(plotDF)
plotDF$correctness = as.logical(plotDF$correctness)
colnames(plotDF) = c("maxPred","maxActual", "correctness")
p = ggplot(plotDF, aes(x=maxActual, y=maxPred, colour = correctness)) + geom_point() +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = .1))


print(p)


#high confidence agreement
highConfSubset = consensusPredictions[consensusPredictions$maxval >= .90,]
confMat = confusionMatrix(factor(highConfSubset$Predicted), factor(highConfSubset$TrueCluster))

