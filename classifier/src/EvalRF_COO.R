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



#set.seed(111)
set.seed(1234)
baselineProbTable = csBaselines(fullDF, connectivityProbs1.0, labels, summedW)
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
print(paste("Number > 90%:",sum(baselineProbTable$maxvals > .90)))
write.table(baselineProbTable, "/DataTables/baselineProbTable.txt", sep='\t', row.names = TRUE, col.names = TRUE)



orderedLabels = labels[rownames(reducedDFCOO),]
reducedDFCOO = cbind(factor(orderedLabels$cluster),reducedDFCOO)
colnames(reducedDFCOO)[1] = "cluster"



trainTable = reducedDFCOO[rownames(reducedDFCOO) %in% trainingSet,]
subset1 = trainTable[trainTable$cluster == 1,]
subset2 = trainTable[trainTable$cluster == 2,]
subset3 = trainTable[trainTable$cluster == 3,]
subset4 = trainTable[trainTable$cluster == 4,]
subset5 = trainTable[trainTable$cluster == 5,]
subset1 = subset1[sample(nrow(subset1)),]
subset2 = subset2[sample(nrow(subset2)),]
subset3 = subset3[sample(nrow(subset3)),]
subset4 = subset4[sample(nrow(subset4)),]
subset5 = subset5[sample(nrow(subset5)),]

v1 = c()
v2 = c()
v3 = c()
v4 = c()
v5 = c()

for(i in 1:nrow(subset1)){
  sample = rownames(subset1)[i]
  if((i %% 5) == 1){
    v1 = c(v1,sample)
  } else if((i %% 5) == 2){
    v2 = c(v2,sample)
  } else if((i %% 5) == 3){
    v3 = c(v3,sample)
  } else if((i %% 5) == 4){
    v4 = c(v4,sample)
  } else if((i %% 5) == 0){
    v5 = c(v5,sample)
  }
}
for(i in 1:nrow(subset2)){
  sample = rownames(subset2)[i]
  if((i %% 5) == 1){
    v1 = c(v1,sample)
  } else if((i %% 5) == 2){
    v2 = c(v2,sample)
  } else if((i %% 5) == 3){
    v3 = c(v3,sample)
  } else if((i %% 5) == 4){
    v4 = c(v4,sample)
  } else if((i %% 5) == 0){
    v5 = c(v5,sample)
  }
}
for(i in 1:nrow(subset3)){
  sample = rownames(subset3)[i]
  if((i %% 5) == 1){
    v1 = c(v1,sample)
  } else if((i %% 5) == 2){
    v2 = c(v2,sample)
  } else if((i %% 5) == 3){
    v3 = c(v3,sample)
  } else if((i %% 5) == 4){
    v4 = c(v4,sample)
  } else if((i %% 5) == 0){
    v5 = c(v5,sample)
  }
}
for(i in 1:nrow(subset4)){
  sample = rownames(subset4)[i]
  if((i %% 5) == 1){
    v1 = c(v1,sample)
  } else if((i %% 5) == 2){
    v2 = c(v2,sample)
  } else if((i %% 5) == 3){
    v3 = c(v3,sample)
  } else if((i %% 5) == 4){
    v4 = c(v4,sample)
  } else if((i %% 5) == 0){
    v5 = c(v5,sample)
  }
}
for(i in 1:nrow(subset5)){
  sample = rownames(subset5)[i]
  if((i %% 5) == 1){
    v1 = c(v1,sample)
  } else if((i %% 5) == 2){
    v2 = c(v2,sample)
  } else if((i %% 5) == 3){
    v3 = c(v3,sample)
  } else if((i %% 5) == 4){
    v4 = c(v4,sample)
  } else if((i %% 5) == 0){
    v5 = c(v5,sample)
  }
}

train1 = trainTable[!(rownames(trainTable) %in% v1),]
train2 = trainTable[!(rownames(trainTable) %in% v2),]
train3 = trainTable[!(rownames(trainTable) %in% v3),]
train4 = trainTable[!(rownames(trainTable) %in% v4),]
train5 = trainTable[!(rownames(trainTable) %in% v5),]

v1 = trainTable[rownames(trainTable) %in% v1,]
v2 = trainTable[rownames(trainTable) %in% v2,]
v3 = trainTable[rownames(trainTable) %in% v3,]
v4 = trainTable[rownames(trainTable) %in% v4,]
v5 = trainTable[rownames(trainTable) %in% v5,]

trainSets = list(train1,train2,train3,train4,train5)
validationSets = list(v1,v2,v3,v4,v5)



xs = vector("list", length = 5)
ys = vector("list", length = 5)
fullPredictions = vector("list", length=5)
fullActuals = vector("list", length = 5)
realProbs = vector("list", length = 5)
predictedProbs = vector("list", length = 5)
#massive loop incoming. It's just k fold, doing the same code in the other file.
for(k in 1:length(trainSets)){
  print(paste("fold",k))
  n = 1
  nforests = 100
  allClusterAssignments = vector("list", length = nforests)
  for(c in 1:nforests){
    expandedClusters <- vector("list", length = nrow(trainSets[[k]])*n)
    for(i in 1:nrow(trainSets[[k]])){
      currSample = rownames(trainSets[[k]])[i]
      currProbs = baselineProbTable[rownames(baselineProbTable) == currSample, 1:5]
      clusters = sample(c(1,2,3,4,5), n, replace = TRUE, prob=currProbs)
      expandedClusters[[i]] = clusters
    }
    expandedClusters = unlist(expandedClusters)
    allClusterAssignments[[c]] = expandedClusters
  }
  
  rfList = vector("list", length = nforests)
  trainSets[[k]]$cluster = allClusterAssignments[[1]]
  for(i in 1:nforests){
    if((i %% floor(nforests/2)) == 0){
      print(i)
    }
    trainSets[[k]]$cluster = allClusterAssignments[[i]]
    rf = randomForest(factor(cluster) ~ ., trainSets[[k]], keep.forest = TRUE, norm.votes = FALSE)
    rfList[[i]] = rf
  }
  
  allPredictions = vector("list", length = nforests)
  prevnames = NA
  for(i in 1:nforests){
    currRF = rfList[[i]]
    predictions = predict(currRF, validationSets[[k]], predict.all = TRUE)
    currRows = names(predictions$aggregate)
    if(i > 1){
      if(!identical(currRows,prevnames)){
        print("Order of predicted different")
      }
    }
    preds = as.integer(as.character(predictions$aggregate))
    names(preds) = currRows
    allPredictions[[i]] = preds
    prevnames = currRows
  }
  allPredictions = as.data.frame(allPredictions)
  colnames(allPredictions) = 1:nforests
  validationProbFrame = baselineProbTable[rownames(baselineProbTable) %in% rownames(validationSets[[k]]),]
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
    assignedCluster = which(max(probVec) == probVec)[1]
    consensusPredictions = rbind(consensusPredictions, c(probVec, maxVal, assignedCluster))
    
  }
  rownames(consensusPredictions) = rownames(allPredictions)
  colnames(consensusPredictions) = c("P1","P2","P3","P4","P5","maxval","Predicted")
  
  consensusPredictions = consensusPredictions[order(consensusPredictions$maxval),]
  validationProbFrame = validationProbFrame[rownames(consensusPredictions),]
  
  truthsP = c()
  predsP = c()
  for(i in 1:nrow(validationProbFrame)){
    correctClus = as.numeric(as.character(validationProbFrame[i,"cluster"]))
    trueP = as.numeric(as.character(validationProbFrame[i, "maxvals"]))
    predP = as.numeric(as.character(consensusPredictions[i,correctClus]))
    truthsP = c(truthsP, trueP)
    predsP = c(predsP, predP)
  }
  
  realProbs[[k]] = truthsP
  predictedProbs[[k]] = predsP
  fullPredictions[[k]] = consensusPredictions$Predicted
  fullActuals[[k]] = validationProbFrame$cluster
  xs[[k]] = consensusPredictions$maxval
  ys[[k]] = validationProbFrame$maxvals
}



plotDF = cbind(unlist(xs), unlist(ys))
absDiffs = c()
for(i in 1:nrow(plotDF)){
  diff = abs(plotDF[i,1] - plotDF[i,2])
  absDiffs = c(absDiffs, diff)
}

plotDF = cbind(plotDF, absDiffs)
plotDF = cbind(plotDF, c(rownames(v1),rownames(v2),rownames(v3),rownames(v4),rownames(v5)))
plotDF = cbind(plotDF, unlist(fullPredictions), unlist(fullActuals))
plotDF = data.frame(plotDF)
colnames(plotDF) = c("maxPred", "maxActual", "absDiffs", "Sample", "Predicted", "Actual")
plotDF = cbind(plotDF, (plotDF$Predicted == plotDF$Actual))
colnames(plotDF)[ncol(plotDF)] = "Correct"
plotDF$maxPred = as.numeric(as.character(plotDF$maxPred))
plotDF$maxActual = as.numeric(as.character(plotDF$maxActual))
plotDF$absDiffs = as.numeric(as.character(plotDF$absDiffs))
#colnames(plotDF) = c("maxPred", "maxActual")




outliers = subset(plotDF, absDiffs > .35)
p = ggplot(plotDF, aes(x=maxActual, y=maxPred)) + geom_point() + geom_smooth(method=lm) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  #geom_text(data=outliers, 
  #(aes(x=maxActual, y=maxPred,label=Sample)), size=3) +
  geom_point(data=subset(plotDF, Correct == FALSE), aes(x=maxActual, y=maxPred), colour = "red", size=1)
jpeg("/Plots/outliersplot1Power2.jpeg", width = 1080, height = 720)
p
dev.off()
p
outliers = reducedDFCOO[outliers$Sample,]
write.table(outliers, "/WrittenFiles/outliersplot1Power2.txt", sep='\t', row.names = TRUE, col.names = TRUE)




#plot the points of Actual probability vs Predicted probability for correct class
plotDF = cbind(unlist(predictedProbs), unlist(realProbs))
absDiffs = c()
for(i in 1:nrow(plotDF)){
  diff = abs(plotDF[i,1] - plotDF[i,2])
  absDiffs = c(absDiffs, diff)
}

plotDF = cbind(plotDF, absDiffs)
plotDF = cbind(plotDF, c(rownames(v1),rownames(v2),rownames(v3),rownames(v4),rownames(v5)))
plotDF = cbind(plotDF, unlist(fullPredictions), unlist(fullActuals))
plotDF = data.frame(plotDF)
colnames(plotDF) = c("Pred", "TrueCluster", "absDiffs", "Sample", "Predicted", "Actual")
plotDF = cbind(plotDF, (plotDF$Predicted == plotDF$Actual))
colnames(plotDF)[ncol(plotDF)] = "Correct"
plotDF$Pred = as.numeric(as.character(plotDF$Pred))
plotDF$TrueCluster = as.numeric(as.character(plotDF$TrueCluster))
plotDF$absDiffs = as.numeric(as.character(plotDF$absDiffs))
#colnames(plotDF) = c("maxPred", "maxActual")



outliers = subset(plotDF, absDiffs > .40)
p = ggplot(plotDF, aes(x=TrueCluster, y=Pred)) + geom_point() + geom_smooth(method=lm) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  #geom_text(data=outliers, 
  #(aes(x=TrueCluster, y=Pred,label=Sample)), size=3) +
  geom_point(data=subset(plotDF, Correct == FALSE), aes(x=TrueCluster, y=Pred), colour = "red", size=1)
jpeg("/Plots/outliersplot2Power2.jpeg", width = 1080, height = 720)
p
dev.off()
p
outliers = reducedDFCOO[outliers$Sample,]
write.table(outliers, "/WrittenFiles/outliersplot2Power2.txt", sep='\t', row.names = TRUE, col.names = TRUE)

cor.test(plotDF$Pred, plotDF$TrueCluster, method = "pearson")

