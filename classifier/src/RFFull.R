
rm(list = ls())
knitr::opts_chunk$set(echo = TRUE)
library(randomForest)
library(reprtree)

source("src/GenerateFullDF.R")
source("src/GenerateTestTrainSets.R")
source("src/ComputeErrors.R")

orderedLabels = labels[rownames(fullDF),]
fullDF = cbind(factor(orderedLabels$cluster),fullDF)
colnames(fullDF)[1] = "cluster"


trainTable = fullDF[rownames(fullDF) %in% trainingSet,]
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
validationSet = rbind(v1,v2,v3,v4,v5)
trainTable = trainTable[!(rownames(trainTable) %in% rownames(validationSet)),]

rf = randomForest(cluster ~ ., trainTable, keep.forest=TRUE)
predicted = predict(rf, validationSet, predict.all = TRUE)
print("Predicted on X, Actual on Y")
confmat = computeErrors(predicted$aggregate,validationSet$cluster)
errors = computeClassErrors(confmat,verbose = TRUE)
confmat

reprtree:::plot.getTree(rf, k=1)
#dummy r chunk
