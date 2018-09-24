
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



orderedLabels = labels[rownames(reducedDFCOO),]
reducedDFCOO = cbind(factor(orderedLabels$cluster),reducedDFCOO)
colnames(reducedDFCOO)[1] = "cluster"



set.seed(123456)
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
validationSet = rbind(v1,v2,v3,v4,v5)
trainTable = trainTable[!(rownames(trainTable) %in% rownames(validationSet)),]



size = 11
rfList <- vector("list", length = size)
allPredictions = vector("list", length = size)
allVotes = vector("list", length = size)
prev = NA

for(i in 1:size){
  set.seed(i*100)
  rf = randomForest(cluster ~ ., trainTable, keep.forest=TRUE, norm.votes=FALSE)
  rfList[[i]] = rf
  predictedClass = predict(rf, validationSet, predict.all = TRUE)$aggregate
  predictedVotes = predict(rf, validationSet, predict.all = TRUE, type="vote", norm.votes = FALSE)$aggregate
  
  allPredictions[[i]] = as.numeric(as.vector(predictedClass))
  allVotes[[i]] = predictedVotes
  prev = names(predictedClass)
}
predDF = as.data.frame(allPredictions)
rownamesPred = prev
predDF = as.matrix(predDF)
rownames(predDF) = rownamesPred
colnames(predDF) = 1:size



predictionMats <- vector("list", length = size)
for(i in 1:ncol(predDF)){
  errMat = computeErrors(predDF[,i],validationSet$cluster)
  predictionMats[[i]] = errMat
  classError = computeClassErrors(errMat,verbose=FALSE)
}



errors = ensembleErrors(predDF,validationSet$cluster,verb=TRUE)
probs = probabilityTable(predDF)


summedVotes = allVotes[[1]]
for(i in 2:length(rfList)){
  summedVotes = summedVotes+allVotes[[i]]
}
summedVotes = apply(summedVotes, 1, function(x) x/sum(x))
tmp = t(summedVotes)
rownames(tmp) = colnames(summedVotes)
colnames(tmp) = rownames(summedVotes)
summedVotes = tmp
maxvals = c()
for(i in 1:nrow(summedVotes)){
  maxvals = c(maxvals, max(summedVotes[i,]))
}
summedVotes = cbind(summedVotes, maxvals)
summedVotes = cbind(summedVotes, errors[[3]])
colnames(summedVotes)[ncol(summedVotes)] = "predicted"
summedVotes = cbind(summedVotes, validationSet$cluster)
colnames(summedVotes)[ncol(summedVotes)] = "actual"
summedVotes = data.frame(summedVotes)
summedVotes = summedVotes[order(summedVotes$maxvals, decreasing = TRUE),]
top70th = head(summedVotes, nrow(summedVotes)*.7)
bottom30th = summedVotes[!(rownames(summedVotes) %in%rownames(top70th)),]



table(top70th$actual == top70th$predicted)
table(bottom30th$actual == bottom30th$predicted)
hconfMat = computeErrors(top70th$predicted, top70th$actual)
lconfMat = computeErrors(bottom30th$predicted, bottom30th$actual)
hConfClass = computeClassErrors(hconfMat, verbose = TRUE)
lconfClass = computeClassErrors(lconfMat, verbose = TRUE)
hconfMat
lconfMat



baselineProbTable = csBaselines(fullDF, connectivityProbs1.0, labels, summedW) 
# n1s = baselineProbTable$n1s
# n1s = data.frame(n1s)
# row.names(n1s) = rownames(baselineProbTable)
baselineProbTable = baselineProbTable[,-1]
maxval = c()
for(i in 1:nrow(baselineProbTable)){
  val = max(baselineProbTable[i,1:5])
  maxval = c(maxval,val)
}
baselineProbTable = cbind(baselineProbTable, maxval)
orderedBPT = baselineProbTable[order(baselineProbTable$maxval, decreasing = TRUE),]



trueRanks = 1:nrow(top70th)
matchedRanks = match(rownames(top70th), rownames(orderedBPT[rownames(orderedBPT) %in% rownames(top70th),]))
cor.test(trueRanks, matchedRanks, method="spearman")



#randomforest heatmap. Samples sorted by predicted cluster and confidence horizontally, features sorted
#by W matrix amplitude summed importance then binned by appropriate cluster by max value
importanceDF = NA
for(i in 1:length(rfList)){
  currImp = rfList[[i]]$importance
  currImp = currImp[order(rownames(currImp)),]
  if(i > 1){
    importanceDF = importanceDF + currImp
  } else{
    importanceDF = currImp
  }
}

c1Features = c("SV.BCL6","B2M.CD70.PD1.FAS.com","BLC10.NOTCH2.A20.SPEN", "SV.TP63", "MYD88.OTHER")
c2Features = c("Chrom.Mod.Enzymes","BCL2.comp","Pi3K.Mod")
c3Features = c("Sum.CNAs.C2","TP53.biallelic","GENOME_DOUBLING","X2P16.1.AMP","X21Q.AMP","X9p21.3.DEL")
c4Features = c("CD79B.MYD88.TBL1XR1.BLIMP.ETV6.ZC3H12A","sumCNAs.C5","COO..GCB.1.unclass.2.ABC.3.na.","X18Q.AMP","MYD88.L265")
c5Features = c("SumC4.mut","NFKBalt.RASalt.JAK.STAT.alt.RHOA.SGK1.KLH6.cd58.cd83","Hist.comp")
featureOrder = c(c1Features,c2Features,c3Features,c4Features,c5Features)



#make heatmap
summedVotes = summedVotes[order(summedVotes$maxvals),]
c1Sub = summedVotes[summedVotes$predicted == 1,]
c2Sub = summedVotes[summedVotes$predicted == 2,]
c3Sub = summedVotes[summedVotes$predicted == 3,]
c4Sub = summedVotes[summedVotes$predicted == 4,]
c5Sub = summedVotes[summedVotes$predicted == 5,]
sampleOrder = c(rownames(c1Sub),rownames(c2Sub),rownames(c3Sub),rownames(c4Sub),rownames(c5Sub))
orderedSummed = summedVotes[sampleOrder,]



prntImage = TRUE
standardizedReduced = standardize(reducedDFCOO)
heatmapMatrix = standardizedReduced
heatmapMatrix = heatmapMatrix[,-1]
heatmapMatrix = heatmapMatrix[rownames(heatmapMatrix) %in% rownames(orderedSummed),rev(featureOrder)]
heatmapMatrix = heatmapMatrix[rownames(orderedSummed),]
#predicted = orderedSummed$predicted
#heatmapMatrix = cbind(predicted,heatmapMatrix)
agreeVec = orderedSummed$predicted == orderedSummed$actual
agreeVec = rep(agreeVec,ncol(heatmapMatrix))
heatmapMatrix = as.matrix(heatmapMatrix)
melted = melt(heatmapMatrix)
melted = cbind(melted,agreeVec)
ggplot(data = melted, aes(x=Var1, y=Var2)) + 
  geom_tile(aes(fill=agreeVec, alpha = value) , colour = "white") +
  scale_y_discrete(expand=c(0,0))+
  #define new breaks on x-axis
  scale_x_discrete(expand=c(0,0), breaks=c(rownames(c1Sub)[[1]]
                                           , rownames(c2Sub)[[1]]
                                           , rownames(c3Sub)[[1]]
                                           , rownames(c4Sub)[[1]]
                                           , rownames(c5Sub)[[1]])
                   , labels = c("C1","C2","C3","C4","C5")) +
  theme(axis.text.y = element_text(face = rev(c('bold', 'plain', 'plain', 'plain', 'plain',
                                                'bold', 'plain','plain',
                                                'bold','plain','plain', 'plain','plain','plain',
                                                'bold','plain','plain','plain','plain',
                                                'bold','plain','plain'))))

jpeg("Plots/heatmapProbRedCOO.jpeg", width = 1080, height = 720)
ggplot(data = melted, aes(x=Var1, y=Var2)) + 
  geom_tile(aes(fill=agreeVec, alpha = value) , colour = "white") +
  scale_y_discrete(expand=c(0,0))+
  #define new breaks on x-axis
  scale_x_discrete(expand=c(0,0), breaks=c(rownames(c1Sub)[[1]]
                                           , rownames(c2Sub)[[1]]
                                           , rownames(c3Sub)[[1]]
                                           , rownames(c4Sub)[[1]]
                                           , rownames(c5Sub)[[1]])
                   , labels = c("C1","C2","C3","C4","C5")) +
  theme(axis.text.y = element_text(face = rev(c('bold', 'plain', 'plain', 'plain', 'plain',
                                                'bold', 'plain','plain',
                                                'bold','plain','plain', 'plain','plain','plain',
                                                'bold','plain','plain','plain','plain',
                                                'bold','plain','plain'))))
dev.off()




#validation set heatmap
subsetDF = baselineProbTable[rownames(baselineProbTable) %in% rownames(validationSet),]
c1Sub = subsetDF[subsetDF$cluster == 1,]
c2Sub = subsetDF[subsetDF$cluster == 2,]
c3Sub = subsetDF[subsetDF$cluster == 3,]
c4Sub = subsetDF[subsetDF$cluster == 4,]
c5Sub = subsetDF[subsetDF$cluster == 5,]
x1 = nrow(c1Sub)+1
x2 = x1+nrow(c2Sub)
x3 = x2+nrow(c3Sub)
x4 = x3+nrow(c4Sub)
y1 = length(featureOrder)-length(c1Features)
y2 = y1-length(c2Features)
y3 = y2-length(c3Features)
y4 = y3-length(c4Features)
sampleOrder = c(rownames(c1Sub),rownames(c2Sub),rownames(c3Sub),rownames(c4Sub),rownames(c5Sub))
id1 = "DLBCL_RICOVER_473"
prntImage = FALSE
standardizedReduced = standardize(reducedDFCOO)
heatmapMatrix = standardizedReduced
heatmapMatrix = heatmapMatrix[,-1]
heatmapMatrix = heatmapMatrix[rownames(heatmapMatrix) %in% rownames(subsetDF),rev(featureOrder)]
heatmapMatrix = heatmapMatrix[sampleOrder,]
heatmapMatrix = as.matrix(heatmapMatrix)
melted = melt(heatmapMatrix)
ggplot(data = melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") +
  scale_y_discrete(expand=c(0,0))+
  #define new breaks on x-axis
  scale_x_discrete(expand=c(0,0), breaks=c(rownames(c1Sub)[[1]]
                                           , rownames(c2Sub)[[1]]
                                           , rownames(c3Sub)[[1]]
                                           , rownames(c4Sub)[[1]]
                                           , rownames(c5Sub)[[1]]
                                           , id1)
                   , labels = c("C1","C2","C3","C4","C5", "S")) +
  theme(axis.text.y = element_text(face = rev(c('bold', 'plain', 'plain', 'plain', 'plain',
                                                'bold', 'plain','plain',
                                                'bold','plain','plain', 'plain','plain','plain',
                                                'bold','plain','plain','plain','plain',
                                                'bold','plain','plain')))) +
  geom_vline(xintercept = c(x1,x2,x3,x4), linetype = "dashed") +
  geom_hline(yintercept = c(y1,y2,y3,y4), linetype = "dashed")

jpeg("Plots/heatmapTrueOrderingCOO.jpeg", width = 1920, height = 1080)
ggplot(data = melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") +
  scale_y_discrete(expand=c(0,0))+
  #define new breaks on x-axis
  scale_x_discrete(expand=c(0,0), breaks=c(rownames(c1Sub)[[1]]
                                           , rownames(c2Sub)[[1]]
                                           , rownames(c3Sub)[[1]]
                                           , rownames(c4Sub)[[1]]
                                           , rownames(c5Sub)[[1]]
                                           , id1)
                   , labels = c("C1","C2","C3","C4","C5","S")) +
  theme(axis.text.y = element_text(face = rev(c('bold', 'plain', 'plain', 'plain', 'plain',
                                                'bold', 'plain','plain',
                                                'bold','plain','plain', 'plain','plain','plain',
                                                'bold','plain','plain','plain','plain',
                                                'bold','plain','plain')), size = 16)) +
  geom_vline(xintercept = c(x1,x2,x3,x4), linetype = "dashed") +
  geom_hline(yintercept = c(y1,y2,y3,y4), linetype = "dashed")
dev.off()





subsetDF = baselineProbTable[rownames(baselineProbTable) %in% rownames(reducedDFCOO),]
c1Sub = subsetDF[subsetDF$cluster == 1,]
c2Sub = subsetDF[subsetDF$cluster == 2,]
c3Sub = subsetDF[subsetDF$cluster == 3,]
c4Sub = subsetDF[subsetDF$cluster == 4,]
c5Sub = subsetDF[subsetDF$cluster == 5,]
x1 = nrow(c1Sub)+1
x2 = x1+nrow(c2Sub)
x3 = x2+nrow(c3Sub)
x4 = x3+nrow(c4Sub)
y1 = length(featureOrder)-length(c1Features)
y2 = y1-length(c2Features)
y3 = y2-length(c3Features)
y4 = y3-length(c4Features)
sampleOrder = c(rownames(c1Sub),rownames(c2Sub),rownames(c3Sub),rownames(c4Sub),rownames(c5Sub))
prntImage = FALSE
standardizedReduced = standardize(reducedDFCOO)
heatmapMatrix = standardizedReduced
heatmapMatrix = heatmapMatrix[,-1]
heatmapMatrix = heatmapMatrix[rownames(heatmapMatrix) %in% rownames(subsetDF),rev(featureOrder)]
heatmapMatrix = heatmapMatrix[sampleOrder,]
heatmapMatrix = as.matrix(heatmapMatrix)
melted = melt(heatmapMatrix)
ggplot(data = melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") +
  scale_y_discrete(expand=c(0,0))+
  #define new breaks on x-axis
  scale_x_discrete(expand=c(0,0), breaks=c(rownames(c1Sub)[[1]]
                                           , rownames(c2Sub)[[1]]
                                           , rownames(c3Sub)[[1]]
                                           , rownames(c4Sub)[[1]]
                                           , rownames(c5Sub)[[1]])
                   , labels = c("C1","C2","C3","C4","C5")) +
  theme(axis.text.y = element_text(face = rev(c('bold', 'plain', 'plain', 'plain', 'plain',
                                                'bold', 'plain','plain',
                                                'bold','plain','plain', 'plain','plain','plain',
                                                'bold','plain','plain','plain','plain',
                                                'bold','plain','plain')))) +
  geom_vline(xintercept = c(x1,x2,x3,x4), linetype = "dashed") +
  geom_hline(yintercept = c(y1,y2,y3,y4), linetype = "dashed")
jpeg("Plots/heatmapAllReducedAlpha.0.9.Power1.2.constantScalar.jpeg", width = 1920, height = 1080)
ggplot(data = melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") +
  scale_y_discrete(expand=c(0,0))+
  #define new breaks on x-axis
  scale_x_discrete(expand=c(0,0), breaks=c(rownames(c1Sub)[[1]]
                                           , rownames(c2Sub)[[1]]
                                           , rownames(c3Sub)[[1]]
                                           , rownames(c4Sub)[[1]]
                                           , rownames(c5Sub)[[1]])
                   , labels = c("C1","C2","C3","C4","C5")) +
  theme(axis.text.y = element_text(face = rev(c('bold', 'plain', 'plain', 'plain', 'plain',
                                                'bold', 'plain','plain',
                                                'bold','plain','plain', 'plain','plain','plain',
                                                'bold','plain','plain','plain','plain',
                                                'bold','plain','plain')), size = 16)) +
  geom_vline(xintercept = c(x1,x2,x3,x4), linetype = "dashed") +
  geom_hline(yintercept = c(y1,y2,y3,y4), linetype = "dashed")
dev.off()






#dummy chunk for run all above button
reprtree:::plot.getTree(rfList[[1]], k=12)


