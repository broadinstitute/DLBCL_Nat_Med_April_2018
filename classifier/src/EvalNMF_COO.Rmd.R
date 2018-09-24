---
  title: "EvalNMF_COO"
author: "Tim"
date: "9/18/2018"
output: html_document
---
  
  {r setup, include=FALSE}
rm(list = ls())
knitr::opts_chunk$set(echo = TRUE)
library(reshape2)



set.seed(123)

source("src/GenerateLabels.R")
source("src/GeneratedReducedCOO.R")
source("src/GenerateTestTrainSets.R")
source("src/reorientG1.R")
source("src/PearsonCustom.R")
source("src/StandardizeReducedFeatures.R")



orderedLabels = labels[rownames(reducedDFCOO),]
reducedDFCOO = cbind(factor(orderedLabels$cluster),reducedDFCOO)
colnames(reducedDFCOO)[1] = "cluster"
validationSize = 4
validationSets = vector("list", validationSize)
trainTable = reducedDFCOO[rownames(reducedDFCOO) %in% trainingSet,]
trainTable = trainTable[sample(nrow(trainTable), nrow(trainTable)),]
subset1 = trainTable[trainTable$cluster == 1,]
subset2 = trainTable[trainTable$cluster == 2,]
subset3 = trainTable[trainTable$cluster == 3,]
subset4 = trainTable[trainTable$cluster == 4,]
subset5 = trainTable[trainTable$cluster == 5,]
validationSets = list()
for(i in 1:validationSize){
  validationSets[[i]] = c(NA)
}
for(i in 0:(nrow(subset1)-1)){
  idx = (i %% validationSize)+1
  validationSets[[idx]] = c(validationSets[[idx]], rownames(subset1)[i+1])
}
for(i in 0:(nrow(subset2)-1)){
  idx = (i %% validationSize)+1
  validationSets[[idx]] = c(validationSets[[idx]], rownames(subset2)[i+1])
}
for(i in 0:(nrow(subset3)-1)){
  idx = (i %% validationSize)+1
  validationSets[[idx]] = c(validationSets[[idx]], rownames(subset3)[i+1])
}
for(i in 0:(nrow(subset4)-1)){
  idx = (i %% validationSize)+1
  validationSets[[idx]] = c(validationSets[[idx]], rownames(subset4)[i+1])
}
for(i in 0:(nrow(subset5)-1)){
  idx = (i %% validationSize)+1
  validationSets[[idx]] = c(validationSets[[idx]], rownames(subset5)[i+1])
}
#remove NAs
for(i in 1:validationSize){
  x = validationSets[[i]]
  validationSets[[i]] = x[!is.na(x)]
}

lowerLabels = labels
rownames(lowerLabels) = tolower(rownames(lowerLabels))

Wmats = c()
Hnorms = data.frame(nrow=5)
orderings = c()

for(k in 1:validationSize){
  trainTableFold = trainTable[!(rownames(trainTable) %in% validationSets[[k]]),]
  #trainTableFold = trainTable
  #trainTableFold = trainTableFold[,-1]
  validationTableFold = trainTable[rownames(trainTable) %in% validationSets[[k]],]
  write.table(trainTableFold, file="/DataTables/currentFoldReducedTrain.txt", sep="\t", row.names = TRUE, col.names = TRUE)
  write.table(validationTableFold, file="/DataTables/currentFoldReducedValidation.txt", sep='\t', row.names = TRUE, col.names = TRUE)
  source("/src/get.classifier_Tim.R")
  g1Fixed = reorient(g1, labels)
  positions = c(g1Fixed[[2]],g1Fixed[[3]],g1Fixed[[4]],g1Fixed[[5]],g1Fixed[[6]])
  g1Fixed = g1Fixed[[1]]
  clus1Position = which(positions == 1)
  clus2Position = which(positions == 2)
  clus3Position = which(positions == 3)
  clus4Position = which(positions == 4)
  clus5Position = which(positions == 5)
  llSub = lowerLabels[rownames(lowerLabels) %in% rownames(g1Fixed),]
  llSub = llSub[rownames(g1Fixed),]
  #print(table(g1Fixed$g1 == llSub$cluster))
  #print(table(g1Fixed$g1, llSub$cluster))
  H1.norm.reoriented = H1.norm[c(clus1Position,clus2Position,clus3Position,clus4Position,clus5Position),]
  H1.norm.reoriented = data.frame(H1.norm.reoriented)
  H1.norm.reoriented = H1.norm.reoriented[rownames(g1Fixed)]
  rownames(H1.norm.reoriented) = c(1,2,3,4,5)
  Hnorms = cbind(Hnorms, H1.norm.reoriented)
  orderings = c(orderings, rownames(g1Fixed))
}
Hnorms = Hnorms[,-1]
Hnorms = rbind(Hnorms, apply(Hnorms, 2, function(x) max(x)))
Hnorms = rbind(Hnorms, apply(Hnorms, 2, function(x) which.max(x)))
Hnorms = Hnorms[,orderings]
Hnorms = t(Hnorms)
Hnorms = data.frame(Hnorms)
colnames(Hnorms) = c("1","2","3","4","5","maxVal","cluster")
Hnorms = Hnorms[order(Hnorms$maxVal, decreasing = TRUE),]
allTrainingLowerLabels = lowerLabels[rownames(lowerLabels) %in% rownames(Hnorms),]
allTrainingLowerLabels = allTrainingLowerLabels[rownames(Hnorms),]
Hnorms$correct = Hnorms$cluster == allTrainingLowerLabels$cluster
Hnorms$trueCluster = allTrainingLowerLabels$cluster



top70th = Hnorms[1:(nrow(Hnorms)*.7),]
top70thConfMat = confusionMatrix(factor(top70th$cluster), factor(top70th$trueCluster))
top70thConfMat



splitList = split(Hnorms, cut(Hnorms$maxVal, seq(0,1,by = 1/20)), drop = FALSE)
plotDF = data.frame()
averages = c()
agreements = c()
upErrs = c()
lowErrs = c()
for(i in 1:length(splitList)){
  correct = sum(splitList[[i]]$correct)
  agreement = correct/(length(splitList[[i]]$correct))
  agreements = c(agreements, agreement)
  avg = mean(splitList[[i]]$maxVal)
  averages = c(averages, avg)
  if(!is.na(correct) && (length(splitList[[i]]$correct) != 0)){
    vals = prop.test(c(correct),length(splitList[[i]]$correct),conf.level=0.84, correct=FALSE)[[6]]
    upper = vals[2]
    lower = vals[1]
    upErrs = c(upErrs,upper)
    lowErrs = c(lowErrs, lower)
  }
}

averages = averages[!is.na(averages)]
agreements = agreements[!is.na(agreements)]
upErrs = upErrs-agreements
lowErrs = agreements-lowErrs

plotDF = cbind(averages, agreements)
plotDF = data.frame(plotDF)
colnames(plotDF) = c("AverageOfBin", "AgreementFraction")
grob = grobTree(textGrob(paste("Pearson Correlation : ", 
                               round(cor.test(plotDF$AverageOfBin, 
                                              plotDF$AgreementFraction, 
                                              method="pearson")[[4]], 4) ),
                         x = 0.05, y = 0.87, 
                         hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))

grob2 = grobTree(textGrob(paste("Weighted Pearson Correlation : ", 
                                round(pearson_custom(averages,agreements,splitList), 4) ),
                          x = 0.05, y = 0.97, 
                          hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))

p = ggplot(plotDF, aes(x=AverageOfBin, y=AgreementFraction)) + geom_point() +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) + 
  annotation_custom(grob) +
  annotation_custom(grob2) +
  geom_pointrange(aes(ymin=AgreementFraction-lowErrs, ymax = AgreementFraction+upErrs), width=.2)

filePath = paste("/Plots/AgreementvsConfidencek",
                 4,"NMF","_COO.jpeg", sep="")
print(filePath)

jpeg(filePath)
print(p)
dev.off()
print(p)


source("src/GeneratedReducedCOO.R")
rownames(reducedDFCOO) = tolower(rownames(reducedDFCOO))
Hnorms = Hnorms[order(Hnorms[,7], Hnorms[,6]),]
bjoernFile = reducedDFCOO[rownames(Hnorms),]
bjoernFile = cbind(Hnorms[6:9], bjoernFile)
write.table(bjoernFile, "WrittenFiles/pNMF_sortedSamples_COO.tsv", sep='\t', row.names = TRUE, col.names = TRUE)



c1Features = c("SV.BCL6","B2M.CD70.PD1.FAS.com","BLC10.NOTCH2.A20.SPEN", "SV.TP63", "MYD88.OTHER")
c2Features = c("Chrom.Mod.Enzymes","BCL2.comp","Pi3K.Mod")
c3Features = c("Sum.CNAs.C2","TP53.biallelic","GENOME_DOUBLING","X2P16.1.AMP","X21Q.AMP","X9p21.3.DEL")
c4Features = c("CD79B.MYD88.TBL1XR1.BLIMP.ETV6.ZC3H12A","sumCNAs.C5","X18Q.AMP","MYD88.L265")
c5Features = c("SumC4.mut","NFKBalt.RASalt.JAK.STAT.alt.RHOA.SGK1.KLH6.cd58.cd83","Hist.comp",
               "COO..GCB.1.unclass.2.ABC.3.na.")
featureOrder = c(c1Features,c2Features,c3Features,c4Features,c5Features)
#make the heatmap
heatmapDF = bjoernFile[,5:ncol(bjoernFile)]
heatmapDF = bjoernFile[,rev(featureOrder)]
heatmapDFMatrix = standardize(heatmapDF)
heatmapDFMatrix = as.matrix(heatmapDFMatrix)


#image(heatmapDFMatrix, useRaster=TRUE, col=gray.colors(20,start=1,end=0,gamma = 2.2, alpha=NULL))
melted = melt(heatmapDFMatrix)
v1 = which(bjoernFile$cluster == 1)[1]
v2 = which(bjoernFile$cluster == 2)[1]
v3 = which(bjoernFile$cluster == 3)[1]
v4 = which(bjoernFile$cluster == 4)[1]
v5 = which(bjoernFile$cluster == 5)[1]

b1 = rownames(bjoernFile)[v1]
b2 = rownames(bjoernFile)[v2]
b3 = rownames(bjoernFile)[v3]
b4 = rownames(bjoernFile)[v4]
b5 = rownames(bjoernFile)[v5]

y1 = length(featureOrder)-length(c1Features)
y2 = y1-length(c2Features)
y3 = y2-length(c3Features)
y4 = y3-length(c4Features)

melted$agreeVec = rep(bjoernFile$correct, ncol(heatmapDFMatrix))

p1 = ggplot(data = melted, aes(x=Var1, y=Var2)) + 
  geom_tile(aes(fill=agreeVec, alpha = value), colour = "white") +
  scale_x_discrete(expand=c(0,0), breaks=c(b1,b2,b3,b4,b5)
                   , labels = c("C1","C2","C3","C4","C5")) +
  scale_y_discrete(expand=c(0,0)) +
  geom_vline(xintercept=c(v2,v3,v4,v5), linetype="solid") +
  geom_hline(yintercept=c(y1,y2,y3,y4), linetype = "solid")
p1
jpeg("/Plots/pNMF_sortedHeatmap_COO.jpeg", width = 1920, height = 1080)
print(p1)
dev.off()
