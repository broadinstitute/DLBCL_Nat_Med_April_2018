
rm(list = ls())
set.seed(200)
knitr::opts_chunk$set(echo = TRUE)
toeval = TRUE
useSV = TRUE
useCNA = TRUE
GD = TRUE
fixedValidationSets = FALSE
fullFeatures = TRUE
if(fullFeatures){
  nFeatures = 161
} else {
  nFeatures = 21
}

source("src/LoadLibraries.R")
source("src/GenerateLabels.R")
source("src/ConstructReducedFeatureMatrix.R")
source("src/GenerateTestTrainSets.R")
source("src/reorientG1.R")
source("src/PearsonCustom.R")
source("src/StandardizeReducedFeatures.R")



downSample = FALSE
orderedLabels = labels[rownames(reducedDF),]
reducedDF = cbind(factor(orderedLabels$cluster),reducedDF)
if(downSample){
  source("src/GeneratedReducedCOO.R")
  reducedDF = reducedDF[rownames(reducedDF) %in% rownames(reducedDFCOO),]
}
colnames(reducedDF)[1] = "cluster"
validationSize = 4
validationSets = vector("list", validationSize)
trainTable = reducedDF[rownames(reducedDF) %in% trainingSet,]
if(fullFeatures){
  source("src/GenerateFullDF2.R")
  orderedLabels = labels[rownames(fullDF),]
  fullDF = cbind(factor(orderedLabels$cluster), fullDF)
  colnames(fullDF)[1] = "cluster"
  trainTable = fullDF[rownames(fullDF) %in% trainingSet,]
}
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
#validationSize = 1

if(fixedValidationSets){
  v1 = read.csv("WrittenFiles/validationSet.1.txt", header = FALSE)
  v2 = read.csv("WrittenFiles/validationSet.2.txt", header = FALSE)
  v3 = read.csv("WrittenFiles/validationSet.3.txt", header = FALSE)
  v4 = read.csv("WrittenFiles/validationSet.4.txt", header = FALSE)
  validationSets = list(v1$V1,v2$V1,v3$V1,v4$V1)
}
lowerLabels = labels
rownames(lowerLabels) = tolower(rownames(lowerLabels))

Wmats = c()
Hnorms = data.frame(nrow=5)
orderings = c()
allTrashSetConfidences = c()
avgMat = NA

for(k in 1:validationSize){
  trainTableFold = trainTable[!(rownames(trainTable) %in% validationSets[[k]]),]
  #trainTableFold = trainTable
  #trainTableFold = trainTableFold[,-1]
  validationTableFold = trainTable[rownames(trainTable) %in% validationSets[[k]],]
  write.table(trainTableFold, file="DataTables/currentFoldReducedTrain.txt", sep="\t", row.names = TRUE, col.names = TRUE)
  write.table(validationTableFold, file="DataTables/currentFoldReducedValidation.txt", sep='\t', row.names = TRUE, col.names = TRUE)
  source("src/get.classifier_Tim.R")
  g1Fixed = reorient(g1, W0, !fullFeatures)
  clus1Position = g1Fixed[[2]]
  clus2Position = g1Fixed[[3]]
  clus3Position = g1Fixed[[4]]
  clus4Position = g1Fixed[[5]]
  clus5Position = g1Fixed[[6]]
  g1Fixed = g1Fixed[[1]]
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
  
  evaluateTrashSet = TRUE
  if(evaluateTrashSet){
    source("src/GenerateTrashReduced.R")
    if(fullFeatures){
      predictDF = t(fullDFTrash)
    } else {
      predictDF = t(reducedDFTrash)
    }
    trashRes <- NMF.W(as.matrix(predictDF),as.matrix(W0),tol,K)
    H1Trash <- trashRes[[2]] #### H1Trash is an association of new samples to clustering
    H1Trash.eps = apply(H1Trash,2,function(x) if(sum(x) == 0){x = x+.1} else {x})
    H1Trash.norm <- apply(H1Trash.eps,2,function(x) x/sum(x))
    H1Trash.norm = H1Trash.norm[c(clus1Position,clus2Position,clus3Position,clus4Position,clus5Position),]
    H1Trash.norm = as.matrix(H1Trash.norm)
    
    maxvalsTrash = c()
    clusters = c()
    for(i in 1:ncol(H1Trash.norm)){
      val = max(H1Trash.norm[,i])
      clus = which(max(H1Trash.norm[,i]) == H1Trash.norm[,i])
      maxvalsTrash = c(maxvalsTrash, val)
      clusters = c(clusters, clus)
    }
    maxvalsTrash = data.frame(maxvalsTrash)
    allTrashSetConfidences = c(allTrashSetConfidences, maxvalsTrash)
    if(k == 1){
      avgMat = H1Trash.norm
    } else {
      avgMat = avgMat + H1Trash.norm
    }
  }
}

avgMat = avgMat/validationSize
avgMat = t(avgMat)
maxVals = c()
cluster = c()
for(i in 1:nrow(avgMat)){
  maxVal = max(avgMat[i,1:5])
  clus = which(avgMat[i,1:5] == maxVal)
  maxVals = c(maxVals, maxVal)
  cluster = c(cluster, clus)
}

avgMat = cbind(avgMat, maxVals, cluster)

if(evaluateTrashSet){
  allTrashSetConfidences = data.frame(avgMat[,"maxVals"])
  colnames(allTrashSetConfidences)[1] = "confidences"
  p = ggplot(allTrashSetConfidences, aes(x=confidences)) +
    geom_histogram(fill="black", alpha=1, position="identity", bins = 100) +
    ggtitle(paste("Confidence in NMF: Null Set, Features ",nFeatures, sep="")) +
    theme(plot.title = element_text(hjust = 0.5),
          title = element_text(size=20)) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by=.1)) +
    theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
          axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))
  fn = paste("Plots/NMF_TrashSet_Confidences_Features",nFeatures,".jpeg")
  jpeg(fn, width = 1080, height = 720)
  print(p)
  dev.off()
  
  clus1Order = c("SV.BCL6", "BCL10", "TNFAIP3", "UBE2A", "CD70", "B2M", "NOTCH2", "TMEM30A", "FAS",
                 "X5p.AMP", "SV.TP63", "ZEB2", "HLA.B", "SPEN", "SV.CD274.PDCD1LG2")
  clus3Order = c("TP53", "X17p.DEL", "X17p11.2.DEL", "X21q.AMP", "X9p21.3.DEL", "X9q21.13.DEL",
                 "X4q35.1.DEL","X1p31.1.DEL", "X1p36.11.DEL", "X1p13.1.DEL", "X4q21.22.DEL", "X14q32.31.DEL",
                 "X3p21.31.DEL", "X2p16.1.AMP", "X16q12.1.DEL", "X1p36.32.DEL", "X3q28.DEL", "X1q23.3.AMP" ,
                 "X18q23.DEL", "X8q24.22.AMP", "X17q24.3.AMP", "X13q14.2.DEL", "X19p13.3.DEL", "X5q.AMP",
                 "X11q.AMP", "X13q31.3.AMP", "X6p.AMP", "X2q22.2.DEL", "X12p13.2.DEL", "X6q.DEL",
                 "X3q28.AMP", "X11q23.3.AMP", "X1q42.12.DEL", "X8q12.1.DEL", "X19q13.32.DEL", "X10q23.31.DEL")
  
  clus2Order = c("BCL2", "SV.BCL2", "CREBBP", "EZH2", "KMT2D", "TNFRSF14", "HVCN1", "IRF8", "GNA13", "MEF2B",
                 "PTEN")
  
  clus5Order = c("SGK1", "HIST1H1E", "NFKBIE","BRAF", "CD83", "NFKBIA", "CD58", "HIST1H2BC", "STAT3",
                 "HIST1H1C", "ZFP36L1", "KLHL6", "HIST1H1D", "HIST1H1B", "ETS1", "TOX", "HIST1H2AM",
                 "HIST1H2BK", "RHOA", "ACTB", "LTB", "SF3B1", "CARD11", "HIST1H2AC")
  
  clus4Order = c("X18q.AMP", "X13q.AMP", "CD79B", "X3p.AMP", "MYD88", "ETV6", "X18p.AMP", "PIM1", 
                 "X17q25.1.DEL", "TBL1XR1", "X19q13.42.AMP", "GRHPR", "ZC3H12A", "X19p13.2.DEL",
                 "X19q.AMP", "HLA.A", "PRDM1", "BTG1", "X18q21.33.BCL2..AMP", "SV.MYC")
  
  rownames(fullDFTrash) = make.names(rep(paste( LETTERS, "row", sep =""), length.out=nrow(fullDFTrash)), unique = TRUE)
  clusterConfidencesDF = data.frame(confidence = allTrashSetConfidences$confidences, predictedCluster = cluster)
  rownames(clusterConfidencesDF) = rownames(fullDFTrash)
  colnames(clusterConfidencesDF)[1] = "confidence"
  #subset by confidence first
  highConfidencesDF = clusterConfidencesDF[order(clusterConfidencesDF$confidence),]
  highConfidencesDF = highConfidencesDF[(nrow(highConfidencesDF)-500):(nrow(highConfidencesDF)),]
  highConfidencesDF = highConfidencesDF[order(highConfidencesDF$predictedCluster,
                                              highConfidencesDF$confidence),]
  
  sortedOrder = rownames(highConfidencesDF)
  
  featureOrder = toupper(c(clus1Order, clus2Order, clus3Order, clus4Order, clus5Order))
  heatmapDF = fullDFTrash[sortedOrder,featureOrder]
  heatmapDF = heatmapDF[,rev(colnames(heatmapDF))]
  #heatmapDF = cbind(highConfidencesDF$predictedCluster, heatmapDF)
  #colnames(heatmapDF)[1] = "clusters"
  heatmapDFMatrix = as.matrix(heatmapDF)
  
  melted = melt(heatmapDFMatrix)
  c1 = which(highConfidencesDF$predictedCluster == 1)[1]
  c2 = which(highConfidencesDF$predictedCluster == 2)[1]
  c3 = which(highConfidencesDF$predictedCluster == 3)[1]
  c4 = which(highConfidencesDF$predictedCluster == 4)[1]
  c5 = which(highConfidencesDF$predictedCluster == 5)[1]
  # 
  b1 = rownames(heatmapDF)[c1]
  b2 = rownames(heatmapDF)[c2]
  b3 = rownames(heatmapDF)[c3]
  b4 = rownames(heatmapDF)[c4]
  b5 = rownames(heatmapDF)[c5]
  
  y1 = length(featureOrder)-length(clus1Order)
  y2 = y1-length(clus2Order)
  y3 = y2-length(clus3Order)
  y4 = y3-length(clus4Order)
  title = paste("High in NMF: Null Set, Top 500, Features ",nFeatures,sep="")
  p1 = ggplot(data = melted, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile(colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") +
    scale_x_discrete(expand=c(0,0), breaks=c(b1,b2,b3,b4,b5)
                     , labels = c("C1", "C2", "C3", "C4", "C5")) +
    geom_vline(xintercept=c(c2,c3,c4,c5), linetype="solid") +
    geom_hline(yintercept=c(y1,y2,y3,y4), linetype = "solid") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5),
          title = element_text(size=18))
  
  fn = paste("Plots/TrashHeatmapNMF_HC_Features",nFeatures,".jpeg", sep="")
  jpeg(fn, width = 1080, height = 720)
  print(p1)
  dev.off()
  
  
  lowConfidencesDF = clusterConfidencesDF[order(clusterConfidencesDF$confidence),]
  lowConfidencesDF = lowConfidencesDF[1:500,]
  lowConfidencesDF = lowConfidencesDF[order(lowConfidencesDF$predictedCluster,
                                            lowConfidencesDF$confidence),]
  
  sortedOrder = rownames(lowConfidencesDF)
  
  featureOrder = toupper(c(clus1Order, clus2Order, clus3Order, clus4Order, clus5Order))
  heatmapDF = fullDFTrash[sortedOrder,featureOrder]
  heatmapDF = heatmapDF[,rev(colnames(heatmapDF))]
  #heatmapDF = cbind(lowConfidencesDF$predictedCluster, heatmapDF)
  #colnames(heatmapDF)[1] = "clusters"
  heatmapDFMatrix = as.matrix(heatmapDF)
  
  melted = melt(heatmapDFMatrix)
  c1 = which(lowConfidencesDF$predictedCluster == 1)[1]
  c2 = which(lowConfidencesDF$predictedCluster == 2)[1]
  c3 = which(lowConfidencesDF$predictedCluster == 3)[1]
  c4 = which(lowConfidencesDF$predictedCluster == 4)[1]
  c5 = which(lowConfidencesDF$predictedCluster == 5)[1]
  # 
  b1 = rownames(heatmapDF)[c1]
  b2 = rownames(heatmapDF)[c2]
  b3 = rownames(heatmapDF)[c3]
  b4 = rownames(heatmapDF)[c4]
  b5 = rownames(heatmapDF)[c5]
  
  y1 = length(featureOrder)-length(clus1Order)
  y2 = y1-length(clus2Order)
  y3 = y2-length(clus3Order)
  y4 = y3-length(clus4Order)
  title = paste("Low Confidence in NMF: Null Set, Bottom 500, Features ",nFeatures,sep="")
  p1 = ggplot(data = melted, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile(colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") +
    scale_x_discrete(expand=c(0,0), breaks=c(b1,b2,b3,b4,b5)
                     , labels = c("C1", "C2", "C3", "C4", "C5")) +
    geom_vline(xintercept=c(c2,c3,c4,c5), linetype="solid") +
    geom_hline(yintercept=c(y1,y2,y3,y4), linetype = "solid") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5),
          title = element_text(size=18))
  
  fn = paste("Plots/TrashHeatmapNMF_LC_Features",nFeatures,".jpeg", sep="")
  jpeg(fn, width = 1080, height = 720)
  print(p1)
  dev.off()
}

Hnorms = Hnorms[,-1]
toBind1 = apply(Hnorms, 2, function(x) max(x))
toBind2 = apply(Hnorms, 2, function(x) which.max(x))
toBind2 = unlist(toBind2)
Hnorms = rbind(Hnorms, toBind1)
Hnorms = rbind(Hnorms, toBind2)
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



splitList = split(Hnorms, cut(seq_along(rownames(Hnorms)), 10, labels = FALSE))
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

weightedPearson = round(pearson_custom(averages,agreements,splitList),4)
grob2 = grobTree(textGrob(paste("Weighted Pearson Correlation : ", 
                                weightedPearson ),
                          x = 0.05, y = 0.97, 
                          hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



p = ggplot(plotDF, aes(x=AverageOfBin, y=AgreementFraction)) + geom_point() +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) + 
  annotation_custom(grob) +
  annotation_custom(grob2) +
  geom_pointrange(aes(ymin=AgreementFraction-lowErrs, ymax = AgreementFraction+upErrs), width=.2)

filePath = paste("Plots/AgreementvsConfidencek",
                 4,"NMF", sep="")
if(downSample){
  filePath = paste(filepath,"_DS", sep="")
}
if (!useSV){
  filePath = paste(filePath,"_noSV", sep="")
}
if (!useCNA){
  filePath = paste(filePath,"_noCNA", sep="")
}
if (GD){
  filePath = paste(filePath,"_GD", sep="")
}
if(fixedValidationSets){
  filePath = paste(filePath,"_FixedSets", sep="")
}
if(fullFeatures){
  filePath = paste(filePath,"_161F", sep="")
}
filePath = paste(filePath, ".jpeg", sep="")

print(filePath)

#jpeg(filePath)
#print(p)
#dev.off()
#print(p)


source("src/ConstructReducedFeatureMatrix.R")
rownames(reducedDF) = tolower(rownames(reducedDF))
Hnorms = Hnorms[order(Hnorms[,7], Hnorms[,6]),]
bjoernFile = reducedDF[rownames(Hnorms),]
bjoernFile = cbind(Hnorms[6:9], bjoernFile)
filename = "WrittenFiles/pNMF_sortedSamples"

if(downSample){
  filename = paste(filename, "_DS", sep = "")
}
if(!useSV){
  filename = paste(filename, "_noSV", sep = "")
}
if(!useCNA){
  filename = paste(filename, "_noCNA", sep = "")
}
if(GD){
  filename = paste(filename, "_GD", sep = "")
}
if(fixedValidationSets){
  filename = paste(filename,"_FixedSets", sep="")
}
if(fullFeatures){
  filename = paste(filename,"_161F", sep="")
}
filename = paste(filename, ".tsv", sep = "")
write.table(bjoernFile, filename, sep='\t', row.names = TRUE, col.names = TRUE)

c1Features = toupper(c("SV.BCL6","B2M.CD70.PD1.FAS.com","BLC10.NOTCH2.A20.SPEN", "SV.TP63", "MYD88.OTHER"))
c2Features = toupper(c("Chrom.Mod.Enzymes","BCL2.comp","Pi3K.Mod"))
c3Features = toupper(c("Sum.CNAs.C2","TP53.biallelic","GENOME_DOUBLING","X2P16.1.AMP","X21Q.AMP","X9p21.3.DEL"))
c4Features = toupper(c("CD79B.MYD88.TBL1XR1.BLIMP.ETV6.ZC3H12A","sumCNAs.C5","X18Q.AMP","MYD88.L265"))
c5Features = toupper(c("SumC4.mut","NFKBalt.RASalt.JAK.STAT.alt.RHOA.SGK1.KLH6.cd58.cd83","Hist.comp"))
c1Features = c1Features[c1Features %in% colnames(bjoernFile)]
c2Features = c2Features[c2Features %in% colnames(bjoernFile)]
c3Features = c3Features[c3Features %in% colnames(bjoernFile)]
c4Features = c4Features[c4Features %in% colnames(bjoernFile)]
c5Features = c5Features[c5Features %in% colnames(bjoernFile)]

featureOrder = c(c1Features,c2Features,c3Features,c4Features,c5Features)
featureOrder = toupper(featureOrder)
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
filename = "Plots/pNMF_sortedHeatmap"

if(downSample){
  filename = paste(filename, "_DS", sep = "")
}
if(!useSV){
  filename = paste(filename, "_noSV", sep = "")
}
if(!useCNA){
  filename = paste(filename, "_noCNA", sep = "")
}
if(GD){
  filename = paste(filename, "_GD", sep = "")
}
if(fixedValidationSets){
  filename = paste(filename,"_FixedSets", sep="")
}
if(fullFeatures){
  filename = paste(filename,"_161F", sep="")
}
filename = paste(filename, ".jpeg", sep = "")

#jpeg(filename, width = 1920, height = 1080)
#print(p1)
#dev.off()

upper = 15
title = paste("Confidence in NMF: All Validation Sets", " Features ",nFeatures, sep="")
p = ggplot(bjoernFile, aes(x=maxVal)) +
  geom_histogram(fill="black", alpha=1, position="identity", bins = 100) +
  ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by=.1)) +
  scale_y_continuous(limits = c(0,upper), breaks = seq(0,upper, by=1)) +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))
fn = paste("Plots/NMF_ValidationSets_Confidences_Features",nFeatures,".jpeg")
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()