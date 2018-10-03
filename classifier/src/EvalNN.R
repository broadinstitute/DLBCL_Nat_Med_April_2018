
rm(list = ls())
toeval = FALSE
knitr::opts_chunk$set(echo = TRUE)
source("src/PearsonCustom.R")
source("src/StandardizeReducedFeatures.R")
source("src/LoadLibraries.R")

downSample = FALSE
useSV = TRUE
useCNA = TRUE
GD = FALSE
useProbs = TRUE
fixedValidationSets = FALSE
fullFeatures = TRUE

filename = "ANN/WrittenFiles/4foldNN"
if(downSample){
  filename = paste(filename, "_DS", sep="")
}
if(!useSV){
  filename = paste(filename, "_noSV", sep="")
}
if(!useCNA){
  filename = paste(filename, "_noCNA", sep="")
}
if(GD){
  filename = paste(filename, "_GD", sep="")
}
if(!useProbs){
  filename = paste(filename, "_BinaryLabels", sep="")
}
if(fixedValidationSets){
  filename = paste(filename, "_FixedSets", sep="")
}
if(fullFeatures){
  filename = paste(filename,"_161F", sep="")
}
filename = paste(filename, ".csv", sep="")

nnDF = read.table(filename, sep=',', header=TRUE)


#numBins = 20
nnDF = nnDF[,-1]
nnDF$correct = nnDF$correct == "True"
splitList = split(nnDF, cut(seq_along(rownames(nnDF)), 10, labels = FALSE))



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

#averages = log(averages)
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
                                weightedPearson),
                          x = 0.05, y = 0.97, 
                          hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))


p = ggplot(plotDF, aes(x=AverageOfBin, y=AgreementFraction)) + geom_point() +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) + 
  annotation_custom(grob) +
  annotation_custom(grob2) +
  geom_pointrange(aes(ymin=AgreementFraction-lowErrs, ymax = AgreementFraction+upErrs), width=.2)

filename = paste("Plots/AgreementvsConfidencek",
                 4,"NN", sep="")
if(downSample){
  filename = paste(filename, "_DS", sep="")
}
if(!useSV){
  filename = paste(filename, "_noSV", sep="")
}
if(!useCNA){
  filename = paste(filename, "_noCNA", sep="")
}
if(GD){
  filename = paste(filename, "_GD", sep="")
}
if(!useProbs){
  filename = paste(filename, "_BinaryLabels", sep="")
}
if(fullFeatures){
  filename = paste(filename,"_161F", sep="")
}
filename = paste(filename, ".jpeg", sep="")
print(filename)

jpeg(filename)
print(p)
dev.off()




#top 70th percentile agreement
percentile = .7
top70th = nnDF[(nrow(nnDF)*(1-percentile)):nrow(nnDF),]
top70thConfMat = confusionMatrix(factor(top70th$predictedCluster), factor(top70th$trueCluster))
top70thConfMat



pearson_custom(averages,agreements,splitList)
source("src/ConstructReducedFeatureMatrix.R")
rownames(nnDF) = tolower(nnDF$Sample)
rownames(reducedDF) = tolower(rownames(reducedDF))

bjoernFile = reducedDF[rownames(nnDF),]
bjoernFile = 
  cbind(nnDF$maxVal,nnDF$predictedCluster, nnDF$correct, nnDF$trueCluster, bjoernFile)
colnames(bjoernFile)[1:4] = c("maxVal", "cluster", "correct", "trueCluster")
bjoernFile = bjoernFile[order(bjoernFile[,2], bjoernFile[,1]),]
filename = "WrittenFiles/NN_sortedSamples"
if(downSample){
  filename = paste(filename, "_DS", sep="")
}
if(!useSV){
  filename = paste(filename, "_noSV", sep="")
}
if(!useCNA){
  filename = paste(filename, "_noCNA", sep="")
}
if(GD){
  filename = paste(filename, "_GD", sep="")
}
if(!useProbs){
  filename = paste(filename, "_BinaryLabels", sep="")
}
if(fullFeatures){
  filename = paste(filename,"_161F", sep="")
}
filename = paste(filename, ".tsv", sep="")
write.table(bjoernFile, filename, sep='\t', row.names = TRUE, col.names = TRUE)


c1Features = toupper(c("SV.BCL6","B2M.CD70.PD1.FAS.com","BLC10.NOTCH2.A20.SPEN", "SV.TP63", "MYD88.OTHER"))
c2Features = toupper(c("Chrom.Mod.Enzymes","BCL2.comp","Pi3K.Mod"))
c3Features = toupper(c("Sum.CNAs.C2","TP53.biallelic","GENOME_DOUBLING","X2P16.1.AMP","X21Q.AMP","X9p21.3.DEL"))
c4Features = toupper(c("CD79B.MYD88.TBL1XR1.BLIMP.ETV6.ZC3H12A","sumCNAs.C5","X18Q.AMP","MYD88.L265"))
c5Features = toupper(c("SumC4.mut","NFKBalt.RASalt.JAK.STAT.alt.RHOA.SGK1.KLH6.cd58.cd83","Hist.comp"))
c1Features = c1Features[c1Features %in% colnames(reducedDF)]
c2Features = c2Features[c2Features %in% colnames(reducedDF)]
c3Features = c3Features[c3Features %in% colnames(reducedDF)]
c4Features = c4Features[c4Features %in% colnames(reducedDF)]
c5Features = c5Features[c5Features %in% colnames(reducedDF)]
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

filename = "Plots/NN_sortedHeatmap"
if(downSample){
  filename = paste(filename, "_DS", sep="")
}
if(!useSV){
  filename = paste(filename, "_noSV", sep="")
}
if(!useCNA){
  filename = paste(filename, "_noCNA", sep="")
}
if(GD){
  filename = paste(filename, "_GD", sep="")
}
if(!useProbs){
  filename = paste(filename, "_BinaryLabels", sep="")
}
if(fullFeatures){
  filename = paste(filename,"_161F", sep="")
}
filename = paste(filename, ".jpeg", sep="")
# jpeg(filename, width = 1920, height = 1080)
# print(p1)
# dev.off()

evaluateTrashSet = TRUE
if(evaluateTrashSet){
  source("src/GenerateTrashReduced.R")
  write.table(reducedDFTrash, "ANN/DataTables/reducedDFTrash.tsv", row.names = TRUE, col.names = TRUE, sep = '\t')
  setwd("ANN")
  system("/Users/twood/anaconda3/bin/python testTrashSet.py")
  nnTrashDF = read.csv("WrittenFiles/NN_TrashSet_Confidences.tsv", header = TRUE, sep='\t')
  nnTrashDF = nnTrashDF[,-1]
  setwd("..")
  
  maxvalsTrash = c()
  for(i in 1:nrow(nnTrashDF)){
    val = max(nnTrashDF[i,])
    maxvalsTrash = c(maxvalsTrash, val)
  }
  nnTrashDF = cbind(nnTrashDF, maxvalsTrash)
  
  p = ggplot(nnTrashDF, aes(x=maxvalsTrash)) +
    geom_histogram(fill="black", alpha=1, position="identity", bins = 100) +
    ggtitle("Confidence in NN: Trash Set") +
    theme(plot.title = element_text(hjust = 0.5),
          title = element_text(size=20)) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by=.1)) +
    theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
          axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))
  jpeg("Plots/NN_TrashSet_Confidences.jpeg", width = 1080, height = 720)
  print(p)
  dev.off()
  
  p = ggplot(bjoernFile, aes(x=maxVal)) +
    geom_histogram(fill="black", alpha=1, position="identity", bins = 100) +
    ggtitle("Confidence in NN: All Validation Sets") +
    theme(plot.title = element_text(hjust = 0.5),
          title = element_text(size=20)) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by=.1)) +
    theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
          axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))
  jpeg("Plots/NN_ValidationSets_Confidences.jpeg", width = 1080, height = 720)
  print(p)
  dev.off()
}


