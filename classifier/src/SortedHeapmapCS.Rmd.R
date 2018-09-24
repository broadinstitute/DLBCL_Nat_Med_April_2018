rm(list = ls())
knitr::opts_chunk$set(echo = TRUE)
library(reshape2)
library(ggplot2)
library(gridExtra)

source("src/GenerateFullDF.R")
source("src/addBaselines.R")
source("src/SumRankWmats.R")
source("src/GenerateProbabilityTable.R")
source("src/addCSBaselinesV2.R")

baselineProbTable = csBaselines(fullDF, connectivityProbs1.0, labels, summedW)
c1Subset = baselineProbTable[baselineProbTable$cluster == 1,]
c2Subset = baselineProbTable[baselineProbTable$cluster == 2,]
c3Subset = baselineProbTable[baselineProbTable$cluster == 3,]
c4Subset = baselineProbTable[baselineProbTable$cluster == 4,]
c5Subset = baselineProbTable[baselineProbTable$cluster == 5,]


#sample order already computed in csBaselines function
#get feature order
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


fullFeatureOrder = c(clus1Order, clus2Order, clus3Order, clus4Order, clus5Order)
fullFeatureOrder = toupper(fullFeatureOrder)

filteredDF = fullDF[rownames(baselineProbTable),fullFeatureOrder]



#make the heatmap
heatmapDF = filteredDF[,rev(colnames(filteredDF))]
heatmapDFMatrix = as.matrix(heatmapDF)

#image(heatmapDFMatrix, useRaster=TRUE, col=gray.colors(20,start=1,end=0,gamma = 2.2, alpha=NULL))
melted = melt(heatmapDFMatrix)
v1 = which(baselineProbTable$cluster == 1)[1]
v2 = which(baselineProbTable$cluster == 2)[1]
v3 = which(baselineProbTable$cluster == 3)[1]
v4 = which(baselineProbTable$cluster == 4)[1]
v5 = which(baselineProbTable$cluster == 5)[1]

b1 = rownames(baselineProbTable)[v1]
b2 = rownames(baselineProbTable)[v2]
b3 = rownames(baselineProbTable)[v3]
b4 = rownames(baselineProbTable)[v4]
b5 = rownames(baselineProbTable)[v5]

y1 = length(fullFeatureOrder)-length(clus1Order)
y2 = y1-length(clus2Order)
y3 = y2-length(clus3Order)
y4 = y3-length(clus4Order)
p1 = ggplot(data = melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") +
  scale_x_discrete(expand=c(0,0), breaks=c(b1,b2,b3,b4,b5, "DLBCL_C_D_1141_NULLPAIR")
                   , labels = c("C1","C2","C3","C4","C5", "S")) +
  geom_vline(xintercept=c(v2,v3,v4,v5), linetype="solid") +
  geom_hline(yintercept=c(y1,y2,y3,y4), linetype = "solid")
print(p1)
jpeg("Plots/heatmapClusterSpecificBaselinesAlpha.0.9.Power1.2.constantScalar.jpeg", width = 1920, height = 1080)
p1
dev.off()

#compute max vals for histogram
maxvals = c()
for(i in 1:nrow(baselineProbTable)){
  val = max(baselineProbTable[i,1:5])
  maxvals = c(maxvals,val)
}
baselineProbTable = cbind(baselineProbTable, maxvals)
colnames(baselineProbTable)[ncol(baselineProbTable)] = "maxval"
p = ggplot(data=baselineProbTable, aes(baselineProbTable$maxval)) + geom_histogram(bins=100)
jpeg("Plots/histogramProbDistsScalar0.5.Power1.2.jpeg")
p
dev.off()
p
print(sum(baselineProbTable$maxval > .90))



#write tables for Bjoern
clus = baselineProbTable[rownames(baselineProbTable) %in% rownames(heatmapDFMatrix),]
clus = clus[rownames(heatmapDFMatrix),]
heatmapDFMatrix = cbind(rownames(heatmapDFMatrix), clus$cluster, maxvals, clus$P1, clus$P2, clus$P3, clus$P4, clus$P5,heatmapDFMatrix)
colnames(heatmapDFMatrix)[1:8] = c("sample","cluster","maxprob", "P1", "P2", "P3", "P4", "P5")
write.table(heatmapDFMatrix,"WrittenFiles/heatmapMatrixAlpha.0.9.Power1.2.constantScalar.txt",sep='\t',row.names = FALSE, col.names = TRUE)
heatmapDFMatrix = data.frame(heatmapDFMatrix)

qp = qplot(seq_along(heatmapDFMatrix$maxprob), as.numeric(as.character(heatmapDFMatrix$maxprob))) + scale_y_continuous()
jpeg("Plots/histogramProbDots.Scalar0.5.Power1.2.jpeg")
qp
dev.off()
