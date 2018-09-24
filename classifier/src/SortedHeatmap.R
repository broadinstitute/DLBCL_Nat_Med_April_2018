rm(list = ls())
knitr::opts_chunk$set(echo = TRUE)
library(reshape2)
library(ggplot2)

source("src/GenerateFullDF.R")
source("src/addBaselines.R")
source("src/SumRankWmats.R")
source("src/GenerateProbabilityTable.R")

baselineProbTable = nBaseline(fullDF, connectivityProbs1.0, summedW)

c1Subset = baselineProbTable[baselineProbTable$cluster == 1,]
c2Subset = baselineProbTable[baselineProbTable$cluster == 2,]
c3Subset = baselineProbTable[baselineProbTable$cluster == 3,]
c4Subset = baselineProbTable[baselineProbTable$cluster == 4,]
c5Subset = baselineProbTable[baselineProbTable$cluster == 5,]

c1Subset = c1Subset[order(c1Subset$maxval),]
c2Subset = c2Subset[order(c2Subset$maxval),]
c3Subset = c3Subset[order(c3Subset$maxval),]
c4Subset = c4Subset[order(c4Subset$maxval),]
c5Subset = c5Subset[order(c5Subset$maxval),]
fullOrdered = rbind(c1Subset, c2Subset, c3Subset, c4Subset, c5Subset)

#get feature order
clus1Order = c("SV.BCL6", "BCL10", "TNFAIP3", "UBE2A", "CD70", "B2M", "NOTCH2", "TMEM30A", "FAS",
               "X5p.AMP", "SV.TP63", "ZEB2", "HLA.B", "SPEN", "SV.CD274.PDCD1LG2")
clus2Order = c("TP53", "X17p.DEL", "X21q.AMP", "X9p21.3.DEL", "X9q21.13.DEL",
               "X4q35.1.DEL","X1p31.1.DEL", "X1p36.11.DEL", "X1p13.1.DEL", "X4q21.22.DEL", "X14q32.31.DEL",
               "X3p21.31.DEL", "X2p16.1.AMP", "X16q12.1.DEL", "X1p36.32.DEL", "X3q28.DEL", "X1q23.3.AMP" ,
               "X18q23.DEL", "X8q24.22.AMP", "X17q24.3.AMP", "X13q14.2.DEL", "X19p13.3.DEL", "X5q.AMP",
               "X11q.AMP", "X13q34.DEL","X11p.AMP", "X13q31.3.AMP", "X6p.AMP", "X2q22.2.DEL", "X12p13.2.DEL", "X6q.DEL",
               "X3q28.AMP", "X11q23.3.AMP", "X1q42.12.DEL", "X8q12.1.DEL", "X19q13.32.DEL", "X10q23.31.DEL")

clus3Order = c("BCL2", "SV.BCL2", "CREBBP", "EZH2", "KMT2D", "TNFRSF14", "HVCN1", "IRF8", "GNA13", "MEF2B",
               "PTEN")

clus4Order = c("SGK1", "HIST1H1E", "NFKBIE","BRAF", "CD83", "NFKBIA", "CD58", "HIST1H2BC", "STAT3",
               "HIST1H1C", "ZFP36L1", "KLHL6", "HIST1H1D", "HIST1H1B", "ETS1", "TOX", "HIST1H2AM",
               "HIST1H2BK", "RHOA", "ACTB", "LTB", "SF3B1", "CARD11", "HIST1H2AC")

clus5Order = c("X18q.AMP", "X13q.AMP", "CD79B", "X3p.AMP", "MYD88", "ETV6", "X18p.AMP", "PIM1", 
               "X17q25.1.DEL", "TBL1XR1", "X19q13.42.AMP", "GRHPR", "ZC3H12A", "X19p13.2.DEL",
               "X19q.AMP", "HLA.A", "PRDM1", "BTG1", "X18q21.33.BCL2..AMP", "SV.MYC")


fullFeatureOrder = c(clus1Order, clus2Order, clus3Order, clus4Order, clus5Order)
fullFeatureOrder = toupper(fullFeatureOrder)

filteredDF = fullDF[rownames(fullOrdered),fullFeatureOrder]


#make the heatmap
heatmapDF = filteredDF[,rev(colnames(filteredDF))]
heatmapDFMatrix = as.matrix(heatmapDF)

#image(heatmapDFMatrix, useRaster=TRUE, col=gray.colors(20,start=1,end=0,gamma = 2.2, alpha=NULL))
melted = melt(heatmapDFMatrix)
ggplot(data = melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") +
  scale_x_discrete(expand=c(0,0), breaks=c(rownames(c1Subset)[[1]]
                                           , rownames(c2Subset)[[1]]
                                           , rownames(c3Subset)[[1]]
                                           , rownames(c4Subset)[[1]]
                                           , rownames(c5Subset)[[1]])
                   , labels = c("C1","C2","C3","C4","C5"))



hc1 = c1Subset[(nrow(c1Subset)*.5):nrow(c1Subset),]
hc2 = c2Subset[(nrow(c2Subset)*.5):nrow(c2Subset),]
hc3 = c3Subset[(nrow(c3Subset)*.5):nrow(c3Subset),]
hc4 = c4Subset[(nrow(c4Subset)*.5):nrow(c4Subset),]
hc5 = c5Subset[(nrow(c5Subset)*.5):nrow(c5Subset),]

lc1 = c1Subset[!(rownames(c1Subset) %in% rownames(hc1)),]
lc2 = c2Subset[!(rownames(c2Subset) %in% rownames(hc2)),]
lc3 = c3Subset[!(rownames(c3Subset) %in% rownames(hc3)),]
lc4 = c4Subset[!(rownames(c4Subset) %in% rownames(hc4)),]
lc5 = c5Subset[!(rownames(c5Subset) %in% rownames(hc5)),]

highconf = rbind(hc1,hc2,hc3,hc4,hc5)
lowconf = rbind(lc1,lc2,lc3,lc4,lc5)

highconfDF = fullDF[rownames(highconf), fullFeatureOrder]
heatmapDF = highconfDF[,rev(colnames(highconfDF))]
heatmapDFMatrix = as.matrix(heatmapDF)
#image(heatmapDFMatrix, useRaster=TRUE, col=gray.colors(20,start=1,end=0,gamma = 2.2, alpha=NULL))
melted = melt(heatmapDFMatrix)
ggplot(data = melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") +
  scale_x_discrete(expand=c(0,0), breaks=c(rownames(hc1)[[1]]
                                           , rownames(hc2)[[1]]
                                           , rownames(hc3)[[1]]
                                           , rownames(hc4)[[1]]
                                           , rownames(hc5)[[1]])
                   , labels = c("C1","C2","C3","C4","C5"))

lowconfDF = fullDF[rownames(lowconf), fullFeatureOrder]
heatmapDF = lowconfDF[,rev(colnames(lowconfDF))]
heatmapDFMatrix = as.matrix(heatmapDF)
#image(heatmapDFMatrix, useRaster=TRUE, col=gray.colors(20,start=1,end=0,gamma = 2.2, alpha=NULL))
melted = melt(heatmapDFMatrix)
ggplot(data = melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") +
  scale_x_discrete(expand=c(0,0), breaks=c(rownames(lc1)[[1]]
                                           , rownames(lc2)[[1]]
                                           , rownames(lc3)[[1]]
                                           , rownames(lc4)[[1]]
                                           , rownames(lc5)[[1]])
                   , labels = c("C1","C2","C3","C4","C5"))
