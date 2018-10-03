rm(list = ls())
library(ggplot2)
nmfFile = "WrittenFiles/NMF_Dist.txt"
annFile = "WrittenFiles/NN_Dist.txt"
rfFile = "WrittenFiles/RF_Dist.txt"
nmfFile_noCNA = "WrittenFiles/NMF_Dist_noCNA.txt"
annFile_noCNA = "WrittenFiles/NN_Dist_noCNA.txt"
rfFile_noCNA = "WrittenFiles/RF_Dist_noCNA.txt"
nmfFile_noCNASV = "WrittenFiles/NMF_Dist_noSV_noCNA.txt"
annFile_noCNASV = "WrittenFiles/NN_Dist_noSV_noCNA.txt"
rfFile_noCNASV = "WrittenFiles/RF_Dist_noSV_noCNA.txt"
nmfFile_GD = "WrittenFiles/NMF_Dist_GD.txt"
annFile_GD = "WrittenFiles/NN_Dist_GD.txt"
rfFile_GD = "WrittenFiles/RF_Dist_GD.txt"
annFile_BinaryLabels = "WrittenFiles/NN_Dist_BinaryLabels.txt"
rfFile_BinaryLabels = "WrittenFiles/RF_Dist_BinaryLabels.txt"
nmfFile_FixedLabels = "WrittenFiles/NMF_Dist_FixedSets.txt"
annFile_FixedLabels = "WrittenFiles/NN_Dist_FixedSets.txt"
rfFile_FixedLabels = "WrittenFiles/RF_Dist_FixedSets.txt"
nmfFile_noSV = "WrittenFiles/NMF_Dist_noSV.txt"
annFile_noSV = "WrittenFiles/NN_Dist_noSV.txt"
rfFile_noSV = "WrittenFiles/RF_Dist_noSV.txt"
nmfFile_full = "WrittenFiles/NMF_Dist_161F.txt"
annFile_full = "WrittenFiles/NN_Dist_161F.txt"
rfFile_full = "WrittenFiles/RF_Dist_161F.txt"

nmfTable = read.csv(nmfFile, header = FALSE, sep = ',')
annTable = read.csv(annFile, header = FALSE, sep = ',')
rfTable = read.csv(rfFile, header = FALSE, sep = ',')
nmfTable_noCNA = read.csv(nmfFile_noCNA, header = FALSE, sep = ',')
annTable_noCNA = read.csv(annFile_noCNA, header = FALSE, sep = ',')
rfTable_noCNA = read.csv(rfFile_noCNA, header = FALSE, sep = ',')
nmfTable_noCNASV = read.csv(nmfFile_noCNASV, header = FALSE, sep = ',')
annTable_noCNASV = read.csv(annFile_noCNASV, header = FALSE, sep = ',')
rfTable_noCNASV = read.csv(rfFile_noCNASV, header = FALSE, sep = ',')
nmfTable_GD = read.csv(nmfFile_GD, header = FALSE, sep = ",")
annTable_GD = read.csv(annFile_GD, header = FALSE, sep = ",")
rfTable_GD = read.csv(rfFile_GD, header = FALSE, sep = ",")
annTable_BinaryLabels = read.csv(annFile_BinaryLabels, header=FALSE, sep=",")
rfTable_BinaryLabels = read.csv(rfFile_BinaryLabels, header=FALSE, sep=",")
nmfTable_FixedLabels = read.csv(nmfFile_FixedLabels, header=FALSE, sep=",")
annTable_FixedLabels = read.csv(annFile_FixedLabels, header=FALSE, sep=",")
rfTable_FixedLabels = read.csv(rfFile_FixedLabels, header=FALSE, sep=",")
nmfTable_noSV = read.csv(nmfFile_noSV, header=FALSE, sep=",")
annTable_noSV = read.csv(annFile_noSV, header=FALSE, sep=",")
rfTable_noSV = read.csv(rfFile_noSV, header=FALSE, sep=",")
nmfTable_full = read.csv(nmfFile_full, header=FALSE, sep=",")
annTable_full = read.csv(annFile_full, header=FALSE, sep=",")
rfTable_full = read.csv(rfFile_full, header=FALSE, sep=",")

colnames(nmfTable) = c("weightedPearson", "accuracy", "performance")
colnames(annTable) = c("weightedPearson", "accuracy", "performance")
colnames(rfTable) = c("weightedPearson", "accuracy", "performance")
colnames(nmfTable_noCNA) = c("weightedPearson", "accuracy", "performance")
colnames(annTable_noCNA) = c("weightedPearson", "accuracy", "performance")
colnames(rfTable_noCNA) = c("weightedPearson", "accuracy", "performance")
colnames(nmfTable_noCNASV) = c("weightedPearson", "accuracy", "performance")
colnames(annTable_noCNASV) = c("weightedPearson", "accuracy", "performance")
colnames(rfTable_noCNASV) = c("weightedPearson", "accuracy", "performance")
colnames(nmfTable_GD) = c("weightedPearson", "accuracy", "performance")
colnames(annTable_GD) = c("weightedPearson", "accuracy", "performance")
colnames(rfTable_GD) = c("weightedPearson", "accuracy", "performance")
colnames(rfTable_BinaryLabels) = c("weightedPearson", "accuracy", "performance")
colnames(annTable_BinaryLabels) = c("weightedPearson", "accuracy", "performance")
colnames(rfTable_FixedLabels) = c("weightedPearson", "accuracy", "performance")
colnames(annTable_FixedLabels) = c("weightedPearson", "accuracy", "performance")
colnames(nmfTable_FixedLabels) = c("weightedPearson", "accuracy", "performance")
colnames(nmfTable_noSV) = c("weightedPearson", "accuracy", "performance")
colnames(annTable_noSV) = c("weightedPearson", "accuracy", "performance")
colnames(rfTable_noSV) = c("weightedPearson", "accuracy", "performance")
colnames(nmfTable_full) = c("weightedPearson", "accuracy", "performance")
colnames(annTable_full) = c("weightedPearson", "accuracy", "performance")
colnames(rfTable_full) = c("weightedPearson", "accuracy", "performance")


nmfTable = nmfTable[order(nmfTable$performance),]
annTable = annTable[order(annTable$performance),]
rfTable = rfTable[order(rfTable$performance),]
nmfTable_noCNA = nmfTable_noCNA[order(nmfTable_noCNA$performance),]
annTable_noCNA = annTable_noCNA[order(annTable_noCNA$performance),]
rfTable_noCNA = rfTable_noCNA[order(rfTable_noCNA$performance),]
nmfTable_noCNASV = nmfTable_noCNASV[order(nmfTable_noCNASV$performance),]
annTable_noCNASV = annTable_noCNASV[order(annTable_noCNASV$performance),]
rfTable_noCNASV = rfTable_noCNASV[order(rfTable_noCNASV$performance),]
nmfTable_GD = nmfTable_GD[order(nmfTable_GD$performance),]
annTable_GD = annTable_GD[order(annTable_GD$performance),]
rfTable_GD = rfTable_GD[order(rfTable_GD$performance),]
rfTable_BinaryLabels = rfTable_BinaryLabels[order(rfTable_BinaryLabels$performance),]
annTable_BinaryLabels = annTable_BinaryLabels[order(annTable_BinaryLabels$performance),]
nmfTable_FixedLabels = nmfTable_FixedLabels[order(nmfTable_FixedLabels$performance),]
annTable_FixedLabels = annTable_FixedLabels[order(annTable_FixedLabels$performance),]
rfTable_FixedLabels = rfTable_FixedLabels[order(rfTable_FixedLabels$performance),]
nmfTable_noSV = nmfTable_noSV[order(nmfTable_noSV$performance),]
annTable_noSV = annTable_noSV[order(annTable_noSV$performance),]
rfTable_full = rfTable_noSV[order(rfTable_noSV$performance),]
nmfTable_full = nmfTable_full[order(nmfTable_full$performance),]
annTable_full = annTable_full[order(annTable_full$performance),]
rfTable_full = rfTable_full[order(rfTable_full$performance),]

confidence = 0.84

medianNMF = median(nmfTable$performance)
medianANN = median(annTable$performance)
medianRF = median(rfTable$performance)
medianNMF_noCNA = median(nmfTable_noCNA$performance)
medianANN_noCNA = median(annTable_noCNA$performance)
medianRF_noCNA = median(rfTable_noCNA$performance)
medianNMF_noCNASV = median(nmfTable_noCNASV$performance)
medianANN_noCNASV = median(annTable_noCNASV$performance)
medianRF_noCNASV = median(rfTable_noCNASV$performance)
medianNMF_GD = median(nmfTable_GD$performance)
medianANN_GD = median(annTable_GD$performance)
medianRF_GD = median(rfTable_GD$performance)
medianANN_BinaryLabels = median(annTable_BinaryLabels$performance)
medianRF_BinaryLabels = median(rfTable_BinaryLabels$performance)
medianNMF_FixedLabels = median(nmfTable_FixedLabels$performance)
medianANN_FixedLabels = median(annTable_FixedLabels$performance)
medianRF_FixedLabels = median(rfTable_FixedLabels$performance)
medianNMF_noSV = median(nmfTable_noSV$performance)
medianANN_noSV = median(annTable_noSV$performance)
medianRF_noSV = median(rfTable_noSV$performance)
medianNMF_full = median(nmfTable_full$performance)
medianANN_full = median(annTable_full$performance)
medianRF_full = median(rfTable_full$performance)

cutoff = nrow(nmfTable)*(1/2*(1-confidence))+1
nmfTable = nmfTable[cutoff:(nrow(nmfTable)-cutoff),]
cutoff = nrow(annTable)*(1/2*(1-confidence))+1
annTable = annTable[cutoff:(nrow(rfTable)-cutoff),]
cutoff = nrow(rfTable)*(1/2*(1-confidence))+1
rfTable = rfTable[cutoff:(nrow(rfTable)-cutoff),]

cutoff = nrow(nmfTable_noCNA)*(1/2*(1-confidence))+1
nmfTable_noCNA = nmfTable_noCNA[cutoff:(nrow(nmfTable_noCNA)-cutoff),]
cutoff = nrow(annTable_noCNA)*(1/2*(1-confidence))+1
annTable_noCNA = annTable_noCNA[cutoff:(nrow(annTable_noCNA)-cutoff),]
cutoff = nrow(rfTable_noCNA)*(1/2*(1-confidence))+1
rfTable_noCNA = rfTable_noCNA[cutoff:(nrow(rfTable_noCNA)-cutoff),]

cutoff = nrow(nmfTable_noCNASV)*(1/2*(1-confidence))+1
nmfTable_noCNASV = nmfTable_noCNASV[cutoff:(nrow(nmfTable_noCNASV)-cutoff),]
cutoff = nrow(annTable_noCNASV)*(1/2*(1-confidence))+1
annTable_noCNASV = annTable_noCNASV[cutoff:(nrow(annTable_noCNASV)-cutoff),]
cutoff = nrow(rfTable_noCNASV)*(1/2*(1-confidence))+1
rfTable_noCNASV = rfTable_noCNASV[cutoff:(nrow(rfTable_noCNASV)-cutoff),]

cutoff = nrow(nmfTable_GD)*(1/2*(1-confidence))+1
nmfTable_GD = nmfTable_GD[cutoff:(nrow(nmfTable_GD)-cutoff),]
cutoff = nrow(annTable_GD)*(1/2*(1-confidence))+1
annTable_GD = annTable_GD[cutoff:(nrow(annTable_GD)-cutoff),]
cutoff = nrow(rfTable_GD)*(1/2*(1-confidence))+1
rfTable_GD = rfTable_GD[cutoff:(nrow(rfTable_GD)-cutoff),]

cutoff = nrow(annTable_BinaryLabels)*(1/2*(1-confidence))+1
annTable_BinaryLabels = annTable_BinaryLabels[cutoff:(nrow(annTable_BinaryLabels)-cutoff),]
cutoff = nrow(rfTable_BinaryLabels)*(1/2*(1-confidence))+1
rfTable_BinaryLabels = rfTable_BinaryLabels[cutoff:(nrow(rfTable_BinaryLabels)-cutoff),]

cutoff = nrow(nmfTable_FixedLabels)*(1/2*(1-confidence))+1
nmfTable_FixedLabels = nmfTable_FixedLabels[cutoff:(nrow(nmfTable_FixedLabels)-cutoff),]
cutoff = nrow(annTable_FixedLabels)*(1/2*(1-confidence))+1
annTable_FixedLabels = annTable_FixedLabels[cutoff:(nrow(annTable_FixedLabels)-cutoff),]
cutoff = nrow(rfTable_FixedLabels)*(1/2*(1-confidence))+1
rfTable_FixedLabels = rfTable_FixedLabels[cutoff:(nrow(rfTable_FixedLabels)-cutoff),]

cutoff = nrow(nmfTable_noSV)*(1/2*(1-confidence))+1
nmfTable_noSV = nmfTable_noSV[cutoff:(nrow(nmfTable_noSV)-cutoff),]
cutoff = nrow(annTable_noSV)*(1/2*(1-confidence))+1
annTable_noSV = annTable_noSV[cutoff:(nrow(annTable_noSV)-cutoff),]
cutoff = nrow(rfTable_noSV)*(1/2*(1-confidence))+1
rfTable_noSV = rfTable_noSV[cutoff:(nrow(rfTable_noSV)-cutoff),]

cutoff = nrow(nmfTable_full)*(1/2*(1-confidence))+1
nmfTable_full = nmfTable_full[cutoff:(nrow(nmfTable_full)-cutoff),]
cutoff = nrow(annTable_full)*(1/2*(1-confidence))+1
annTable_full = annTable_full[cutoff:(nrow(annTable_full)-cutoff),]
cutoff = nrow(rfTable_full)*(1/2*(1-confidence))+1
rfTable_full = rfTable_full[cutoff:(nrow(rfTable_full)-cutoff),]

maxNMF = max(nmfTable$performance)
minNMF = min(nmfTable$performance)
maxANN = max(annTable$performance)
minANN = min(annTable$performance)
maxRF = max(rfTable$performance)
minRF = min(rfTable$performance)

maxNMF_noCNA = max(nmfTable_noCNA$performance)
minNMF_noCNA = min(nmfTable_noCNA$performance)
maxANN_noCNA = max(annTable_noCNA$performance)
minANN_noCNA = min(annTable_noCNA$performance)
maxRF_noCNA = max(rfTable_noCNA$performance)
minRF_noCNA = min(rfTable_noCNA$performance)

maxNMF_noCNASV = max(nmfTable_noCNASV$performance)
minNMF_noCNASV = min(nmfTable_noCNASV$performance)
maxANN_noCNASV = max(annTable_noCNASV$performance)
minANN_noCNASV = min(annTable_noCNASV$performance)
maxRF_noCNASV = max(rfTable_noCNASV$performance)
minRF_noCNASV = min(rfTable_noCNASV$performance)

maxNMF_GD = max(nmfTable_GD$performance)
minNMF_GD = min(nmfTable_GD$performance)
maxANN_GD = max(annTable_GD$performance)
minANN_GD = min(annTable_GD$performance)
maxRF_GD = max(rfTable_GD$performance)
minRF_GD = min(rfTable_GD$performance)

maxANN_BinaryLabels = max(annTable_BinaryLabels$performance)
minANN_BinaryLabels = min(annTable_BinaryLabels$performance)
maxRF_BinaryLabels = max(rfTable_BinaryLabels$performance)
minRF_BinaryLabels = min(rfTable_BinaryLabels$performance)

maxNMF_FixedLabels = max(nmfTable_FixedLabels$performance)
minNMF_FixedLabels = min(nmfTable_FixedLabels$performance)
maxANN_FixedLabels = max(annTable_FixedLabels$performance)
minANN_FixedLabels = min(annTable_FixedLabels$performance)
maxRF_FixedLabels = max(rfTable_FixedLabels$performance)
minRF_FixedLabels = min(rfTable_FixedLabels$performance)

maxNMF_noSV = max(nmfTable_noSV$performance)
minNMF_noSV = min(nmfTable_noSV$performance)
maxANN_noSV = max(annTable_noSV$performance)
minANN_noSV = min(annTable_noSV$performance)
maxRF_noSV = max(rfTable_noSV$performance)
minRF_noSV = min(rfTable_noSV$performance)

maxNMF_full = max(nmfTable_full$performance)
minNMF_full = min(nmfTable_full$performance)
maxANN_full = max(annTable_full$performance)
minANN_full = min(annTable_full$performance)
maxRF_full = max(rfTable_full$performance)
minRF_full = min(rfTable_full$performance)




uppers = c(maxNMF_noCNASV, maxRF_noCNASV, maxANN_noCNASV, 
           maxNMF_noCNA, maxRF_noCNA, maxANN_noCNA,
           maxNMF_GD, maxRF_GD, maxANN_GD,
           maxNMF, maxRF, maxANN)

uppers_Default = c(maxNMF, maxRF, maxANN)
uppers_GD = c(maxNMF_GD, maxRF_GD, maxANN_GD)
uppers_noCNA = c(maxNMF_noCNA, maxRF_noCNA, maxANN_noCNA)
uppers_noCNASV = c(maxNMF_noCNASV, maxRF_noCNASV, maxANN_noCNASV)
uppers_BinaryLabels = c(maxRF_BinaryLabels, maxANN_BinaryLabels)
uppers_FixedLabels = c(maxNMF_FixedLabels, maxRF_FixedLabels, maxANN_FixedLabels)
uppers_noSV = c(maxNMF_noSV, maxRF_noSV, maxANN_noSV)
uppers_full = c(maxNMF_full, maxRF_full, maxANN_full)

lowers = c(minNMF_noCNASV, minRF_noCNASV, minANN_noCNASV,
           minNMF_noCNA, minRF_noCNA, minANN_noCNA,
           minNMF_GD, minRF_GD, minANN_GD,
           minNMF, minRF, minANN)

lowers_Default = c(minNMF, minRF, minANN)
lowers_GD = c(minNMF_GD, minRF_GD, minANN_GD)
lowers_noCNA = c(minNMF_noCNA, minRF_noCNA, minANN_noCNA)
lowers_noCNASV = c(minNMF_noCNASV, minRF_noCNASV, minANN_noCNASV)
lowers_BinaryLabels = c(minRF_BinaryLabels, minANN_BinaryLabels)
lowers_FixedLabels = c(minNMF_FixedLabels, minRF_FixedLabels, minANN_FixedLabels)
lowers_noSV = c(minNMF_noSV, minRF_noSV, minANN_noSV)
lowers_full = c(minNMF_full, minRF_full, minANN_full)

medians_Default = c(medianNMF, medianRF, medianANN)
medians_noCNASV = c(medianNMF_noCNASV, medianRF_noCNASV, medianANN_noCNASV)
medians_noCNA = c(medianNMF_noCNA, medianRF_noCNA, medianANN_noCNA)
medians_GD = c(medianNMF_GD, medianRF_GD, medianANN_GD)
medians_BinaryLabels = c(medianRF_BinaryLabels, medianANN_BinaryLabels)
medians_FixedLabels = c(medianNMF_FixedLabels, medianRF_FixedLabels, medianANN_FixedLabels)
medians_noSV = c(medianNMF_noSV, medianRF_noSV, medianANN_noSV)
medians_full = c(medianNMF_full, medianRF_full, medianANN_full)


uppers_Default = uppers_Default - medians_Default
uppers_GD = uppers_GD - medians_GD
uppers_noCNA = uppers_noCNA - medians_noCNA
uppers_noCNASV = uppers_noCNASV - medians_noCNASV
uppers_BinaryLabels = uppers_BinaryLabels - medians_BinaryLabels
uppers_FixedLabels = uppers_FixedLabels - medians_FixedLabels
uppers_noSV = uppers_noSV - medians_noSV
uppers_full = uppers_full - medians_full

lowers_Default = medians_Default - lowers_Default
lowers_GD = medians_GD - lowers_GD
lowers_noCNA = medians_noCNA - lowers_noCNA
lowers_noCNASV = medians_noCNASV - lowers_noCNASV
lowers_BinaryLabels = medians_BinaryLabels - lowers_BinaryLabels
lowers_FixedLabels = medians_FixedLabels - lowers_FixedLabels
lowers_noSV = medians_noSV - lowers_noSV
lowers_full = medians_full - lowers_full

data_Default = data.frame(cbind(medians_Default, c("NMF","RF","ANN")))
data_GD = data.frame(cbind(medians_GD, c("NMF_GD", "RF_GD", "ANN_GD")))
data_noCNA = data.frame(cbind(medians_noCNA, c("NMF_noCNA","RF_noCNA","ANN_noCNA")))
data_noCNASV = data.frame(cbind(medians_noCNASV, c("NMF_noCNASV","RF_noCNASV","ANN_noCNASV")))
data_BinaryLabels = data.frame(cbind(medians_BinaryLabels, c("RF_BinaryLabels", "ANN_BinaryLabels")))
data_FixedLabels = data.frame(cbind(medians_FixedLabels, c("NMF_FixedLabels", "RF_FixedLabels", "ANN_FixedLabels")))
data_noSV = data.frame(cbind(medians_noSV, c("NMF_noSV", "RF_noSV", "ANN_noSV")))
data_full = data.frame(cbind(medians_full, c("NMF_161", "RF_161" , "ANN_161")))

colnames(data_Default)[2] = "Model"
colnames(data_GD)[2] = "Model"
colnames(data_noCNA)[2] = "Model"
colnames(data_noCNASV)[2] = "Model"
colnames(data_BinaryLabels)[2] = "Model"
colnames(data_FixedLabels)[2] = "Model"
colnames(data_noSV)[2] = "Model"
colnames(data_full)[2] = "Model"

data_Default$medians_Default = as.numeric(as.character(data_Default$medians))
data_GD$medians_GD = as.numeric(as.character(data_GD$medians))
data_noCNA$medians_noCNA = as.numeric(as.character(data_noCNA$medians))
data_noCNASV$medians_noCNASV = as.numeric(as.character(data_noCNASV$medians))
data_BinaryLabels$medians_BinaryLabels = as.numeric(as.character(data_BinaryLabels$medians))
data_FixedLabels$medians_FixedLabels = as.numeric(as.character(data_FixedLabels$medians))
data_noSV$medians_noSV = as.numeric(as.character(data_noSV$medians))
data_full$medians_full = as.numeric(as.character(data_full$medians))

p = ggplot(data=data_Default, aes(x=Model, y=medians_Default)) + scale_y_continuous(limit = c(0,1), breaks=seq(0,1,by=.1)) + 
  scale_x_discrete(limits = c("NMF","RF","ANN")) +
  geom_point() +
  geom_errorbar(aes(ymax = medians_Default+uppers_Default, ymin = medians_Default-lowers_Default), width=.1) +
  ggtitle("Default: Performance Over 100 Iterations") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Model") +
  ylab("Performance") +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))
print(p)
fn = "Plots/PerformancePlot_Default.jpeg"
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()

p = ggplot(data=data_BinaryLabels, aes(x=Model, y=medians_BinaryLabels)) + scale_y_continuous(limit = c(0,1), breaks=seq(0,1,by=.1)) + 
  scale_x_discrete(limits = c("RF_BinaryLabels","ANN_BinaryLabels")) +
  geom_point() +
  geom_errorbar(aes(ymax = medians_BinaryLabels+uppers_BinaryLabels, ymin = medians_BinaryLabels-lowers_BinaryLabels), width=.1) +
  ggtitle("Binary Labels: Performance (weighted pearson * accuracy) Over 100 Iterations") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Model") +
  ylab("Performance") +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))
print(p)
fn = "Plots/PerformancePlot_BinaryLabels.jpeg"
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()

p = ggplot(data=data_GD, aes(x=Model, y=medians_GD)) + scale_y_continuous(limit = c(0,1), breaks=seq(0,1,by=.1)) + 
  scale_x_discrete(limits = c("NMF_GD","RF_GD", "ANN_GD")) +
  geom_point() +
  geom_errorbar(aes(ymax = medians_GD+uppers_GD, ymin = medians_GD-lowers_GD), width=.1) +
  ggtitle("With Genome Doubling: Performance Over 100 Iterations") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Model") +
  ylab("Performance") +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))
print(p)
fn = "Plots/PerformancePlot_GD.jpeg"
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()

p = ggplot(data=data_noCNA, aes(x=Model, y=medians_noCNA)) + scale_y_continuous(limit = c(0,1), breaks=seq(0,1,by=.1)) + 
  scale_x_discrete(limits = c("NMF_noCNA","RF_noCNA", "ANN_noCNA")) +
  geom_point() +
  geom_errorbar(aes(ymax = medians_noCNA+uppers_noCNA, ymin = medians_noCNA-lowers_noCNA), width=.1) +
  ggtitle("no CNA: Performance Over 100 Iterations") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Model") +
  ylab("Performance") +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))
print(p)
fn = "Plots/PerformancePlot_noCNA.jpeg"
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()

p = ggplot(data=data_noCNASV, aes(x=Model, y=medians_noCNASV)) + scale_y_continuous(limit = c(0,1), breaks=seq(0,1,by=.1)) + 
  scale_x_discrete(limits = c("NMF_noCNASV","RF_noCNASV", "ANN_noCNASV")) +
  geom_point() +
  geom_errorbar(aes(ymax = medians_noCNASV+uppers_noCNASV, ymin = medians_noCNASV-lowers_noCNASV), width=.1) +
  ggtitle("no CNA and no SV: Performance Over 100 Iterations") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Model") +
  ylab("Performance") +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))
print(p)
fn = "Plots/PerformancePlot_noCNASV.jpeg"
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()

##combined GD and default

combinedDF = data_GD
colnames(combinedDF) = c("Medians", "Model")
toAdd = data_Default
colnames(toAdd) = c("Medians", "Model")
combinedDF = rbind(combinedDF, toAdd)
combinedUppers = c(uppers_GD, uppers_Default)
combinedLowers = c(lowers_GD, lowers_Default)
combinedDF$Experiment = c("Genome Doubling","Genome Doubling","Genome Doubling",
                          "Default","Default","Default")

p = ggplot(data=combinedDF, aes(x=Model, y=Medians)) + scale_y_continuous(limit = c(0,1), breaks=seq(0,1,by=.1)) + 
  scale_x_discrete(limits = c("NMF_GD", "RF_GD", "ANN_GD", "NMF", "RF", "ANN")) +
  geom_point(aes(color=Experiment)) +
  geom_errorbar(aes(ymax = Medians+combinedUppers, ymin = Medians-combinedLowers, color=Experiment), width=.1) +
  ggtitle("GD and Default: Performance Over 100 Iterations") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Model") +
  ylab("Performance") +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))
print(p)

fn = "Plots/PerformancePlot_GD_Default.jpeg"
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()

combinedDF_Fixed = data_FixedLabels
colnames(combinedDF_Fixed) = c("Medians", "Model")
toAdd = data_Default
colnames(toAdd) = c("Medians", "Model")
combinedDF_Fixed = rbind(combinedDF_Fixed, toAdd)
combinedUppers = c(uppers_FixedLabels, uppers_Default)
combinedLowers = c(lowers_FixedLabels, lowers_Default)
combinedDF_Fixed$Experiment = c("Fixed Validation","Fixed Validation","Fixed Validation",
                                "Default","Default","Default")

p = ggplot(data=combinedDF_Fixed, aes(x=Model, y=Medians)) + scale_y_continuous(limit = c(0,1), breaks=seq(0,1,by=.1)) + 
  scale_x_discrete(limits = c("NMF_FixedLabels", "RF_FixedLabels", "ANN_FixedLabels", "NMF", "RF", "ANN")) +
  geom_point(aes(color=Experiment)) +
  geom_errorbar(aes(ymax = Medians+combinedUppers, ymin = Medians-combinedLowers, color=Experiment), width=.1) +
  ggtitle("Fixed Validation Sets and Default: Performance Over 100 Iterations") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Model") +
  ylab("Performance") +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))


fn = "Plots/PerformancePlot_FixedValidation_Default.jpeg"
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()


combinedDF_noSV = data_noSV
colnames(combinedDF_noSV) = c("Medians", "Model")
toAdd = data_Default
colnames(toAdd) = c("Medians", "Model")
combinedDF_noSV = rbind(combinedDF_noSV, toAdd)
combinedUppers = c(uppers_FixedLabels, uppers_Default)
combinedLowers = c(lowers_FixedLabels, lowers_Default)
combinedDF_noSV$Experiment = c("No SV","No SV","No SV",
                               "Default","Default","Default")

p = ggplot(data=combinedDF_noSV, aes(x=Model, y=Medians)) + scale_y_continuous(limit = c(0,1), breaks=seq(0,1,by=.1)) + 
  scale_x_discrete(limits = c("NMF_noSV", "RF_noSV", "ANN_noSV", "NMF", "RF", "ANN")) +
  geom_point(aes(color=Experiment)) +
  geom_errorbar(aes(ymax = Medians+combinedUppers, ymin = Medians-combinedLowers, color=Experiment), width=.1) +
  ggtitle("No SV and Default: Performance Over 100 Iterations") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Model") +
  ylab("Performance") +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))


fn = "Plots/PerformancePlot_noSV_Default.jpeg"
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()

combinedDF_full = data_full
colnames(combinedDF_full) = c("Medians", "Model")
toAdd = data_Default
colnames(toAdd) = c("Medians", "Model")
combinedDF_full = rbind(combinedDF_full, toAdd)
combinedUppers = c(uppers_FixedLabels, uppers_Default)
combinedLowers = c(lowers_FixedLabels, lowers_Default)
combinedDF_full$Experiment = c("161 Features","161 Features","161 Features",
                               "Default","Default","Default")

p = ggplot(data=combinedDF_full, aes(x=Model, y=Medians)) + scale_y_continuous(limit = c(0,1), breaks=seq(0,1,by=.1)) + 
  scale_x_discrete(limits = c("NMF_161", "RF_161", "ANN_161", "NMF", "RF", "ANN")) +
  geom_point(aes(color=Experiment)) +
  geom_errorbar(aes(ymax = Medians+combinedUppers, ymin = Medians-combinedLowers, color=Experiment), width=.1) +
  ggtitle("Full 161 and Default: Performance Over 100 Iterations") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Model") +
  ylab("Performance") +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))


fn = "Plots/PerformancePlot_full_Default.jpeg"
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()