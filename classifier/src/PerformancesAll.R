rm(list = ls())
library(ggplot2)
nmfFile = "WrittenFiles/NMF_Dist.txt"
annFile = "WrittenFiles/NN_Dist.txt"
rfFile = "WrittenFiles/RF_Dist.txt"
nmfFile_noCNA = "WrittenFiles/NMF_noCNA_Dist.txt"
annFile_noCNA = "WrittenFiles/NN_noCNA_Dist.txt"
rfFile_noCNA = "WrittenFiles/RF_noCNA_Dist.txt"
nmfFile_noCNASV = "WrittenFiles/NMF_noSV_noCNA_Dist.txt"
annFile_noCNASV = "WrittenFiles/NN_noSV_noCNA_Dist.txt"
rfFile_noCNASV = "WrittenFiles/RF_noSV_noCNA_Dist.txt"

nmfTable = read.csv(nmfFile, header = FALSE, sep = ',')
annTable = read.csv(annFile, header = FALSE, sep = ',')
rfTable = read.csv(rfFile, header = FALSE, sep = ',')
nmfTable_noCNA = read.csv(nmfFile_noCNA, header = FALSE, sep = ',')
annTable_noCNA = read.csv(annFile_noCNA, header = FALSE, sep = ',')
rfTable_noCNA = read.csv(rfFile_noCNA, header = FALSE, sep = ',')
nmfTable_noCNASV = read.csv(nmfFile_noCNASV, header = FALSE, sep = ',')
annTable_noCNASV = read.csv(annFile_noCNASV, header = FALSE, sep = ',')
rfTable_noCNASV = read.csv(rfFile_noCNASV, header = FALSE, sep = ',')

colnames(nmfTable) = c("weightedPearson", "accuracy", "performance")
colnames(annTable) = c("weightedPearson", "accuracy", "performance")
colnames(rfTable) = c("weightedPearson", "accuracy", "performance")
colnames(nmfTable_noCNA) = c("weightedPearson", "accuracy", "performance")
colnames(annTable_noCNA) = c("weightedPearson", "accuracy", "performance")
colnames(rfTable_noCNA) = c("weightedPearson", "accuracy", "performance")
colnames(nmfTable_noCNASV) = c("weightedPearson", "accuracy", "performance")
colnames(annTable_noCNASV) = c("weightedPearson", "accuracy", "performance")
colnames(rfTable_noCNASV) = c("weightedPearson", "accuracy", "performance")

nmfTable = nmfTable[order(nmfTable$performance),]
annTable = annTable[order(annTable$performance),]
rfTable = rfTable[order(rfTable$performance),]
nmfTable_noCNA = nmfTable_noCNA[order(nmfTable_noCNA$performance),]
annTable_noCNA = annTable_noCNA[order(annTable_noCNA$performance),]
rfTable_noCNA = rfTable_noCNA[order(rfTable_noCNA$performance),]
nmfTable_noCNASV = nmfTable_noCNASV[order(nmfTable_noCNASV$performance),]
annTable_noCNASV = annTable_noCNASV[order(annTable_noCNASV$performance),]
rfTable_noCNASV = rfTable_noCNASV[order(rfTable_noCNASV$performance),]

confidence = 0.84
cutoff = nrow(nmfTable)*(1/2*(1-confidence))+1

medianNMF = median(nmfTable$performance)
medianANN = median(annTable$performance)
medianRF = median(rfTable$performance)
medianNMF_noCNA = median(nmfTable_noCNA$performance)
medianANN_noCNA = median(annTable_noCNA$performance)
medianRF_noCNA = median(rfTable_noCNA$performance)
medianNMF_noCNASV = median(nmfTable_noCNASV$performance)
medianANN_noCNASV = median(annTable_noCNASV$performance)
medianRF_noCNASV = median(rfTable_noCNASV$performance)


nmfTable = nmfTable[cutoff:(nrow(nmfTable)-cutoff),]
annTable = annTable[cutoff:(nrow(annTable)-cutoff),]
rfTable = rfTable[cutoff:(nrow(rfTable)-cutoff),]
nmfTable_noCNA = nmfTable_noCNA[cutoff:(nrow(nmfTable_noCNA)-cutoff),]
annTable_noCNA = annTable_noCNA[cutoff:(nrow(annTable_noCNA)-cutoff),]
rfTable_noCNA = rfTable_noCNA[cutoff:(nrow(rfTable_noCNA)-cutoff),]
nmfTable_noCNASV = nmfTable_noCNASV[cutoff:(nrow(nmfTable_noCNASV)-cutoff),]
annTable_noCNASV = annTable_noCNASV[cutoff:(nrow(annTable_noCNASV)-cutoff),]
rfTable_noCNASV = rfTable_noCNASV[cutoff:(nrow(rfTable_noCNASV)-cutoff),]

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

uppers = c(maxNMF_noCNASV, maxRF_noCNASV, maxANN_noCNASV, 
           maxNMF_noCNA, maxRF_noCNA, maxANN_noCNA,
           maxNMF, maxRF, maxANN)

lowers = c(minNMF_noCNASV, minRF_noCNASV, minANN_noCNASV,
           minNMF_noCNA, minRF_noCNA, minANN_noCNA,
           minNMF, minRF, minANN)

medians = c(medianNMF_noCNASV, medianRF_noCNASV, medianANN_noCNASV,
            medianNMF_noCNA, medianRF_noCNA, medianANN_noCNA,
            medianNMF, medianRF, medianANN)

uppers = uppers-medians
lowers = medians-lowers

data = data.frame(cbind(medians, c("NMF_noCNASV","RF_noCNASV","ANN_noCNASV",
                                   "NMF_noCNA", "RF_noCNA", "ANN_noCNA",
                                   "NMF","RF","ANN")))

colnames(data)[2] = "Model"
data$medians = as.numeric(as.character(data$medians))

p = ggplot(data=data, aes(x=Model, y=medians)) + scale_y_continuous(limit = c(-.2,1), breaks=seq(-.2,1,by=.1)) + 
  scale_x_discrete(limits = c("NMF_noCNASV","RF_noCNASV","ANN_noCNASV",
                              "NMF_noCNA","RF_noCNA","ANN_noCNA",
                              "NMF","RF","ANN")) +
  geom_point() +
  geom_errorbar(aes(ymax = medians+uppers, ymin = medians-lowers), width=.1) +
  ggtitle("Performance (weighted pearson * accuracy) Over 200 Iterations") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Model") +
  ylab("Performance") +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=45,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=25,angle=0,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))
print(p)
fn = "Plots/PerformancePlot.jpeg"
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()
