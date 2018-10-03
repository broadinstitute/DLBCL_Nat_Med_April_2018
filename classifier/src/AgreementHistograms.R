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

vec = c(nmfTable$accuracy, nmfTable_noCNA$accuracy, nmfTable_noCNASV$accuracy, nmfTable_GD$accuracy)
nmf_DF = data.frame("Accuracy" = vec)
nmf_DF = cbind(nmf_DF, c(rep("Default", 200), rep("noCNA", 200), rep("noCNASV", 200), rep("GenomeDoubling", 200)))
colnames(nmf_DF)[2] = "Model"

#p = ggplot(nmf_DF, aes(x=Accuracy, color=Model)) +
#  geom_histogram(fill="white", binwidth = .001)
# Overlaid histograms
p = ggplot(nmf_DF, aes(x=Accuracy, color=Model)) +
  geom_histogram(fill="white", alpha=0.01, position="identity", bins = 100) +
  ggtitle("Agreement in NMF") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20))

fn = "Plots/AgreementsNMF.jpeg"
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()

vec = c(rfTable$accuracy, 
        rfTable_noCNA$accuracy, 
        rfTable_noCNASV$accuracy, 
        rfTable_GD$accuracy, 
        rfTable_BinaryLabels$accuracy)

rf_DF = data.frame("Accuracy" = vec)
rf_DF = cbind(rf_DF, 
              c(rep("Default", 200), 
                rep("noCNA", 200), 
                rep("noCNASV", 200), 
                rep("GenomeDoubling", 200),
                rep("BinaryLabels", 200)))
colnames(rf_DF)[2] = "Model"

#p = ggplot(nmf_DF, aes(x=Accuracy, color=Model)) +
#  geom_histogram(fill="white", binwidth = .001)
# Overlaid histograms
p = ggplot(rf_DF, aes(x=Accuracy, color=Model)) +
  geom_histogram(fill="white", alpha=0.1, position="identity", bins = 100) +
  ggtitle("Agreement in RF") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20))

fn = "Plots/AgreementsRF.jpeg"
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()

vec = c(annTable$accuracy, 
        annTable_noCNA$accuracy, 
        annTable_noCNASV$accuracy, 
        annTable_GD$accuracy, 
        annTable_BinaryLabels$accuracy)

ann_DF = data.frame("Accuracy" = vec)
ann_DF = cbind(ann_DF, 
              c(rep("Default", 200), 
                rep("noCNA", 200), 
                rep("noCNASV", 200), 
                rep("GenomeDoubling", 200),
                rep("BinaryLabels", 200)))
colnames(ann_DF)[2] = "Model"

#p = ggplot(nmf_DF, aes(x=Accuracy, color=Model)) +
#  geom_histogram(fill="white", binwidth = .001)
# Overlaid histograms
p = ggplot(ann_DF, aes(x=Accuracy, color=Model)) +
  geom_histogram(fill="white", alpha=0.1, position="identity", bins = 100) +
  ggtitle("Agreement in ANN") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20))

fn = "Plots/AgreementsANN.jpeg"
jpeg(fn, width = 1080, height = 720)
print(p)
dev.off()