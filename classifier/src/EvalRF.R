rm(list = ls())
#set.seed(199)
knitr::opts_chunk$set(echo = TRUE)

useSV = TRUE
useCNA = TRUE
source("LoadLibraries.R")
source("src/ConstructReducedFeatureMatrix.R")
source("src/GenerateFullDF.R")
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
source("src/PearsonCustom.R")



#set.seed(1111)
#set.seed(111)

fullDF = fullDF[rownames(fullDF) %in% rownames(connectivityProbs1.0),]
baselineProbTable = csBaselines(fullDF, connectivityProbs1.0, labels, summedW)
c1Subset = baselineProbTable[baselineProbTable$cluster == 1,]
c2Subset = baselineProbTable[baselineProbTable$cluster == 2,]
c3Subset = baselineProbTable[baselineProbTable$cluster == 3,]
c4Subset = baselineProbTable[baselineProbTable$cluster == 4,]
c5Subset = baselineProbTable[baselineProbTable$cluster == 5,]
#compute max vals for histogram
maxvals = c()
for(i in 1:nrow(baselineProbTable)){
  val = max(baselineProbTable[i,2:6])
  maxvals = c(maxvals,val)
}
baselineProbTable = cbind(baselineProbTable, maxvals)
#print(paste("Number > 90%:",sum(baselineProbTable$maxvals > .90)))

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

#generate clusters from probability distribution
n = 1
nforests = 10
useProbabilities = TRUE
if(useProbabilities){
  allClusterAssignments = vector("list", length = nforests)
  for(c in 1:nforests){
    expandedClusters <- vector("list", length = nrow(trainTable)*n)
    for(i in 1:nrow(trainTable)){
      currSample = rownames(trainTable)[i]
      currProbs = baselineProbTable[rownames(baselineProbTable) == currSample, 2:6]
      clusters = sample(c(1,2,3,4,5), n, replace = TRUE, prob=currProbs)
      expandedClusters[[i]] = clusters
    }
    expandedClusters = unlist(expandedClusters)
    allClusterAssignments[[c]] = expandedClusters
  }
} else {
  allClusterAssignments = vector("list", length = nforests)
  sortedProbs = baselineProbTable[rownames(trainTable),]
  for(c in 1:nforests){
    allClusterAssignments[[c]] = sortedProbs$cluster
  }
}

#train nforests
rfKList = vector("list", length = validationSize)

for (k in 1:validationSize){
  trainTableFold = trainTable[!(rownames(trainTable) %in% validationSets[[k]]),]
  indices = match(rownames(trainTableFold), rownames(trainTable))
  rfList = vector("list", length = nforests)
  
  for(i in 1:nforests){
    if((i %% floor(nforests/10)) == 0){
      #print(i)
    }
    trainTableFold$cluster = trainTable[indices,"cluster"]
    rf = randomForest(factor(cluster) ~ ., trainTableFold, keep.forest = TRUE, norm.votes = FALSE)
    rfList[[i]] = rf
  }
  rfKList[[k]] = rfList
}

useAggregate = FALSE
allKPredictions = vector("list", length = validationSize)
if(useAggregate){
  for(k in 1:validationSize){
    validationTable = trainTable[(rownames(trainTable) %in% validationSets[[k]]),]
    allPredictions = vector("list", length = nforests)
    prevnames = NA
    for(i in 1:nforests){
      currRF = rfKList[[k]][[i]]
      predictions = predict(currRF, validationTable, predict.all = TRUE)
      currRows = names(predictions$aggregate)
      if(i > 1){
        if(!identical(currRows,prevnames)){
          print("Order of predicted different")
        }
      }
      preds = as.integer(as.character(predictions$aggregate))
      #allPredictions[[i]] = predictions$aggregate
      names(preds) = currRows
      allPredictions[[i]] = preds
      prevnames = currRows
    }
    allPredictions = as.data.frame(allPredictions)
    colnames(allPredictions) = 1:nforests
    allKPredictions[[k]] = allPredictions
  }
} else {
  for(k in 1:validationSize){
    validationTable = trainTable[(rownames(trainTable) %in% validationSets[[k]]),]
    allPredictions = vector("list", length = nforests)
    prevnames = NA
    summedProbs = NA
    for(i in 1:nforests){
      currRF = rfKList[[k]][[i]]
      predictions = predict(currRF, validationTable, predict.all = TRUE, norm.votes = TRUE, type="prob")
      currRows = names(predictions$aggregate)
      if(i > 1){
        if(!identical(currRows,prevnames)){
          print("Order of predicted different")
        }
      }
      preds = as.matrix(predictions$aggregate)
      if(i == 1){
        summedProbs = preds
      } else {
        summedProbs = summedProbs+preds
      }
      #allPredictions[[i]] = predictions$aggregate
      prevnames = currRows
    }
    summedProbs = t(apply(summedProbs, 1, function(x) x/sum(x)))
    allPredictions = as.data.frame(summedProbs)
    colnames(allPredictions) = c("P1", "P2", "P3", "P4", "P5")
    allKPredictions[[k]] = allPredictions
  }
}



kConsensusPredictions = vector("list", length = validationSize)
if(useAggregate){
  for(k in 1:validationSize){
    consensusPredictions = data.frame()
    allPredictions = allKPredictions[[k]]
    for(i in 1:nrow(allPredictions)){
      probVec = vector("integer",5)
      for(j in 1:length(allPredictions[i,])){
        val = allPredictions[i,j]
        if(val == "1"){
          probVec[1] = probVec[1] + 1
        } else if (val == "2"){
          probVec[2] = probVec[2] + 1
        } else if (val == "3"){
          probVec[3] = probVec[3] + 1
        } else if (val == "4"){
          probVec[4] = probVec[4] + 1
        } else if (val == "5"){
          probVec[5] = probVec[5] + 1
        }
      }
      probVec = probVec/sum(probVec)
      maxVal = max(probVec)
      assignedCluster = which(max(probVec) == probVec)[1]
      trueCluster = labels[rownames(labels) == rownames(allPredictions)[i],"cluster"]
      consensusPredictions = rbind(consensusPredictions, c(probVec, maxVal, assignedCluster, trueCluster))
    }
    rownames(consensusPredictions) = rownames(allPredictions)
    colnames(consensusPredictions) = c("P1","P2","P3","P4","P5","maxval","Predicted","TrueCluster")
    kConsensusPredictions[[k]] = consensusPredictions
  }
} else {
  for(k in 1:validationSize){
    consensusPredictions = data.frame()
    allPredictions = allKPredictions[[k]]
    for(i in 1:nrow(allPredictions)){
      maxval = max(allPredictions[i,])
      predictedCluster = which(maxval == allPredictions[i,])[1]
      sample = rownames(allPredictions)[i]
      trueCluster = labels[rownames(labels) == sample, "cluster"]
      consensusPredictions = rbind(consensusPredictions, 
                                   unlist(c(allPredictions[i,], maxval, predictedCluster, trueCluster)))
    }
    rownames(consensusPredictions) = rownames(allPredictions)
    consensusPredictions = cbind(rownames(consensusPredictions), consensusPredictions)
    colnames(consensusPredictions) = c("Sample","P1","P2","P3","P4","P5","maxval","Predicted","TrueCluster")
    kConsensusPredictions[[k]] = consensusPredictions
  }
}






plots = c()
if(useAggregate){
  for(k in 1:validationSize){
    consensusPredictions = kConsensusPredictions[[k]][order(kConsensusPredictions[[k]]$maxval),]
    validationProbFrame = baselineProbTable[rownames(baselineProbTable) %in% validationSets[[k]],]
    validationProbFrame = validationProbFrame[rownames(consensusPredictions),]
    
    plotDF = data.frame()
    correctness = consensusPredictions$Predicted == consensusPredictions$TrueCluster
    plotDF = cbind(consensusPredictions$maxval, validationProbFrame$maxvals, correctness)
    plotDF = data.frame(plotDF)
    plotDF$correctness = as.logical(plotDF$correctness)
    colnames(plotDF) = c("maxPred","maxActual", "correctness")
    p = ggplot(plotDF, aes(x=maxActual, y=maxPred, colour = correctness)) + geom_point() +
      scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) +
      scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = .1))
    print(p)
  }
} else {
  numBins = 20
  combined = data.frame()
  for(k in 1:validationSize){
    combined = rbind(combined, kConsensusPredictions[[k]])
  }
  #remove sample column
  combined = combined[,-1]
  combined = combined[order(combined$maxval),]
  
  trueProbs = c()
  for(i in 1:nrow(combined)){
    sample = rownames(combined)[i]
    prob = baselineProbTable[rownames(baselineProbTable) == sample, "maxvals"]
    trueProbs = c(trueProbs, prob)
  }
  
  plotDF = data.frame()
  correctness = combined$Predicted == combined$TrueCluster
  plotDF = data.frame(cbind(combined$maxval, trueProbs ,correctness))
  
  plotDF$correctness = as.logical(plotDF$correctness)
  colnames(plotDF) = c("maxPred","maxActual", "correctness")
  splitList = split(plotDF, cut(plotDF$maxPred, seq(0,1,by = .05)), drop = FALSE)
  
  plotDF = data.frame()
  averages = c()
  agreements = c()
  upErrs = c()
  lowErrs = c()
  for(i in 1:length(splitList)){
    correct = sum(splitList[[i]]$correctness)
    agreement = correct/(length(splitList[[i]]$correctness))
    agreements = c(agreements, agreement)
    
    avg = mean(splitList[[i]]$maxPred)
    averages = c(averages, avg)
    if(!is.na(correct) && (length(splitList[[i]]$correctness) != 0)){
      vals = prop.test(c(correct),length(splitList[[i]]$correctness),conf.level=0.84, correct=FALSE)[[6]]
      upper = vals[2]
      lower = vals[1]
      upErrs = c(upErrs,upper)
      lowErrs = c(lowErrs, lower)
    }
  }
  agreements = agreements[!is.na(agreements)]
  averages = averages[!is.na(averages)]
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
  
  weightedPearson = round(pearson_custom(averages,agreements,splitList), 4)
  grob2 = grobTree(textGrob(paste("Weighted Pearson Correlation : ", 
                                  weightedPearson),
                            x = 0.05, y = 0.97, 
                            hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
  
  p = ggplot(plotDF, aes(x=AverageOfBin, y=AgreementFraction)) + geom_point() +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) + 
    annotation_custom(grob) +
    annotation_custom(grob2) +
    geom_pointrange(aes(ymin=AgreementFraction-lowErrs, ymax = AgreementFraction+upErrs))
  
  filePath = paste("Plots/AgreementvsConfidencek",
                   validationSize,".forests", nforests,"probLabels",useProbabilities,".jpeg", sep="")
  if(!useSV && !useCNA){
    filePath = paste("Plots/AgreementvsConfidencek",
                     validationSize,".forests", nforests,"probLabels",useProbabilities,"_noSV_noCNA.jpeg", sep="")
  } else if(!useSV){
    filePath = paste("Plots/AgreementvsConfidencek",
                     validationSize,".forests", nforests,"probLabels",useProbabilities,"_noSV.jpeg", sep="")
  } else if(!useCNA){
    filePath = paste("Plots/AgreementvsConfidencek",
                     validationSize,".forests", nforests,"probLabels",useProbabilities,"_noCNA.jpeg", sep="")
  }
  if(downSample){
    filePath = paste("Plots/AgreementvsConfidencek",
                     validationSize,".forests", nforests,"probLabels",useProbabilities,"_DS.jpeg", sep="")
  }
  #print(filePath)
  
  #jpeg(filePath)
  #print(p)
  #dev.off()
  #print(p)
  
}



cor.test(plotDF$AverageOfBin, plotDF$AgreementFraction)


#high confidence agreement
fullPredictions = c()
fullTruths = c()
for(k in 1:validationSize){
  highConfSubset = kConsensusPredictions[[k]][kConsensusPredictions[[k]]$maxval >= .70,]
  fullPredictions = c(fullPredictions, highConfSubset$Predicted)
  fullTruths = c(fullTruths, highConfSubset$TrueCluster)
}

hConfMat = confusionMatrix(factor(fullPredictions), factor(fullTruths))


#top 70th percentile agreement
percentile = .7
top70th = combined[(nrow(combined)*(1-percentile)):nrow(combined),]
top70thConfMat = confusionMatrix(factor(top70th$Predicted), factor(top70th$TrueCluster))


source("src/ConstructReducedFeatureMatrix.R")
sortedFrame = rbind(kConsensusPredictions[[1]], kConsensusPredictions[[2]], kConsensusPredictions[[3]], kConsensusPredictions[[4]])
sortedFrame = sortedFrame[order(sortedFrame[,9], sortedFrame[,7]),]
rownames(sortedFrame) = tolower(rownames(sortedFrame))
rownames(reducedDF) = tolower(rownames(reducedDF))

bjoernFile = reducedDF[rownames(sortedFrame),]
bjoernFile = 
  cbind(sortedFrame$maxval,sortedFrame$Predicted, sortedFrame$Predicted == sortedFrame$TrueCluster,sortedFrame$TrueCluster, bjoernFile)
colnames(bjoernFile)[1:4] = c("maxVal", "cluster", "correct", "trueCluster")
filename = "WrittenFiles/RF_sortedSamples.tsv"
if(downSample){
  filename = "WrittenFiles/RF_sortedSamples_DS.tsv"
}
if(!useSV && !useCNA){
  filename = "WrittenFiles/RF_sortedSamples_noSV_noCNA.tsv"
} else if(!useSV){
  filename = "WrittenFiles/RF_sortedSamples_noSV.tsv"
} else if(!useCNA){
  filename = "WrittenFiles/RF_sortedSamples_noCNA.tsv"
}
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
v1 = which(bjoernFile$trueCluster == 1)[1]
v2 = which(bjoernFile$trueCluster == 2)[1]
v3 = which(bjoernFile$trueCluster == 3)[1]
v4 = which(bjoernFile$trueCluster == 4)[1]
v5 = which(bjoernFile$trueCluster == 5)[1]

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

filename = "Plots/RF_sortedHeatmap.jpeg"
if(downSample){
  filename = "Plots/RF_sortedHeatmap_DS.jpeg"
}
if(!useSV && !useCNA){
  filename = "Plots/RF_sortedHeatmap_noSV_noCNA.jpeg"
} else if(!useSV){
  filename = "Plots/RF_sortedHeatmap_noSV.jpeg"
} else if(!useCNA){
  filename = "Plots/RF_sortedHeatmap_noCNA.jpeg"
}
#jpeg(filename, width = 1920, height = 1080)
#print(p1)
#dev.off()
