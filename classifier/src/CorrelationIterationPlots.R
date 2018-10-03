rm(list = ls())
source("src/PearsonCustom.R")
source("src/LoadLibraries.R")
useSV = TRUE
useCNA = TRUE
GD = FALSE
FixedSets = FALSE
useProbs = TRUE
fullFeatures = FALSE

model = 2
if(model == 1){
  model = "ANN"
  modelAdjusted = "NN"
} else if(model == 2){
  model = "RF"
  modelAdjusted = "RF"
} else {
  model = "NMF"
  modelAdjusted = "NMF"
}

folder = paste("WrittenFiles/",model,sep="")
fileFormat = paste(modelAdjusted,"_Dist",sep="")

if(!useSV){
  fileFormat = paste(fileFormat, "_noSV",sep="")
}
if(!useCNA){
  fileFormat = paste(fileFormat, "_noCNA",sep="")
}
if(GD){
  fileFormat = paste(fileFormat, "_GD", sep="")
}
if(!useProbs){
  fileFormat = paste(fileFormat, "_BinaryLabels", sep="")
}
if(FixedSets){
  fileFormat = paste(fileFormat, "_FixedSets", sep="")
}
if(fullFeatures){
  fileFormat = paste(fileFormat, "_161F", sep="")
}
fileFormat = paste(fileFormat, "_iteration", sep="")

files = list.files(path = folder, pattern = fileFormat)
sampleIndex = sample(1:100, 1)
exampleFrame = NA

for(i in 1:length(files)){
  toOpen = paste(folder,"/",files[[i]],sep="")
  currTable = read.csv(toOpen, header = TRUE)
  rownames(currTable) = c()
  if(i == 1){
    combinedDF = currTable
  } else {
    combinedDF = rbind(combinedDF, currTable, make.row.names=TRUE)
  }
  if(i == sampleIndex){
    exampleFrame = currTable
  }
}

combinedDF = combinedDF[order(combinedDF$maxVal),]
splitListAll = split(combinedDF, cut(seq_along(rownames(combinedDF)), 10, labels = FALSE))
averages = c()
agreements = c()
upErrs = c()
lowErrs = c()

for(i in 1:length(splitListAll)){
  numCorrect = sum(splitListAll[[i]]$correct)
  agreement = numCorrect/(length(splitListAll[[i]]$correct))
  agreements = c(agreements, agreement)
  avg = mean(splitListAll[[i]]$maxVal)
  averages = c(averages, avg)
  
  vals = prop.test(c(numCorrect),length(splitListAll[[i]]$correct),conf.level=0.84, correct=FALSE)[[6]]
  upper = vals[2]
  lower = vals[1]
  upErrs = c(upErrs, upper)
  lowErrs = c(lowErrs, lower)
}

agreements = agreements[!is.na(agreements)]
averages = averages[!is.na(averages)]
upErrs = upErrs-agreements
lowErrs = agreements-lowErrs

plotDF = cbind(averages, agreements)
plotDF = data.frame(plotDF)

colnames(plotDF) = c("AverageOfBin", "AgreementFraction")

weightedPearson = round(pearson_custom(averages,agreements,splitListAll), 4)
grob2 = grobTree(textGrob(paste("Weighted Pearson Correlation : ", 
                                weightedPearson),
                          x = 0.05, y = 0.97, 
                          hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))

p = ggplot(plotDF, aes(x=AverageOfBin, y=AgreementFraction)) + geom_point() +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) + 
  annotation_custom(grob2) +
  geom_pointrange(aes(ymin=AgreementFraction-lowErrs, ymax = AgreementFraction+upErrs), shape=20)

fn = paste("Plots/OverallCorrelation_",modelAdjusted,sep="")
if(!useSV){
  fn = paste(fn, "_noSV", sep="")
}
if(!useCNA){
  fn = paste(fn, "_noCNA", sep="")
}
if(GD){
  fn = paste(fn, "_GD", sep="")
}
if(!useProbs){
  fn = paste(fn, "_BinaryLabels", sep="")
}
if(FixedSets){
  fn = paste(fn, "_FixedSets", sep="")
}
if(fullFeatures){
  fn = paste(fn, "_161F", sep="")
}
fn = paste(fn, ".jpeg", sep="")

jpeg(fn)
print(p)
dev.off()

exampleFrame = exampleFrame[order(exampleFrame$maxVal),]
splitListEx = split(exampleFrame, cut(seq_along(rownames(exampleFrame)), 10, labels = FALSE))
averages = c()
agreements = c()
upErrs = c()
lowErrs = c()

for(i in 1:length(splitListEx)){
  numCorrect = sum(splitListEx[[i]]$correct)
  agreement = numCorrect/(length(splitListEx[[i]]$correct))
  agreements = c(agreements, agreement)
  avg = mean(splitListEx[[i]]$maxVal)
  averages = c(averages, avg)
  
  vals = prop.test(c(numCorrect),length(splitListEx[[i]]$correct),conf.level=0.84, correct=FALSE)[[6]]
  upper = vals[2]
  lower = vals[1]
  upErrs = c(upErrs, upper)
  lowErrs = c(lowErrs, lower)
}

agreements = agreements[!is.na(agreements)]
averages = averages[!is.na(averages)]
upErrs = upErrs-agreements
lowErrs = agreements-lowErrs

plotDF = cbind(averages, agreements)
plotDF = data.frame(plotDF)

colnames(plotDF) = c("AverageOfBin", "AgreementFraction")

weightedPearson = round(pearson_custom(averages,agreements,splitListEx), 4)
grob2 = grobTree(textGrob(paste("Weighted Pearson Correlation : ", 
                                weightedPearson),
                          x = 0.05, y = 0.97, 
                          hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))

p = ggplot(plotDF, aes(x=AverageOfBin, y=AgreementFraction)) + geom_point() +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) + 
  annotation_custom(grob2) +
  geom_pointrange(aes(ymin=AgreementFraction-lowErrs, ymax = AgreementFraction+upErrs), shape=20)

fn = paste("Plots/ExampleCorrelation_",modelAdjusted,sep="")
if(!useSV){
  fn = paste(fn, "_noSV", sep="")
}
if(!useCNA){
  fn = paste(fn, "_noCNA", sep="")
}
if(GD){
  fn = paste(fn, "_GD", sep="")
}
if(!useProbs){
  fn = paste(fn, "_BinaryLabels", sep="")
}
if(FixedSets){
  fn = paste(fn, "_FixedSets", sep="")
}
if(fullFeatures){
  fn = paste(fn, "_161F", sep="")
}
fn = paste(fn, ".jpeg", sep="")
jpeg(fn)
print(p)
dev.off()

exampleConfMat = confusionMatrix(factor(exampleFrame$predictedCluster), factor(exampleFrame$trueCluster))
overallConfMat = confusionMatrix(factor(combinedDF$predictedCluster), factor(combinedDF$trueCluster))

top70thExample = exampleFrame[((nrow(exampleFrame)*.3):nrow(exampleFrame)),]
top70thExampleConfMat = confusionMatrix(factor(top70thExample$predictedCluster), factor(top70thExample$trueCluster))

top70thOverall = combinedDF[((nrow(combinedDF)*.3):nrow(combinedDF)),]
top70thOverallConfMat = confusionMatrix(factor(top70thOverall$predictedCluster), factor(top70thOverall$trueCluster))
