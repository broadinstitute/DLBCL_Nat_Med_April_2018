source("src/GenerateLabels.R")
source("src/GenerateTestTrainSets.R")

validationSize = 4
validationSets = vector("list", validationSize)
trainTable = labels[rownames(labels) %in% trainingSet,]

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

for(k in 1:validationSize){
  fn = paste("WrittenFiles/validationSet.",k,".txt",sep="")
  fn2 = paste("ANN/DataTables/validationSet.",k,".txt",sep="")
  write.table(validationSets[[k]], fn, row.names = FALSE, col.names = FALSE, sep = '\t')
  write.table(validationSets[[k]], fn2, row.names = FALSE, col.names = FALSE, sep = '\t')
}