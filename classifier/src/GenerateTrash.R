rm(list = ls())
seed = 100
set.seed(seed)
source("src/LoadLibraries.R")
source("src/GenerateFullDF2.R")
fullDF2 = fullDF
source("src/GenerateFullDF.R")
toremove = rownames(fullDF2)[!rownames(fullDF2) %in% rownames(fullDF)]
fullDF = fullDF2[!rownames(fullDF2) %in% toremove, ]

rowsums = as.vector(rowSums(fullDF != 0))
p = ggplot() +
  geom_histogram(aes(x = rowsums), bins=max(rowsums)) +
  scale_x_continuous(breaks = seq(min(rowsums), max(rowsums), by=1)) +
  ggtitle("Row Sums") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20))
jpeg("Plots/TrueDist.jpeg")
print(p)
dev.off()

cdf = ecdf(rowsums)
distvals = c()

for(i in 1:max(rowsums)){
  distvals = c(distvals, cdf(i))
}

sampleSize = 6000
vals = runif(sampleSize,0,1)

idxs = c()
for(i in 1:length(vals)){
  for(j in length(distvals):1){
    if(j == 1){
      bin = 1
    } else if(vals[[i]] > distvals[[j]]){
      bin = j+1
      break
    }
  }
  idxs = c(idxs, bin)
}

p2 = ggplot() +
  geom_histogram(aes(x = idxs), bins=max(idxs)) +
  scale_x_continuous(breaks = seq(min(idxs), max(idxs), by=1)) +
  ggtitle("Reconstructed Dist") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20))
jpeg("Plots/ReconstructedDist.jpeg")
print(p2)
dev.off()

# mutationCDFs = vector("list", length = ncol(fullDF))
# valsCDFs = vector("list", length = ncol(fullDF))
# for(i in 1:ncol(fullDF)){
#   mutationCDFs[[i]] = ecdf(fullDF[,i])
#   currvec = c()
#   for(k in 0:max(fullDF[,i])){
#     currvec = c(currvec, mutationCDFs[[i]](k))
#   }
# }

mutationFreqs = vector("list", length = ncol(fullDF))
featureOrder = colnames(fullDF)

for(i in 1:ncol(fullDF)){
  currTable = table(fullDF[,i])/nrow(fullDF)
  mutationFreqs[[i]] = currTable
}

mainDF = data.frame()

for(i in 1:sampleSize){
  if((i %% (sampleSize/10)) == 0){
    print(i)
  }
  numMuts = idxs[[i]]
  currMuts = 0
  shuffleOrder = sample(seq(1:length(featureOrder)))
  shuffledMutationFreqs = mutationFreqs[shuffleOrder]
  shuffledFeatureOrder = featureOrder[shuffleOrder]
  currSample = matrix(0, nrow=1, ncol=ncol(fullDF))
  currSample = data.frame(currSample)
  colnames(currSample) = shuffledFeatureOrder
  
  featureIndex = 1
  while(currMuts < numMuts){
    currFeature = shuffledFeatureOrder[[featureIndex]]
    currProbs = as.vector(shuffledMutationFreqs[[featureIndex]])
    possibleValues = as.integer(names(shuffledMutationFreqs[[featureIndex]]))
    currentCall = sample(possibleValues, 1, prob = currProbs, replace = TRUE)
    if(currentCall != 0){
      if(currSample[,featureIndex] == 0){
        currMuts = currMuts + 1
        currSample[,featureIndex] = currentCall
      }
      if(featureIndex == length(shuffledFeatureOrder)){
        featureIndex = 1
      } else {
        featureIndex = featureIndex + 1
      }
    } else {
      if(featureIndex == length(shuffledFeatureOrder)){
        featureIndex = 1
      } else {
        featureIndex = featureIndex + 1
      }
    }
  }
  currSample = currSample[featureOrder]
  mainDF = rbind(mainDF, currSample)
}

for(i in 1:ncol(fullDF)){
  if(colnames(fullDF)[[i]] != colnames(mainDF)[[i]]){
    print("Column names are misaligned, breaking")
    break
  }
  tmpDF <- rbind( data.frame(Distribution="Constructed", obs=table(mainDF[,i])/nrow(mainDF)),
                  data.frame(Distribution="True", obs=table(fullDF2[,i])/nrow(fullDF)))
  
  title = paste("Constructed vs True: ", colnames(fullDF)[[i]])
  p = ggplot(tmpDF, aes(x=obs.Var1, y=obs.Freq, fill=Distribution)) +
    geom_bar(stat="identity", position="dodge") +
    xlab("Mutation Value")+
    ylab("Frequency")
  
  fn = paste("Plots/FeatureDistributions/Seed",seed,"_",colnames(fullDF)[[i]],".jpeg",sep="")
  jpeg(fn)
  print(p)
  dev.off()
}

fn = "DataTables/junkSet"
fn = paste(fn, "_seed", seed, sep="")
write.table(mainDF, fn, sep="\t")



