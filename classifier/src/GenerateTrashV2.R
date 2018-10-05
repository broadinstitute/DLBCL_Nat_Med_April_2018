rm(list = ls())
seed = 200
set.seed(seed)
source("src/LoadLibraries.R")
source("src/GenerateFullDF2.R")
fullDF2 = fullDF
source("src/GenerateFullDF.R")
toremove = rownames(fullDF2)[!rownames(fullDF2) %in% rownames(fullDF)]
fullDF = fullDF2[!rownames(fullDF2) %in% toremove, ]
sampleSize = 5000

#the data frame is on a gene by gene basis. geneVec refers to the specific gene being constructed.
# geneVec first argument is the range of values it can be, pulled from the full data frame
# second argument is how many we want to generate (sample size)
# third argument is replacement
# fourth argument is the probability vector for drawing the numbers. This is calculated from the full data frame
constructedDF = data.frame()
for(i in 1:ncol(fullDF)){
  #generate the appropriate gene vector
  geneVec = sample(as.integer(names(table(fullDF[,i])/nrow(fullDF))), 
                   sampleSize, replace=TRUE, 
                   prob = as.vector(table(fullDF[,i]))/nrow(fullDF))
  
  #convert it to a matrix
  geneVec = matrix(geneVec, ncol = 1, nrow = sampleSize)
  #small conditional that appends if it exists, and if it doesn't exist then just references
  if(length(constructedDF) == 0){
    constructedDF = geneVec
  } else {
    constructedDF = cbind(constructedDF, geneVec)
  }
}

#write the data table to the corresponding file
constructedDF = data.frame(constructedDF)
colnames(constructedDF) = colnames(fullDF)
write.table(constructedDF, "DataTables/trashV2.txt", col.names = TRUE, row.names = TRUE, sep="\t")


#convert the dataframe to binary values, then do one sided exact fisher test for all columns against all other columns
binaryDF = ifelse(constructedDF > 0, 1, 0)
#initiailize an empty list to store the pvalues in
PValsAll = vector("list", length = ncol(binaryDF))

#begin loop, choose a column
for(i in 1:(ncol(binaryDF))){
  print(i)
  pvals = c()
  #check column i against all other columns
  for(j in 1:ncol(binaryDF)){
    #if we are looking at the same column, i==j, skip this iteration
    if(i == j){
      next
    } else {
      #otherwise, do an exact fisher test. Workspace refers to memory.
      #alternatve argument refers to alternative hypothesis, and can be of values two.sided, greater, or less
      #I think greater refers to right tail.
      p = fisher.test(x=factor(binaryDF[,i]), y=factor(binaryDF[,j]),
                      alternative = "greater",
                      workspace=99999999)[[1]]
      pvals = c(pvals, p)
    }
  }
  #append pvalues for gene i to the pvalsall list
  PValsAll[[i]] = pvals
}


#computation is done. This just generates plots. If there's any bias somewhere, it's above this line.
for(i in 1:length(PValsAll)){
  feature = colnames(binaryDF)[[i]]
  title = paste("Pvals for ",feature,sep="")
  df = data.frame(PValsAll[[i]])
  colnames(df)[1] = "Pvals"
  df = df[order(df$Pvals),]
  df = data.frame(df)
  xAxis = seq(1/(ncol(binaryDF)-1),1,by=1/(ncol(binaryDF)-1))
  df = cbind(df, xAxis)
  colnames(df) = c("Pvals","ExpectedP")
  df$Pvals = -log10(df$Pvals)
  df$ExpectedP = -log10(df$ExpectedP)
  
  p = ggplot(df, aes(x=ExpectedP, y=Pvals)) +
    geom_point() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5),
          title = element_text(size=20))+
    theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="grey20",size=15,angle=0,vjust=0,face="plain"),  
          axis.title.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain")) +
    geom_abline() +
    xlab("- Log10 Expected PVal") +
    ylab("- Log10 PVal")
  
  fn = paste("Plots/PValuesNullModel/QQPlots/QQPlot_",colnames(constructedDF)[i],".jpeg",sep="")
  jpeg(fn, width = 1080, height = 720)
  print(p)
  dev.off()
}
fn = "Plots/PValuesNullModel/AllPValues_2.jpeg"
jpeg(fn, width = 1080, height = 720)
print(hist(unlist(PValsAll)))
dev.off()

