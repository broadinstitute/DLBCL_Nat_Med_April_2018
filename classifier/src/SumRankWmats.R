fileList = list.files("Wmatrices")
summedW = NA

for(i in 1:length(fileList)){
  filePath = paste("Wmatrices/",fileList[[i]], sep="")
  currW = read.table(filePath, sep = '\t', header = TRUE, row.names = 1)
  if(i > 2){
    summedW = summedW + currW
  } else{
    summedW = currW
  }
}

summedAmps = c()
for(i in 1:nrow(summedW)){
  summedAmps = c(summedAmps, sum(summedW[i,]))
}

summedW = cbind(summedW, summedAmps)
summedW = summedW[order(summedW$summedAmps),]
summedW = data.frame(summedW)
rownames(summedW) = make.names(rownames(summedW))
rownames(summedW) = toupper(rownames(summedW))