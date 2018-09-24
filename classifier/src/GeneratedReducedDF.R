source("src/GenerateFullDF.R")
print(getwd())

reducedDF = read.table("DataTables/Gene_Sample_Collapsed_BiologBC.1-May-2018-noSignature.addREL.CDKN2A-different.C5mutVector.txt",
                       sep = '\t', header=TRUE, row.names = 1)
reducedDF = reducedDF[,colnames(reducedDF) %in% rownames(fullDF)]
tmpDF = t(reducedDF)
rownames(tmpDF) = colnames(reducedDF)
colnames(tmpDF) = rownames(reducedDF)
reducedDF = data.frame(tmpDF)
remove(tmpDF)
