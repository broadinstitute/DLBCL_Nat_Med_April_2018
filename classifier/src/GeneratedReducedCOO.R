source("src/GenerateFullDF.R")

reducedDFCOO = read.table("DataTables/Gene_Sample_Collapsed_BiologBC.8-Aug-2018-noSignature.addREL.CDKN2A-differentorgC5mutVector.COO.txt",
                                                  sep = '\t', header=TRUE, row.names = 1)
coorow = reducedDFCOO[nrow(reducedDFCOO),]

reducedDFCOO = read.table("DataTables/Gene_Sample_Collapsed_BiologBC.1-May-2018-noSignature.addREL.CDKN2A-different.C5mutVector.txt",
                          sep = '\t', header=TRUE, row.names = 1)
reducedDFCOO = rbind(reducedDFCOO, coorow)
reducedDFCOO = reducedDFCOO[,colnames(reducedDFCOO) %in% rownames(fullDF)]
rownames(reducedDFCOO)[nrow(reducedDFCOO)] = "COO..GCB.1.unclass.2.ABC.3.na."
tmpDF = t(reducedDFCOO)
rownames(tmpDF) = colnames(reducedDFCOO)
colnames(tmpDF) = rownames(reducedDFCOO)
reducedDFCOO = data.frame(tmpDF)
reducedDFCOO = reducedDFCOO[complete.cases(reducedDFCOO),]
reducedDFCOO = reducedDFCOO[reducedDFCOO$COO..GCB.1.unclass.2.ABC.3.na. != "na",]
rown = rownames(reducedDFCOO)
reducedDFCOO <- data.frame(lapply(reducedDFCOO, as.character), stringsAsFactors=FALSE)
reducedDFCOO <- data.frame(lapply(reducedDFCOO, as.integer), stringsAsFactors=FALSE)
reducedDFCOO$COO..GCB.1.unclass.2.ABC.3.na. = as.factor(reducedDFCOO$COO..GCB.1.unclass.2.ABC.3.na.)
rownames(reducedDFCOO) = rown
write.table(reducedDFCOO, "DataTables/fixedReducedCOO.txt", sep="\t", row.names = TRUE, col.names = TRUE)
remove(tmpDF)
