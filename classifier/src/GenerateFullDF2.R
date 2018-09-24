source("src/GenerateLabels.R")
library(caret)

fullDF = read.table("DataTables/Gene_Sample_Matrix_Genome_Doubling_MYD88_Subtypes.08-Jan-2018.txt", sep='\t', header=TRUE,
                    row.names = 1)
fullDF = fullDF[!grepl("CCF",rownames(fullDF)),]
fullDF = fullDF[!grepl("Clonal",rownames(fullDF)),]
fullDF = fullDF[ , colSums(is.na(fullDF)) == 0]
colnames(fullDF) = toupper(colnames(fullDF))
fullDF = fullDF[,colnames(fullDF) %in% rownames(labels)]
fullDF = fullDF[complete.cases(fullDF),]
tmpDF = t(fullDF)
rownames(tmpDF) = colnames(fullDF)
colnames(tmpDF) = rownames(fullDF)
fullDF = data.frame(tmpDF)
colnames(fullDF) = toupper(colnames(fullDF))
remove(tmpDF)