source("src/GenerateLabels.R")

fullDF = read.table("DataTables/DLBCL_mutation_scna_sv_matrix.no_21q22.3.21_July_2017.txt", sep='\t', header=TRUE,
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