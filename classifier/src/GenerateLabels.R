labels = read.table("DataTables/DLBCL_no_21q22.3.bestclus.txt", header=TRUE, row.names = 1)
labels = data.frame(labels,check.names = T)
rownames(labels) = toupper(rownames(labels))
