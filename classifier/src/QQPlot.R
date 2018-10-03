fullDFTrash = read.csv("DataTables/junkSet_seed100", header = TRUE, sep="\t")
fullDFTrashBinary = ifelse(fullDFTrash > 0, 1, 0)

PValsAll = vector("list", length = ncol(fullDFTrashBinary))
for(i in 1:(ncol(fullDFTrashBinary))){
  print(i)
  pvals = c()
  for(j in 1:ncol(fullDFTrashBinary)){
    if(i == j){
      next
    } else {
      p = fisher.test(x=factor(fullDFTrashBinary[,i]), y=factor(fullDFTrashBinary[,j]), 
                      workspace=99999999)[[1]]
      pvals = c(pvals, p)
    }
  }
  PValsAll[[i]] = pvals
}

for(i in 1:length(PValsAll)){
  feature = colnames(fullDFTrashBinary)[[i]]
  title = paste("Pvals for ",feature,sep="")
  df = data.frame(PValsAll[[i]])
  colnames(df)[1] = "Pvals"
  p = ggplot(df, aes(x=Pvals)) +
    geom_histogram(fill="black", alpha=1, position="identity", bins=25) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5),
          title = element_text(size=20))+
    theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="grey20",size=15,angle=0,vjust=0,face="plain"),  
          axis.title.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
  
  fn = paste("Plots/PValuesNullModel/PValues_",feature,".jpeg",sep="")
  jpeg(fn, width = 1080, height = 720)
  print(p)
  dev.off()
}
