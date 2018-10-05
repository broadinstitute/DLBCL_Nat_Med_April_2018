#This passes the sanity test
rm(list = ls())
r = 6000
c = 161
m0 <- matrix(sample(0:1,r*c, replace=TRUE),r,c)
colnames(m0) = seq(1:ncol(m0))

PValsNull = vector("list", length = ncol(m0))
for(i in 1:(ncol(m0))){
  print(i)
  pvals = c()
  for(j in 1:ncol(m0)){
    if(i == j){
      next
    } else {
      p = fisher.test(x=factor(m0[,i]), y=factor(m0[,j]), 
                      workspace=99999999)[[1]]
      pvals = c(pvals, p)
    }
  }
  PValsNull[[i]] = pvals
}

for(i in 1:length(PValsNull)){
  feature = colnames(m0)[[i]]
  title = paste("Pvals for ",feature,sep="")
  df = data.frame(PValsNull[[i]])
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
  
  fn = paste("Plots/PValuesNullModel/PValues_RandomTest",feature,".jpeg",sep="")
  jpeg(fn, width = 1080, height = 720)
  print(p)
  dev.off()
}