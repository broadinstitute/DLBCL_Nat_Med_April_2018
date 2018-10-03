#README
#You need to create and set useSV and useCNA for this to work as a standalone script.
#If you want it to print out the data tables it creates, there is a block at the bottom with the "verify" variable, and
#the writeFiles variable inside of that block. Set both to true, and update the paths.
#

useSV = TRUE
useCNA = TRUE
GD = TRUE
fullDFTrash = read.csv("DataTables/junkSet_seed100", header = TRUE, sep="\t")

c2F = toupper(c("X17p.DEL", "X21q.AMP", "X9p21.3.DEL", "X9q21.13.DEL",
                "X4q35.1.DEL","X1p31.1.DEL", "X1p36.11.DEL", "X1p13.1.DEL", "X4q21.22.DEL", "X14q32.31.DEL",
                "X3p21.31.DEL", "X2p16.1.AMP", "X16q12.1.DEL", "X1p36.32.DEL", "X3q28.DEL", "X1q23.3.AMP" ,
                "X18q23.DEL", "X8q24.22.AMP", "X17q24.3.AMP", "X13q14.2.DEL", "X19p13.3.DEL", "X5q.AMP",
                "X11q.AMP", "X13q34.DEL","X11p.AMP", "X13q31.3.AMP", "X6p.AMP", "X2q22.2.DEL", 
                "X12p13.2.DEL", "X6q.DEL", "X3q28.AMP", "X11q23.3.AMP", "X1q42.12.DEL", 
                "X8q12.1.DEL", "X19q13.32.DEL", "X10q23.31.DEL", "X19Q.AMP", "X5P.AMP"))

BLC10.NOTCH2.A20.SPEN = fullDFTrash$BCL10+fullDFTrash$NOTCH2+fullDFTrash$SPEN+fullDFTrash$TNFAIP3
B2M.CD70.PD1.FAS.com = fullDFTrash$B2M+fullDFTrash$CD70+fullDFTrash$FAS
TP53.biallelic = fullDFTrash$TP53
BCL2.comp = fullDFTrash$BCL2
Chrom.Mod.Enzymes = fullDFTrash$CREBBP+fullDFTrash$KMT2D+fullDFTrash$EZH2
Pi3k.Mod = fullDFTrash$TNFRSF14+fullDFTrash$HVCN1+fullDFTrash$GNA13+fullDFTrash$PTEN
Hist.comp = fullDFTrash$HIST1H1E+fullDFTrash$HIST1H2BC+fullDFTrash$HIST1H1C+fullDFTrash$HIST1H1D+fullDFTrash$HIST1H1B+fullDFTrash$HIST1H2AM+fullDFTrash$HIST1H2BK+fullDFTrash$HIST1H2AC
NFKBALT.RASALT.JAK.STAT.ALT.RHOA.SGK1.KLH6.cd58.cd83 =  
  fullDFTrash$NFKBIA+fullDFTrash$NFKBIE+fullDFTrash$CARD11+fullDFTrash$BRAF+fullDFTrash$STAT3+
  fullDFTrash$RHOA+fullDFTrash$SGK1+fullDFTrash$KLHL6+fullDFTrash$CD58+fullDFTrash$CD83

SumC4.mut = fullDFTrash$SGK1+fullDFTrash$HIST1H1E+fullDFTrash$NFKBIE+fullDFTrash$BRAF+fullDFTrash$CD83+fullDFTrash$NFKBIA+fullDFTrash$CD58+
  fullDFTrash$HIST1H2BC+fullDFTrash$STAT3+fullDFTrash$HIST1H1C+fullDFTrash$ZFP36L1+fullDFTrash$KLHL6+fullDFTrash$HIST1H1D+
  fullDFTrash$HIST1H1B+fullDFTrash$ETS1+fullDFTrash$TOX+fullDFTrash$HIST1H2AM+fullDFTrash$HIST1H2BK+fullDFTrash$RHOA+
  fullDFTrash$ACTB+fullDFTrash$LTB+fullDFTrash$SF3B1+fullDFTrash$CARD11+fullDFTrash$HIST1H2AC

MYD88.L265 = fullDFTrash$MYD88.L265

MYD88.OTHER = fullDFTrash$MYD88.OTHER

CD79B.MYD88.TBL1XR1.BLIMP.ETV6.ZC3H12A = fullDFTrash$CD79B+fullDFTrash$MYD88+fullDFTrash$TBL1XR1+fullDFTrash$PRDM1+fullDFTrash$ETV6+fullDFTrash$ZC3H12A

mutationsOnly = data.frame( BLC10.NOTCH2.A20.SPEN, Chrom.Mod.Enzymes, Hist.comp, 
                            NFKBALT.RASALT.JAK.STAT.ALT.RHOA.SGK1.KLH6.cd58.cd83, SumC4.mut, MYD88.L265, MYD88.OTHER,
                            CD79B.MYD88.TBL1XR1.BLIMP.ETV6.ZC3H12A)

SVonlydf = data.frame()
if(useSV){
  SV.TP63 = fullDFTrash$SV.TP63
  SV.BCL6 = fullDFTrash$SV.BCL6
  #combined features
  B2M.CD70.PD1.FAS.com = B2M.CD70.PD1.FAS.com+fullDFTrash$SV.CD274.PDCD1LG2
  BCL2.comp = BCL2.comp + fullDFTrash$SV.BCL2
  SVonlydf = data.frame(SV.TP63, SV.BCL6)
}


CNAsonlydf = data.frame()
if(useCNA){
  X21Q.AMP = fullDFTrash$X21Q.AMP
  
  #sumCNAsC2 not resolved
  Sum.CNAs.C2 = fullDFTrash$X17P.DEL+fullDFTrash$X21Q.AMP+fullDFTrash$X9P21.3.DEL+fullDFTrash$X9Q21.13.DEL+
    fullDFTrash$X4Q35.1.DEL+fullDFTrash$X1P31.1.DEL+fullDFTrash$X1P36.11.DEL+fullDFTrash$X1P13.1.DEL+fullDFTrash$X4Q21.22.DEL+fullDFTrash$X14Q32.31.DEL+
    fullDFTrash$X3P21.31.DEL+fullDFTrash$X2P16.1.AMP+fullDFTrash$X16Q12.1.DEL+fullDFTrash$X1P36.32.DEL+fullDFTrash$X3Q28.DEL+fullDFTrash$X1Q23.3.AMP+
    fullDFTrash$X18Q23.DEL+fullDFTrash$X8Q24.22.AMP+fullDFTrash$X17Q24.3.AMP+fullDFTrash$X13Q14.2.DEL+fullDFTrash$X19P13.3.DEL+fullDFTrash$X5Q.AMP+
    fullDFTrash$X11Q.AMP+fullDFTrash$X13Q34.DEL+fullDFTrash$X11P.AMP+fullDFTrash$X13Q31.3.AMP+fullDFTrash$X6P.AMP+fullDFTrash$X2Q22.2.DEL+
    fullDFTrash$X12P13.2.DEL+fullDFTrash$X6Q.DEL+fullDFTrash$X3Q28.AMP+fullDFTrash$X11Q23.3.AMP+fullDFTrash$X1Q42.12.DEL+
    fullDFTrash$X8Q12.1.DEL+fullDFTrash$X19Q13.32.DEL+fullDFTrash$X10Q23.31.DEL
  
  sumCNAs.C5 =  fullDFTrash$X18Q.AMP+fullDFTrash$X3Q.AMP+fullDFTrash$X3P.AMP+fullDFTrash$X18P.AMP+fullDFTrash$X17Q25.1.DEL+fullDFTrash$X19Q13.42.AMP+
    fullDFTrash$X19P13.2.DEL+fullDFTrash$X19Q.AMP
  
  X2P16.1.AMP = fullDFTrash$X2P16.1.AMP
  
  X9P21.3.DEL = fullDFTrash$X9P21.3.DEL
  
  X18Q.AMP = fullDFTrash$X18Q.AMP
  
  GENOME_DOUBLING = fullDFTrash$GENOME_DOUBLING
  
  #combined features
  TP53.biallelic = TP53.biallelic + fullDFTrash$X17P.DEL
  Pi3k.Mod = Pi3k.Mod+fullDFTrash$X10Q23.31.DEL
  
  if(GD){
    CNAsonlydf = data.frame(X21Q.AMP, Sum.CNAs.C2, sumCNAs.C5, X2P16.1.AMP, X9P21.3.DEL, X18Q.AMP, GENOME_DOUBLING)
  } else {
    CNAsonlydf = data.frame(X21Q.AMP, Sum.CNAs.C2, sumCNAs.C5, X2P16.1.AMP, X9P21.3.DEL, X18Q.AMP)
  }
}

combinedDF = data.frame(B2M.CD70.PD1.FAS.com, BCL2.comp, TP53.biallelic, Pi3k.Mod)
if(!useSV && !useCNA){
  constructedDF = cbind(combinedDF, mutationsOnly) 
} else if(!useSV){
  constructedDF = cbind(CNAsonlydf, combinedDF, mutationsOnly)
} else if (!useCNA){
  constructedDF = cbind(SVonlydf, combinedDF, mutationsOnly)
} else {
  #constructedDF = reducedDF
  constructedDF = cbind(combinedDF, mutationsOnly, SVonlydf, CNAsonlydf)
}
rownames(constructedDF) = rownames(fullDFTrash)
colnames(constructedDF) = toupper(colnames(constructedDF))
# verify = FALSE
# if(verify){
#   colnames(reducedDF) = toupper(colnames(reducedDF))
#   constructedDF = constructedDF[,colnames(reducedDF)]
#   notSame = rownames(constructedDF)[constructedDF[,7] != reducedDF[,7]]
#   writeFiles = FALSE
#   if(writeFiles){
#     fn = "WrittenFiles/inconsistentSamples.tsv"
#     write.table(notSame, fn, row.names = FALSE, col.names = FALSE, sep='\t')
#     fn = "DataTables/fullDFTrash161.tsv"
#     write.table(t(fullDFTrash), fn, row.names = TRUE, col.names = NA, sep='\t')
#     fn = "DataTables/reducedDFConstructed.tsv"
#     write.table(t(constructedDF), fn, row.names = TRUE, col.names = NA, sep='\t')
#     fn = "DataTables/reducedDFOriginal.tsv"
#     write.table(t(reducedDF), fn, row.names = TRUE, col.names = NA, sep='\t')
#   }
#   for(i in 1:ncol(constructedDF)){
#     if(!identical(constructedDF[,i], reducedDF[,i])){
#       print(i)
#     }
#   }
# }
reducedDFTrash = constructedDF