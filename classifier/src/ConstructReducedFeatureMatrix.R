#README
#You need to create and set useSV and useCNA for this to work as a standalone script.
#If you want it to print out the data tables it creates, there is a block at the bottom with the "verify" variable, and
#the writeFiles variable inside of that block. Set both to true, and update the paths.
#
source("src/GeneratedReducedDF.R")
source("src/GenerateFullDF2.R")
fullDF = fullDF[rownames(fullDF) %in% rownames(reducedDF),]

c2F = toupper(c("X17p.DEL", "X21q.AMP", "X9p21.3.DEL", "X9q21.13.DEL",
  "X4q35.1.DEL","X1p31.1.DEL", "X1p36.11.DEL", "X1p13.1.DEL", "X4q21.22.DEL", "X14q32.31.DEL",
  "X3p21.31.DEL", "X2p16.1.AMP", "X16q12.1.DEL", "X1p36.32.DEL", "X3q28.DEL", "X1q23.3.AMP" ,
  "X18q23.DEL", "X8q24.22.AMP", "X17q24.3.AMP", "X13q14.2.DEL", "X19p13.3.DEL", "X5q.AMP",
  "X11q.AMP", "X13q34.DEL","X11p.AMP", "X13q31.3.AMP", "X6p.AMP", "X2q22.2.DEL", 
  "X12p13.2.DEL", "X6q.DEL", "X3q28.AMP", "X11q23.3.AMP", "X1q42.12.DEL", 
  "X8q12.1.DEL", "X19q13.32.DEL", "X10q23.31.DEL", "X19Q.AMP", "X5P.AMP"))

BLC10.NOTCH2.A20.SPEN = fullDF$BCL10+fullDF$NOTCH2+fullDF$SPEN+fullDF$TNFAIP3
B2M.CD70.PD1.FAS.com = fullDF$B2M+fullDF$CD70+fullDF$FAS
TP53.biallelic = fullDF$TP53
BCL2.comp = fullDF$BCL2
Chrom.Mod.Enzymes = fullDF$CREBBP+fullDF$KMT2D+fullDF$EZH2
Pi3k.Mod = fullDF$TNFRSF14+fullDF$HVCN1+fullDF$GNA13+fullDF$PTEN
Hist.comp = fullDF$HIST1H1E+fullDF$HIST1H2BC+fullDF$HIST1H1C+fullDF$HIST1H1D+fullDF$HIST1H1B+fullDF$HIST1H2AM+fullDF$HIST1H2BK+fullDF$HIST1H2AC
NFKBALT.RASALT.JAK.STAT.ALT.RHOA.SGK1.KLH6.cd58.cd83 =  
                                              fullDF$NFKBIA+fullDF$NFKBIE+fullDF$CARD11+fullDF$BRAF+fullDF$STAT3+
                                              fullDF$RHOA+fullDF$SGK1+fullDF$KLHL6+fullDF$CD58+fullDF$CD83

SumC4.mut = fullDF$SGK1+fullDF$HIST1H1E+fullDF$NFKBIE+fullDF$BRAF+fullDF$CD83+fullDF$NFKBIA+fullDF$CD58+
            fullDF$HIST1H2BC+fullDF$STAT3+fullDF$HIST1H1C+fullDF$ZFP36L1+fullDF$KLHL6+fullDF$HIST1H1D+
            fullDF$HIST1H1B+fullDF$ETS1+fullDF$TOX+fullDF$HIST1H2AM+fullDF$HIST1H2BK+fullDF$RHOA+
            fullDF$ACTB+fullDF$LTB+fullDF$SF3B1+fullDF$CARD11+fullDF$HIST1H2AC

MYD88.L265 = fullDF$MYD88.L265

MYD88.OTHER = fullDF$MYD88.OTHER

CD79B.MYD88.TBL1XR1.BLIMP.ETV6.ZC3H12A = fullDF$CD79B+fullDF$MYD88+fullDF$TBL1XR1+fullDF$PRDM1+fullDF$ETV6+fullDF$ZC3H12A

mutationsOnly = data.frame( BLC10.NOTCH2.A20.SPEN, Chrom.Mod.Enzymes, Hist.comp, 
                            NFKBALT.RASALT.JAK.STAT.ALT.RHOA.SGK1.KLH6.cd58.cd83, SumC4.mut, MYD88.L265, MYD88.OTHER,
                            CD79B.MYD88.TBL1XR1.BLIMP.ETV6.ZC3H12A)

SVonlydf = data.frame()
if(useSV){
  SV.TP63 = fullDF$SV.TP63
  SV.BCL6 = fullDF$SV.BCL6
  #combined features
  B2M.CD70.PD1.FAS.com = B2M.CD70.PD1.FAS.com+fullDF$SV.CD274.PDCD1LG2
  BCL2.comp = BCL2.comp + fullDF$SV.BCL2
  SVonlydf = data.frame(SV.TP63, SV.BCL6)
}


CNAsonlydf = data.frame()
if(useCNA){
  X21Q.AMP = fullDF$X21Q.AMP
  
  #sumCNAsC2 not resolved
  Sum.CNAs.C2 = fullDF$X17P.DEL+fullDF$X21Q.AMP+fullDF$X9P21.3.DEL+fullDF$X9Q21.13.DEL+
    fullDF$X4Q35.1.DEL+fullDF$X1P31.1.DEL+fullDF$X1P36.11.DEL+fullDF$X1P13.1.DEL+fullDF$X4Q21.22.DEL+fullDF$X14Q32.31.DEL+
    fullDF$X3P21.31.DEL+fullDF$X2P16.1.AMP+fullDF$X16Q12.1.DEL+fullDF$X1P36.32.DEL+fullDF$X3Q28.DEL+fullDF$X1Q23.3.AMP+
    fullDF$X18Q23.DEL+fullDF$X8Q24.22.AMP+fullDF$X17Q24.3.AMP+fullDF$X13Q14.2.DEL+fullDF$X19P13.3.DEL+fullDF$X5Q.AMP+
    fullDF$X11Q.AMP+fullDF$X13Q34.DEL+fullDF$X11P.AMP+fullDF$X13Q31.3.AMP+fullDF$X6P.AMP+fullDF$X2Q22.2.DEL+
    fullDF$X12P13.2.DEL+fullDF$X6Q.DEL+fullDF$X3Q28.AMP+fullDF$X11Q23.3.AMP+fullDF$X1Q42.12.DEL+
    fullDF$X8Q12.1.DEL+fullDF$X19Q13.32.DEL+fullDF$X10Q23.31.DEL
    
    sumCNAs.C5 =  fullDF$X18Q.AMP+fullDF$X3Q.AMP+fullDF$X3P.AMP+fullDF$X18P.AMP+fullDF$X17Q25.1.DEL+fullDF$X19Q13.42.AMP+
                  fullDF$X19P13.2.DEL+fullDF$X19Q.AMP
    
    X2P16.1.AMP = fullDF$X2P16.1.AMP
    
    X9P21.3.DEL = fullDF$X9P21.3.DEL
    
    X18Q.AMP = fullDF$X18Q.AMP
    
    GENOME_DOUBLING = reducedDF$GENOME_DOUBLING
    
    #combined features
    TP53.biallelic = TP53.biallelic + fullDF$X17P.DEL
    Pi3k.Mod = Pi3k.Mod+fullDF$X10Q23.31.DEL
    
    CNAsonlydf = data.frame(X21Q.AMP, Sum.CNAs.C2, sumCNAs.C5, X2P16.1.AMP, X9P21.3.DEL, X18Q.AMP, GENOME_DOUBLING)
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
rownames(constructedDF) = rownames(reducedDF)
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
#     fn = "DataTables/fullDF161.tsv"
#     write.table(t(fullDF), fn, row.names = TRUE, col.names = NA, sep='\t')
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
reducedDF = constructedDF

