# reorient = function(g1, labels){
#   g1 = data.frame(g1)
#   g1 = cbind(rownames(g1),g1)
#   colnames(g1)[1] = "sample"
#   rownames(labels) = tolower(rownames(labels))
#   trueg1 = g1[g1$sample %in% rownames(labels),]
#   g1c1 = trueg1[trueg1$g1 == 1,]
#   g1c2 = trueg1[trueg1$g1 == 2,]
#   g1c3 = trueg1[trueg1$g1 == 3,]
#   g1c4 = trueg1[trueg1$g1 == 4,]
#   g1c5 = trueg1[trueg1$g1 == 5,]
#   
#   tab1 = table(labels[rownames(labels) %in% rownames(g1c1),'cluster'])
#   
#   tab2 = table(labels[rownames(labels) %in% rownames(g1c2),'cluster'])
#   
#   tab3 = table(labels[rownames(labels) %in% rownames(g1c3),'cluster'])
#   
#   tab4 = table(labels[rownames(labels) %in% rownames(g1c4),'cluster'])
#   
#   tab5 = table(labels[rownames(labels) %in% rownames(g1c5),'cluster'])
#   
#   c = c(1,2,3,4,5)
#   available = c(1,2,3,4,5)
#   ptr = 1
#   while(length(c) > 0){
#     marker = ""
#     if(c[ptr] == 1){
#       currTab = tab1
#       marker = 1
#     } else if(c[ptr] == 2){
#       currTab = tab2
#       marker = 2
#     } else if(c[ptr] == 3){
#       currTab = tab3
#       marker = 3
#     } else if(c[ptr] == 4){
#       currTab = tab4
#       marker= 4
#     } else if(c[ptr] == 5){
#       currTab = tab5
#       marker = 5
#     }
#     if(length(currTab) == 0){
#       if(length(c) == 1){
#         if(marker == 1){
#           pos1Clus = available[1]
#         } else if(marker == 2){
#           pos2Clus = available[1]
#         } else if(marker == 3){
#           pos3Clus = available[1]
#         } else if(marker == 4){
#           pos4Clus = available[1]
#         } else if(marker == 5){
#           pos5Clus = available[1]
#         }
#         c = c[-ptr]
#         next
#       } else {
#         ptr = ptr+1
#         next
#       }
#     } 
#     if((length(which(currTab==max(currTab))) > 1)){
#       if(length(c) == 1){
#         if(marker == 1){
#           pos1Clus = available[1]
#         } else if(marker == 2){
#           pos2Clus = available[1]
#         } else if(marker == 3){
#           pos3Clus = available[1]
#         } else if(marker == 4){
#           pos4Clus = available[1]
#         } else if(marker == 5){
#           pos5Clus = available[1]
#         }
#         c = c[-ptr]
#       } else {
#         if(ptr > length(c)){
#           ptr = 1
#         } else {
#           ptr = ptr+1
#         }
#       }
#     } else {
#       if(marker == 1){
#         val = as.integer(names(which(currTab==max(currTab))))
#         pos1Clus = val
#       } else if(marker == 2){
#         val = as.integer(names(which(currTab==max(currTab))))
#         pos2Clus = val
#       } else if(marker == 3){
#         val = as.integer(names(which(currTab==max(currTab))))
#         pos3Clus = val
#       } else if(marker == 4){
#         val = as.integer(names(which(currTab==max(currTab))))
#         pos4Clus = val
#       } else if(marker == 5){
#         val = as.integer(names(which(currTab==max(currTab))))
#         pos5Clus = val
#       }
#       c = c[-ptr]
#       available = available[!available %in% val]
#       if(ptr > length(c)){
#         ptr = 1
#       }
#     }
#   }
#   
#   if(nrow(g1c1) > 0){
#     g1c1[,'g1'] = pos1Clus
#   }
#   if(nrow(g1c2) > 0){
#     g1c2[,'g1'] = pos2Clus
#   }
#   if(nrow(g1c3) > 0){
#     g1c3[,'g1'] = pos3Clus
#   }
#   if(nrow(g1c4) > 0){
#     g1c4[,'g1'] = pos4Clus
#   }
#   if(nrow(g1c5) > 0){
#     g1c5[,'g1'] = pos5Clus
#   }
#   
#   
#   g1fixed = rbind(g1c1,g1c2,g1c3,g1c4,g1c5)
#   g1fixed = g1fixed[order(rownames(g1fixed)),]
#   return(list(g1fixed,pos1Clus,pos2Clus,pos3Clus,pos4Clus,pos5Clus))
# }


reorient = function(g1, W0, reduced){
  if(!reduced){
    SVBCL6pos = which(max(W0[rownames(W0) == "SV.BCL6",]) == W0[rownames(W0) == "SV.BCL6",])
    BCL10pos = which(max(W0[rownames(W0) == "BCL10",]) == W0[rownames(W0) == "BCL10",])
    if(SVBCL6pos != BCL10pos){
      print("Not a strong enough association of BCL10/6 to reorient cluster 1")
    } else {
      pos1Clus = as.integer(SVBCL6pos)
    }
    
    TP53pos = which(max(W0[rownames(W0) == "TP53",]) == W0[rownames(W0) == "TP53",])
    X17PDEL = which(max(W0[rownames(W0) == "X17P.DEL",]) == W0[rownames(W0) == "X17P.DEL",])
    if(TP53pos != X17PDEL){
      print("Not a strong enough association of TP53/17PDEL to reorient cluster 2")
    } else {
      pos2Clus = as.integer(TP53pos)
    }
    
    BCL2pos = which(max(W0[rownames(W0) == "BCL2",]) == W0[rownames(W0) == "BCL2",])
    SVBCL2pos = which(max(W0[rownames(W0) == "SV.BCL2",]) == W0[rownames(W0) == "SV.BCL2",])
    if(BCL2pos != SVBCL2pos){
      print("Not a strong enough association of BCL2 to reorient cluster 3")
    } else {
      pos3Clus = as.integer(BCL2pos)
    }
    
    SGK1pos = which(max(W0[rownames(W0) == "SGK1",]) == W0[rownames(W0) == "SGK1",])
    HIST1H1Epos = which(max(W0[rownames(W0) == "HIST1H1E",]) == W0[rownames(W0) == "HIST1H1E",])
    NFKBIEpos = which(max(W0[rownames(W0) == "NFKBIE",]) == W0[rownames(W0) == "NFKBIE",])
    BRAFpos = which(max(W0[rownames(W0) == "BRAF",]) == W0[rownames(W0) == "BRAF",])
    CD83pos = which(max(W0[rownames(W0) == "CD83",]) == W0[rownames(W0) == "CD83",])
    
    tab = table(c(SGK1pos, HIST1H1Epos, NFKBIEpos, BRAFpos, CD83pos))
    pos4Clus = as.integer(names(which(max(tab) == tab))[[1]])
    
    X18Q = which(max(W0[rownames(W0) == "X18Q.AMP",]) == W0[rownames(W0) == "X18Q.AMP",])
    X3Q.AMP = which(max(W0[rownames(W0) == "X3Q.AMP",]) == W0[rownames(W0) == "X3Q.AMP",])
    
    if(X18Q != X3Q.AMP){
      print("Not a strong enough association of 18Q/3Q to reorient cluster 4")
    } else {
      pos5Clus = as.integer(X18Q)
    }
    
    vec = c(pos1Clus, pos2Clus, pos3Clus, pos4Clus, pos5Clus)
    g1fixed = vec[as.factor(g1)]
    g1fixed = data.frame(g1fixed)
    rownames(g1fixed) = names(g1)
    
    return(list(g1fixed,pos1Clus,pos2Clus,pos3Clus,pos4Clus,pos5Clus))
  } else {
    clusters = c(1,2,3,4,5)
    SVBCL6pos = which(max(W0[rownames(W0) == "SV.BCL6",]) == W0[rownames(W0) == "SV.BCL6",])
    pos1Clus = as.integer(SVBCL6pos)
    
    TP53BI = which(max(W0[rownames(W0) == "TP53.BIALLELIC",]) == W0[rownames(W0) == "TP53.BIALLELIC",])
    pos2Clus = as.integer(TP53BI)
    
    BCL2comp = which(max(W0[rownames(W0) == "BCL2.COMP",]) == W0[rownames(W0) == "BCL2.COMP",])
    pos3Clus = as.integer(BCL2comp)
    
    C5Features = which(max(W0[rownames(W0) == "CD79B.MYD88.TBL1XR1.BLIMP.ETV6.ZC3H12A",]) 
                       == W0[rownames(W0) == "CD79B.MYD88.TBL1XR1.BLIMP.ETV6.ZC3H12A",])
    pos5Clus = as.integer(C5Features)
    
    pos4Clus = clusters[!(clusters %in% c(pos1Clus, pos2Clus, pos3Clus, pos5Clus))]
    
    vec = c(pos1Clus, pos2Clus, pos3Clus, pos4Clus, pos5Clus)
    g1fixed = vec[as.factor(g1)]
    g1fixed = data.frame(g1fixed)
    rownames(g1fixed) = names(g1)
    
    return(list(g1fixed,pos1Clus,pos2Clus,pos3Clus,pos4Clus,pos5Clus))
  }
}