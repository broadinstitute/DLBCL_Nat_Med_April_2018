#source("src/get.classifier_Tim.R")
#source("src/GenerateLabels.R")
reorient = function(g1, labels){
  g1 = data.frame(g1)
  g1 = cbind(rownames(g1),g1)
  colnames(g1)[1] = "sample"
  rownames(labels) = tolower(rownames(labels))
  trueg1 = g1[g1$sample %in% rownames(labels),]
  g1c1 = trueg1[trueg1$g1 == 1,]
  g1c2 = trueg1[trueg1$g1 == 2,]
  g1c3 = trueg1[trueg1$g1 == 3,]
  g1c4 = trueg1[trueg1$g1 == 4,]
  g1c5 = trueg1[trueg1$g1 == 5,]
  
  tab1 = table(labels[rownames(labels) %in% rownames(g1c1),'cluster'])
  
  tab2 = table(labels[rownames(labels) %in% rownames(g1c2),'cluster'])
  
  tab3 = table(labels[rownames(labels) %in% rownames(g1c3),'cluster'])
  
  tab4 = table(labels[rownames(labels) %in% rownames(g1c4),'cluster'])
  
  tab5 = table(labels[rownames(labels) %in% rownames(g1c5),'cluster'])
  
  c = c(1,2,3,4,5)
  available = c(1,2,3,4,5)
  ptr = 1
  while(length(c) > 0){
    marker = ""
    if(c[ptr] == 1){
      currTab = tab1
      marker = 1
    } else if(c[ptr] == 2){
      currTab = tab2
      marker = 2
    } else if(c[ptr] == 3){
      currTab = tab3
      marker = 3
    } else if(c[ptr] == 4){
      currTab = tab4
      marker= 4
    } else if(c[ptr] == 5){
      currTab = tab5
      marker = 5
    }
    if(length(currTab) == 0){
      if(length(c) == 1){
        if(marker == 1){
          pos1Clus = available[1]
        } else if(marker == 2){
          pos2Clus = available[1]
        } else if(marker == 3){
          pos3Clus = available[1]
        } else if(marker == 4){
          pos4Clus = available[1]
        } else if(marker == 5){
          pos5Clus = available[1]
        }
        c = c[-ptr]
        next
      } else {
        ptr = ptr+1
        next
      }
    } 
    if((length(which(currTab==max(currTab))) > 1)){
      if(length(c) == 1){
        if(marker == 1){
          pos1Clus = available[1]
        } else if(marker == 2){
          pos2Clus = available[1]
        } else if(marker == 3){
          pos3Clus = available[1]
        } else if(marker == 4){
          pos4Clus = available[1]
        } else if(marker == 5){
          pos5Clus = available[1]
        }
        c = c[-ptr]
      } else {
        if(ptr > length(c)){
          ptr = 1
        } else {
          ptr = ptr+1
        }
      }
    } else {
      if(marker == 1){
        val = as.integer(names(which(currTab==max(currTab))))
        pos1Clus = val
      } else if(marker == 2){
        val = as.integer(names(which(currTab==max(currTab))))
        pos2Clus = val
      } else if(marker == 3){
        val = as.integer(names(which(currTab==max(currTab))))
        pos3Clus = val
      } else if(marker == 4){
        val = as.integer(names(which(currTab==max(currTab))))
        pos4Clus = val
      } else if(marker == 5){
        val = as.integer(names(which(currTab==max(currTab))))
        pos5Clus = val
      }
      c = c[-ptr]
      available = available[!available %in% val]
      if(ptr > length(c)){
        ptr = 1
      }
    }
  }
  
  if(nrow(g1c1) > 0){
    g1c1[,'g1'] = pos1Clus
  }
  if(nrow(g1c2) > 0){
    g1c2[,'g1'] = pos2Clus
  }
  if(nrow(g1c3) > 0){
    g1c3[,'g1'] = pos3Clus
  }
  if(nrow(g1c4) > 0){
    g1c4[,'g1'] = pos4Clus
  }
  if(nrow(g1c5) > 0){
    g1c5[,'g1'] = pos5Clus
  }
  
  
  g1fixed = rbind(g1c1,g1c2,g1c3,g1c4,g1c5)
  g1fixed = g1fixed[order(rownames(g1fixed)),]
  return(list(g1fixed,pos1Clus,pos2Clus,pos3Clus,pos4Clus,pos5Clus))
}

#g1Fixed = reorient(g1, labels)
#rownames(labels) = tolower(rownames(labels))
#labelsSubset = labels[rownames(labels) %in% rownames(g1Fixed),]
#labelsSubset = labelsSubset[rownames(g1Fixed),]