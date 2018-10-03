NMF.W <- function(X,W,tol,K) {
        n.run <- 1
        n.iter <- 1000000
        eps <- 1.e-50
        N <- dim(X)[1]
        M <- dim(X)[2]
        meanX <- mean(X,na.rm=T)
        eps <- 1.e-50
        for (j in 1:n.run) {
                H <- matrix(runif(K * M)*meanX,ncol=M)
                X.ap <- W %*% H
                error.EU <- sum((X-X.ap)^2)
                #error.KL <- sum(X*log((X+eps)/(X.ap+eps))+X.ap-X)
                del <- 1
                count <- 1
                while (del >= tol & count < n.iter) {
                        H <- H * (t(W) %*% X) / (t(W)%*%(W%*%H) + eps)
                        X.ap <- W %*% H
                        del <- abs(error.EU-sum((X-X.ap)^2))
                        error.EU <- sum((X-X.ap)^2)
                        #error.KL <- sum(X*log((X+eps)/(X.ap+eps))+X.ap-X)
                        #if (count %% 1000 == 0) cat(count,error.EU,del,'\n')
                        count <- count+1
                }
        }
        return(list(W,H))
}

NMF.H <- function(X,H,tol,K) {
        n.run <- 1
        n.iter <- 1000000
        eps <- 1.e-50
        N <- dim(X)[1]
        M <- dim(X)[2]
        meanX <- mean(X,na.rm=T)
        eps <- 1.e-50
        for (j in 1:n.run) {
                W <- matrix(runif(N * K)*meanX,ncol=K)
                X.ap <- W %*% H
                error.EU <- sum((X-X.ap)^2)
                #error.KL <- sum(X*log((X+eps)/(X.ap+eps))+X.ap-X)
                del <- 1
                count <- 1
                while (del >= tol & count < n.iter) {
                        W <- W * (X %*% t(H)) / ((W%*%H) %*% t(H)+eps)
                        X.ap <- W %*% H
                        del <- abs(error.EU-sum((X-X.ap)^2))
                        error.EU <- sum((X-X.ap)^2)
                        if (count %% 1000 == 0) cat(count,error.EU,del,'\n')
                        count <- count+1
                }
        }
        return(list(W,H))
}

NMF.FN <- function(X,tol,K) {
        n.run <- 1
        n.iter <- 1000000
        eps <- 1.e-50
        N <- dim(X)[1]
        M <- dim(X)[2]
        meanX <- mean(X,na.rm=T)
        eps <- 1.e-50
        for (j in 1:n.run) {
                W <- matrix(runif(N * K)*sqrt(meanX),ncol=K)
                H <- matrix(runif(K * M)*sqrt(meanX),ncol=M)
                X.ap <- W %*% H
                error.EU <- sum((X-X.ap)^2)
                error.KL <- sum(X*log((X+eps)/(X.ap+eps))+X.ap-X)
                del <- 1
                count <- 1
                while (del >= tol & count < n.iter) {
                        W <- W * (X %*% t(H)) / ((W%*%H) %*% t(H)+eps)
                        H <- H * (t(W) %*% X) / (t(W)%*%(W%*%H) + eps)
                        X.ap <- W %*% H
                        del <- abs(error.EU-sum((X-X.ap)^2))
                        error.EU <- sum((X-X.ap)^2)
                        if (count %% 10000 == 0) cat(count,error.EU,del,'\n')
                        count <- count+1
                }
        }
        return(list(W,H,X.ap))
}

#set.seed(1)
trainingOnly = TRUE
fullDFfile = "DataTables/DLBCL_mutation_scna_sv_matrix.no_21q22.3.21_July_2017.txt"
DLBCL <- read.delim(paste(fullDFfile,sep=""),header=T,sep='\t',as.is=T,comment="#")
rownames(DLBCL) <- DLBCL[,1]
DLBCL <- DLBCL[,2:ncol(DLBCL)]
DLBCL.original <- DLBCL
DLBCL[is.na(DLBCL)] <- 0
DLBCL <- DLBCL[,colSums(DLBCL)!=0]
colnames(DLBCL) <- tolower(colnames(DLBCL))

consensus <- read.delim(paste("DataTables/DLBCL.k5.connectivity.matrix.with_headers.txt",sep=""),header=T,sep='\t',as.is=T,comment="#")
colnames(consensus) <- tolower(colnames(consensus))
rownames(consensus) <- colnames(consensus)

comm <- intersect(colnames(DLBCL),colnames(consensus))
DLBCL.comm <- DLBCL[,match(comm,colnames(DLBCL),nomatch=0)]
consensus.comm <- consensus[match(comm,colnames(consensus),nomatch=0),match(comm,colnames(consensus),nomatch=0)]
if(trainingOnly){
  #trainingSet <- read.csv("DataTables/DLBCL_train_test_sets_01May2018.txt",sep='\t')
  trainingSet <- read.csv("DataTables/currentFoldreducedTrain.txt",sep='\t')
  comm <- intersect(colnames(DLBCL), tolower(rownames(trainingSet)))
  #trainingSet$pair_id = tolower(trainingSet$pair_id)
  #comm <- intersect(colnames(DLBCL), trainingSet[trainingSet$train.test == "train","pair_id"])
  DLBCL.comm <- DLBCL[,match(comm,colnames(DLBCL),nomatch=0)]
  consensus.comm <- consensus[match(comm,colnames(consensus),nomatch=0),match(comm,colnames(consensus),nomatch=0)]
  ordering = tolower(rownames(trainingSet))
  consensus.comm <- consensus.comm[ordering,ordering]
}

############## determining H matrix based on consensus matrix
K <- 5
tol <- 1.e-05
res1 <- NMF.FN(as.matrix(consensus.comm),tol,K)
H0 <- res1[[2]]
rownames(H0) <- paste("G",seq(1:nrow(H0)),sep="")
g0 <- apply(H0,2,function(x) which.max(x)) ### original clustering membership based on the consensus matrix
H0.norm <- apply(H0,2,function(x) x/sum(x))

############## determining W matrix conditioned on H0.norm to best approximate DLBCL.comm
reduced = TRUE
tmpX <- as.matrix(DLBCL.comm)
if(reduced){
  reducedDFclassifierfile = "DataTables/currentFoldReducedTrain.txt"
  reducedDFclassifier = read.table(file = reducedDFclassifierfile, sep = '\t', header = TRUE)
  reducedDFclassifier = data.frame((reducedDFclassifier), check.names = T, check.rows = T)
  reducedDFclassifier = reducedDFclassifier[,-1]
  rownames(reducedDFclassifier) = make.names(rownames(reducedDFclassifier))
  colnames(reducedDFclassifier) = make.names(colnames(reducedDFclassifier))
  rownames(reducedDFclassifier) = tolower(rownames(reducedDFclassifier))
  reducedDFclassifier = subset(reducedDFclassifier, rownames(reducedDFclassifier) %in% comm)
  for(i in 1:ncol(reducedDFclassifier)){
    reducedDFclassifier[,i] = as.numeric(as.character(reducedDFclassifier[,i]))
  }
  tmpX <- as.matrix(t(reducedDFclassifier))
  colnames(tmpX) = rownames(reducedDFclassifier)
  rownames(tmpX) = colnames(reducedDFclassifier)
  reducedDFclassifier = data.frame(tmpX)
  reducedDFclassifier = reducedDFclassifier[,rownames(consensus.comm)]
}
tmpH <- H0.norm
res2 <- NMF.H(as.matrix(tmpX),as.matrix(tmpH),tol,K)
W0 <- res2[[1]]

############## classification
tmpX <- DLBCL.comm ### feature matrix for new samples.
if(reduced){
  reducedValidationFile = "DataTables/currentFoldReducedValidation.txt"
  reducedDFclassifier = read.table(file = reducedValidationFile, sep = '\t', header = TRUE)
  reducedDFclassifier = data.frame((reducedDFclassifier), check.names = T, check.rows = T)
  reducedDFclassifier = reducedDFclassifier[,-1]
  rownames(reducedDFclassifier) = make.names(rownames(reducedDFclassifier))
  colnames(reducedDFclassifier) = make.names(colnames(reducedDFclassifier))
  rownames(reducedDFclassifier) = tolower(rownames(reducedDFclassifier))
  for(i in 1:ncol(reducedDFclassifier)){
    reducedDFclassifier[,i] = as.numeric(as.character(reducedDFclassifier[,i]))
  }
  tmpX <- as.matrix(t(reducedDFclassifier))
  colnames(tmpX) = rownames(reducedDFclassifier)
  rownames(tmpX) = colnames(reducedDFclassifier)
}
res3 <- NMF.W(as.matrix(tmpX),as.matrix(W0),tol,K)
H1 <- res3[[2]] #### H1 is an association of new samples to clustering
H1.eps = apply(H1,2,function(x) if(sum(x) == 0){x = x+.1} else {x})
H1.norm <- apply(H1.eps,2,function(x) x/sum(x))
g1 <- apply(H1,2,function(x) which.max(x)) ### hard-partitioned clustering membership for new samples
#print(any(is.na(H1.norm)))

