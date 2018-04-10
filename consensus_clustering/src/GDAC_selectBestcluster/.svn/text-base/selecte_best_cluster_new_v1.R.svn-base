# <R> <libdir>selecte_best_cluster_new_v1.R main -m<measure> -u<inputexp> -v<output> -w<inputallexp> -a<file.clu.2> -b<file.clu.3> -c<file.clu.4> -d<file.clu.5> -e<file.clu.6> -f<file.clu.7> -g<file.clu.8>

#### gplots function was used to plot heatmap with color key#############
#library(gplots)
sink(stderr(), type="message")
sink(stdout(), type="message")

library(cluster)
library(Cairo)


options(warn=-1)

main=function(...) {
	args <- list(...)
	for(i in 1:length(args)) {
		flag <- substring(args[[i]], 0, 2)
		value <- substring(args[[i]], 3, nchar(args[[i]]))
		if(flag=='-m') {
			measure=value
		}
		else if(flag=='-u'){
			inputexp=value
		}
		else if(flag=='-v'){
			output=value
		}
		else if(flag=='-w'){
			inputallexp=value
		}
		else if(flag=='-a'){
			file.clu.2=value
		}
		else if(flag=='-b'){
			file.clu.3=value
		}
		else if(flag=='-c'){
			file.clu.4=value
		}
		else if(flag=='-d'){
			file.clu.5=value
		}
		else if(flag=='-e'){
			file.clu.6=value
		}
		else if(flag=='-f'){
			file.clu.7=value
		}
		else if(flag=='-g'){
			file.clu.8=value
		}
	}

	suppressPackageStartupMessages(loader())

	exp <- read.data(inputallexp)
	colnames(exp) <- gsub("\\.","-",colnames(exp))
	
	if (all(is.na(as.numeric(exp[,1])))) {
		exp=exp[, -1,drop=F]
		if (all(is.na(as.numeric(exp[,1])))) {
			exp=exp[, -1,drop=F]
		}
	
	}
	#print(dim(exp)[1])
	exp[exp == -Inf] <- NA
	selectBestcluster(exp,measure,inputexp,output,file.clu.2,file.clu.3,file.clu.4,file.clu.5,file.clu.6,file.clu.7,file.clu.8)
}


## for the expression file ##

selectBestcluster <- function(allexpression,measure,inputexp,output,file.clu.2,file.clu.3,file.clu.4,file.clu.5,file.clu.6,file.clu.7,file.clu.8) {
	#get all cluter groups
	#allfile <- list.files("../",pattern = "\\.clu", rec = TRUE,full.names=T)

	#curdir <- getwd()
	#print(curdir)
	tmmm <- ls()
	print(tmmm)
	#x<-strsplit(curdir,"\\/")[[1]]
	#write.table(x,file="test.txt")
	#str <- as.numeric(rev(x)[1])-1
	#print(str)
	#dir <-paste("../",str,sep="")
	#print(dir)
	#print(allfile)
	#print(file.clu.2)
	#print("file.clu.2")
	allfile <- c(file.clu.2,file.clu.3,file.clu.4,file.clu.5,file.clu.6,file.clu.7,file.clu.8)
	#print(allfile)
	#allfile <- ls(pattern = "\\.clu")
	#allfile <- as.factor(allfile)
	#calculate dissimilarity matrix
	dism <-read.delim(file=inputexp,sep="\t",skip=2,as.is=T,header=T,check.names = F)
	row.names(dism) <- dism[,1]
	dism <- as.matrix(dism[,-c(1:2)])

	if (measure == "Pearson") {
		fdist <- cor(dism,method="pearson")
		cormatrix <- fdist
		fdist <- 1-fdist
		
	}else if (measure == "Euclidean") {
		fdist <- dist(t(dism),method="euclidean")
		fdist <- as.matrix(fdist)
		cormatrix <- fdist
	}
		
	#calculate silhouette value for selecting Best cluster
		#snms <- NULL
		sil <- NULL
		sil_all <- NULL
		allclu <- NULL
		for (i in 1:length(allfile)) {
			clufile <-allfile[i]
			clu <- convertfile(clufile)
			#snms=sort(union(snms, rownames(clu)))
			#clu=clu[snms,,drop=F]
			clu <- clu[sort(names(clu))]
			allclu <- rbind(allclu,clu)
			fdist <- fdist[names(clu),names(clu)]
			clu_si <- silhouette(clu,fdist)
			rownames(clu_si) <- names(clu)
			sil <- cbind(sil,mean(clu_si[,"sil_width"]))
			sil_all <- cbind(sil_all,clu_si[,"sil_width"])
		}
		t <- which(sil == max(sil[-1]))
		kclus <- as.numeric(t)
		if (kclus > 5){kclus <- 1} 
                tname <- paste("K=",c(seq(2,length(allfile)+1)),sep="")
		CairoPNG(paste(output,".silfig.png",sep=""),width=1200,height=1200)
		#CairoPNG(paste(output,".silfig.png",sep=""))
		par(mai = c(0.9, 0.9, 0.2, 0.2), mfrow = c(1, 2))
		boxplot(sil,names=tname,main="Average silhouette width in each cluster",ylab="Silhouette value")
		points(kclus,sil[kclus],pch=18,col="red",bg="red",cex=2.5)
		#dev.off()
		#line1 <- as.numeric(t)+1
		#allclu <- cbind("Clustes"=seq(2,length(allfile)+1),allclu)
		#line1 <- paste("We identified ",line1," subtypes according to silhouette width.",sep="")
		#print("finished")
		clubest <- allclu[kclus,]
		silbest <- sil_all[,kclus]
		allbest <- cbind(SampleName = names(clubest), cluster = clubest,silhouetteValue = silbest)
		#allbest <- allbest[order(allbest[,2]),]
		fdist <- fdist[names(clubest),names(clubest)]
		clu_si <- silhouette(clubest,fdist)
		#CairoPNG(paste(output,".bestsilhouettefig.png",sep=""),width=600,height=1200)
	
		groupCol <- getColor(allbest[,2])
		plot(clu_si,col=groupCol,main="Silhouette width for each sample in best cluster")
		dev.off()

		#find markers
		allbest <- allbest[order(allbest[,2]),]

		cormatrix <- cormatrix[as.character(allbest[,1]),as.character(allbest[,1])]
		CairoPNG(paste(output,".cormatrix.png",sep=""),width=800,height=800)
		
		#Cairo(600, 600, file=paste(output,".cormatrix.jpg",sep=""), type="jpg", bg="white", quality=100 )


		groupCol <- getColor(allbest[,2])

		blueWhiteRedGradient <-colorRampPalette( c("blue", "white", "red" ) )

		#heatmap(cormatrix,col=rainbow(100),Rowv = NA,Colv = "Rowv",scale="none",labRow="",labCol="",xlab="Samples",ylab="Samples",ColSideColors=groupCol)

		heatmap.2(cormatrix,col=blueWhiteRedGradient(100),Rowv = FALSE,Colv = "Rowv",scale="none",labRow="",labCol="",xlab="Samples",ylab="Samples",ColSideColors=groupCol,key = TRUE,dendrogram="none",density.info="none",trace="none",breaks=seq(-1,1,0.02),symbreaks=TRUE,symkey=TRUE)

		dev.off()

		#dism <- dism[,colnames(cormatrix)]
		#CairoPNG(paste(output,".geneheatmap.png",sep=""),width=600,height=600)

		#heatmap(dism,col=rainbow(100),Rowv = NA,Colv = "Rowv",scale="none",labRow="",labCol="",xlab="Samples",ylab="Genes",ColSideColors=groupCol)
		#dev.off()


		#colnames(clubest) <- colnames(allclu)
		#tname <- paste("K=",c(seq(2,length(allfile)+1)),sep="")
		allclu <- t(allclu)
		colnames(allclu) <- tname
		allclu <- cbind("TCGA_ID"=rownames(allclu),allclu)
		line1 <- t(colnames(allbest))
		write.table(line1,file=paste(output,".bestclus.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
		write.table(allbest,file=paste(output,".bestclus.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T,append=T)
		write.table(allclu,file=paste(output,".allclusters.txt",sep=""),sep="\t",quote=F,row.names=F)
		
		tpositive <- which(as.numeric(allbest[,3]) >= 0)
		allbest_o <- allbest
		allbest <- allbest[tpositive,]
		vectors <- unique(allbest[,2])
		allttest <- NULL
		Sigttest <- NULL
		for (i in 1:length(vectors)) {
			#print(i)
			tid <- which(allbest[,2] %in% vectors[i])
			if(length(tid) > 5) {
				tumor_id <- which(colnames(allexpression) %in% allbest[tid,1])
				normal_id <- which(colnames(allexpression) %in% allbest[-tid,1])
				ttest <- getttest(allexpression,tumor_id,normal_id)
				ttest <- t(ttest)
				tna <- which(ttest[,3] %in% "NA")
				if (length(tna) > 0 ) {
				 ttest <- ttest[-tna,]
				}
				if(dim(ttest)[1] > 300 ) {
					ttest <- cbind(ttest,p.adjust(as.numeric(ttest[,2]),"BH"))
					ts <- which(as.numeric(ttest[,2]) <= 0.05 & as.numeric(ttest[,4]) <= 0.25 )
					if (length(ts) >0 ){  
						ttest <- ttest[ts,]
					}      
					ttest <- ttest[order(as.numeric(ttest[,3]),decreasing=T),]
					
					pos <- which(as.numeric(ttest[,3]) >= 0)
					
					if(length(pos) > 0 & length(pos) < dim(ttest)[1]) {
						mpos <- matrix(ttest[pos,],ncol=4)
						mneg <- matrix(ttest[-pos,],ncol=4)
						neg <- dim(mneg)[1]
						if (length(pos) >= 150) {mpos <- mpos[1:150,]}
						if(neg >= 150) {mneg <- mneg[1:150,]}
						sttest <- rbind(mpos,mneg)
						sttest <- cbind(sttest,rep(vectors[i],dim(sttest)[1]))
						ttest <- cbind(ttest,rep(vectors[i],dim(ttest)[1]))				
						allttest <- rbind(allttest,ttest)
						Sigttest <- rbind(Sigttest,sttest)
							
					} else {
						ttest <- cbind(ttest,rep(vectors[i],dim(ttest)[1]))				
						allttest <- rbind(allttest,ttest)
						Sigttest <- allttest
					}
										
							
				}else {
					ttest <- cbind(ttest,p.adjust(as.numeric(ttest[,2]),"BH"))
					ts <- which(as.numeric(ttest[,2]) <= 0.05 & as.numeric(ttest[,4]) <= 0.25 )
					if (length(ts) >0 ){
						ttest <- ttest[ts,]
					}
					ttest <-  matrix(ttest,ncol=4)	
		                	ttest <- ttest[order(as.numeric(ttest[,3]),decreasing=T),]
                    			ttest <-  matrix(ttest,ncol=4)
					ttest <- cbind(ttest,rep(vectors[i],dim(ttest)[1]))
					allttest <- rbind(allttest,ttest)
					Sigttest <- allttest

				}

			}
		}

	colnames(allttest) <- c("Hybridization REF","p","difference","q","subclass")
	colnames(Sigttest) <- c("Hybridization REF","p","difference","q","subclass")
	vv <- c("Composite Element REF","p","difference","q","subclass")
	allttest <- rbind(vv,allttest)
	write.table(allttest,file=paste(output,".subclassmarkers.txt",sep=""),sep="\t",quote=F,row.names=F)
	write.table(Sigttest,file=paste(output,".seclectedSubclassmarkers.txt",sep=""),sep="\t",quote=F,row.names=F)

	qnum <- quantile(dism,c(0.1,0.9),na.rm=T)
	dism[dism >= qnum[2]] <- qnum[2]
	dism[dism <= qnum[1]] <- qnum[1]
	dism <- dism[,allbest_o[,1]]
	
	mapexp <- allexpression[Sigttest[,1],allbest_o[,1]]
	 temp <- apply(mapexp,2,as.numeric)
     row.names(temp) <- row.names(mapexp)
    mapexp <- temp

	
	qnum <- quantile(mapexp,c(0.1,0.9),na.rm=T)
	mapexp[mapexp >= qnum[2]] <- qnum[2]
	mapexp[mapexp <= qnum[1]] <- qnum[1]
	

	#install gdata, gplot and gtools package

		#if(!require("gtools",lib.loc=libdir,quietly = TRUE)){
		#	junk <- capture.output(install.packages(paste(libdir,"gtools_2.6.2.tar.gz",sep=""),lib=libdir, quiet = TRUE))
		#	require("gtools",lib.loc=libdir, quietly =TRUE)
		#}

		#if(!require("gdata",lib.loc=libdir,quietly = TRUE)){
		#	junk <- capture.output(install.packages(paste(libdir,"gdata_2.8.2.tar.gz",sep=""),lib=libdir, quiet = TRUE))
		#	require("gdata",lib.loc=libdir, quietly =TRUE)
		#}

		#if(!require("gplots",lib.loc=libdir,quietly = TRUE)){
		#	junk <- capture.output(install.packages(paste(libdir,"gplots_2.8.0.tar.gz",sep=""),lib=libdir, quiet = TRUE))
		#	require("gplots",lib.loc=libdir, quietly =TRUE)
		#}

##### Plot heatmap figure #############
	#verticalBorder <- 150
	#horizontalBorder <- 150

	#CairoPNG(paste(output,".geneheatmap.png",sep=""),width=dim( mapexp )[2] + horizontalBorder, height=dim( mapexp )[1] + verticalBorder )

	CairoPNG(paste(output,".geneheatmap.png",sep=""),width=1800, height=1200)

	#Cairo(1800, 1200, file=paste(output,".geneheatmap.jpg",sep=""), type="jpg", bg="white", quality=100 )

###########The following two command is used to plot heatmap with color key by heatmap.2 function

	rowside <- getColor(as.character(Sigttest[,5]))
	maxValue <- max(abs(mapexp), na.rm = TRUE)
	heatmap.2(mapexp,col=blueWhiteRedGradient(100),Rowv = FALSE,Colv = "Rowv",scale="none",labRow="",labCol="",xlab="Samples",ylab="Genes",ColSideColors=groupCol,RowSideColors=rowside,key = TRUE,dendrogram="none",density.info="none",trace="none",
				na.rm = TRUE, breaks = seq(-maxValue, maxValue, (2 * maxValue) / 100))

############The following two command is used to plot heatmap without color key by heatmap function
	#rowside <- getColor(sort(Sigttest[,5],decreasing=T))
	#heatmap(mapexp,col=rainbow(100),Rowv = NA,Colv = "Rowv",scale="row",labRow="",labCol="",xlab="Samples",ylab="Genes",ColSideColors=groupCol,RowSideColors=rowside)
	
	dev.off()


########### plot heatmap figure with 1500 genes and all samples #######################
	CairoPNG(paste(output,".geneheatmaptopgenes.png",sep=""),width=2500, height=2000)

		heatmap.2(dism,col=blueWhiteRedGradient(100),Rowv = TRUE,Colv = TRUE,scale="row",labRow="",labCol="",xlab="Samples",ylab="Genes",key = TRUE,dendrogram="both",density.info="none",trace="none")
	dev.off()


}

convertfile <- function(clufile){
	#get the column names
	#print(clufile)
	exp <-read.delim(clufile,header=F,sep="\t",as.is = TRUE)
	texp <- t(exp)
	clus_index <- grep(":",texp)
	clus_number <- length(clus_index)
	clus <- NULL
	for (i in 2:clus_number) {
		clus <- c(clus,c(rep(i-1,clus_index[i]-clus_index[i-1]-1)))
	}
	clus <- c(clus,c(rep(clus_number,length(texp)-clus_index[clus_number])))
	names <- texp[-clus_index]
	#res <- t(clus)
	names(clus) <- names
	#clus <- as.data.frame(clus)
	return(clus)
}

getttest <- function(expression,tumor_id,normal_id){
		expressedGenes <- NULL
		temp <- sapply(1:nrow(expression),function(i){
			#print(i)
				na.tumor <- sum(is.na(expression[i, tumor_id]))
				na.normal <- sum(is.na(expression[i, normal_id]))
				l.tumor <- length(tumor_id)
				l.normal <- length(normal_id)
				
				#print(i)
				if (na.tumor > 0.7*l.tumor | (l.tumor-na.tumor) <3 | na.normal > 0.7*l.normal | (l.normal - na.normal) < 3) {
			#print(i)
				expressedGenes <- rbind(expressedGenes, c(GeneName=as.character(rownames(expression)[i]),
							    p = "NA",
						  difference = "NA" ))

			} else {
				#print(i)
				tTest <- t.test(as.numeric(expression[i,tumor_id]), as.numeric(expression[i,normal_id]),na.rm=TRUE)
				disstance = tTest[["estimate"]]["mean of x"]-tTest[["estimate"]]["mean of y"]
				expressedGenes <- rbind(expressedGenes, c(GeneName=as.character(rownames(expression)[i]),
							    p = tTest[["p.value"]],
						  difference =disstance ))

			}
		return(expressedGenes)
		})
	return(temp)
}

read.data=function(fn, suffix.dup=F) {
	check=readLines(fn, n=1)
	check0=strsplit(check, '\t')[[1]][1]
	if (check0=='#1.2') {
		check=readLines(fn, n=3)[3]
		mat=read.delim(fn, as.is=T, skip=3, header=F)
		gnms=toupper(mat[,1])
		mat=mat[, -(1:2),drop=F]
		dnms=toupper(strsplit(check, '\t')[[1]][-(1:2)])
	} else {
		check2=readLines(fn, n=2)[2]
		check2=toupper(strsplit(check2, '\t')[[1]])
		nskip=1
		if (grepl('COMPOSITE *ELEMENT REF', check2[1])) nskip=2
		mat=read.delim(fn, as.is=T, skip=nskip, header=F)
		gnms=toupper(mat[,1])
		mat=mat[, -1,drop=F]
		dnms=toupper(strsplit(check, '\t')[[1]][-1])
	}
	ind=grep('^\\d+-SEP$', gnms)
	if (length(ind)>0) {
		tmp=gnms[ind]
		tmp=sub('-SEP','',tmp)
		tmp=sprintf('SEPT%s', tmp)
		gnms[ind]=tmp
	}
	
	ind=grep('^\\d+-MAR$', gnms)
	if (length(ind)>0) {
		tmp=gnms[ind]
		tmp=sub('-MAR','',tmp)
		tmp=sprintf('MARCH%s', tmp)
		gnms[ind]=tmp
	}
	
	ind=grep('^\\d+-APR$', gnms)
	if (length(ind)>0) {
		tmp=gnms[ind]
		tmp=sub('-APR','',tmp)
		tmp=sprintf('APR-%s', tmp)
		gnms[ind]=tmp
	}
	
	mat=as.matrix(mat)
	if (suffix.dup) {
		dup=rep(0, nrow(mat))
		gnm2=gnms
		ii=which(duplicated(gnm2))
		dup[ii]=dup[ii]+1
		while(length(ii)>0) {
			gnm2=paste(gnms, dup, sep='__')
			ii=which(duplicated(gnm2))
			dup[ii]=dup[ii]+1
		}
		gnm2=sub('\\__0$', '', gnm2)
		rownames(mat)=gnm2
	} else rownames(mat)=gnms
	colnames(mat)=dnms
	return(mat)
}

qvalue <-  function (p = NULL, lambda = seq(0, 0.9, 0.05), pi0.method = "smoother", 
    fdr.level = NULL, robust = FALSE, gui = FALSE, smooth.df = 3, 
    smooth.log.pi0 = FALSE) {
    if (is.null(p)) {
        qvalue.gui()
        return("Launching point-and-click...")
    }
    if (gui & !interactive()) 
        gui = FALSE
    if (min(p) < 0 || max(p) > 1) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: p-values not in valid range.", 
                "\n"))), parent.frame())
        else print("ERROR: p-values not in valid range.")
        return(0)
    }
    if (length(lambda) > 1 && length(lambda) < 4) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: If length of lambda greater than 1, you need at least 4 values.", 
                "\n"))), parent.frame())
        else print("ERROR: If length of lambda greater than 1, you need at least 4 values.")
        return(0)
    }
    if (length(lambda) > 1 && (min(lambda) < 0 || max(lambda) >= 
        1)) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).", 
                "\n"))), parent.frame())
        else print("ERROR: Lambda must be within [0, 1).")
        return(0)
    }
    m <- length(p)
    if (length(lambda) == 1) {
        if (lambda < 0 || lambda >= 1) {
            if (gui) 
                eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).", 
                  "\n"))), parent.frame())
            else print("ERROR: Lambda must be within [0, 1).")
            return(0)
        }
        pi0 <- mean(p >= lambda)/(1 - lambda)
        pi0 <- min(pi0, 1)
    }
    else {
        pi0 <- rep(0, length(lambda))
        for (i in 1:length(lambda)) {
            pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
        }
        if (pi0.method == "smoother") {
            if (smooth.log.pi0) 
                pi0 <- log(pi0)
            spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
            pi0 <- predict(spi0, x = max(lambda))$y
            if (smooth.log.pi0) 
                pi0 <- exp(pi0)
            pi0 <- min(pi0, 1)
        }
        else if (pi0.method == "bootstrap") {
            minpi0 <- min(pi0)
            mse <- rep(0, length(lambda))
            pi0.boot <- rep(0, length(lambda))
            for (i in 1:100) {
                p.boot <- sample(p, size = m, replace = TRUE)
                for (i in 1:length(lambda)) {
                  pi0.boot[i] <- mean(p.boot > lambda[i])/(1 - 
                    lambda[i])
                }
                mse <- mse + (pi0.boot - minpi0)^2
            }
            pi0 <- min(pi0[mse == min(mse)])
            pi0 <- min(pi0, 1)
        }
        else {
            print("ERROR: 'pi0.method' must be one of 'smoother' or 'bootstrap'.")
            return(0)
        }
    }
    if (pi0 <= 0) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.", 
                "\n"))), parent.frame())
        else print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.")
        return(0)
    }
    if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 
        1)) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: 'fdr.level' must be within (0, 1].", 
                "\n"))), parent.frame())
        else print("ERROR: 'fdr.level' must be within (0, 1].")
        return(0)
    }
    u <- order(p)
    qvalue.rank <- function(x) {
        idx <- sort.list(x)
        fc <- factor(x)
        nl <- length(levels(fc))
        bin <- as.integer(fc)
        tbl <- tabulate(bin)
        cs <- cumsum(tbl)
        tbl <- rep(cs, tbl)
        tbl[idx] <- tbl
        return(tbl)
    }
    v <- qvalue.rank(p)
    qvalue <- pi0 * m * p/v
    if (robust) {
        qvalue <- pi0 * m * p/(v * (1 - (1 - p)^m))
    }
    qvalue[u[m]] <- min(qvalue[u[m]], 1)
    for (i in (m - 1):1) {
        qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 1)
    }
    if (!is.null(fdr.level)) {
        retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, 
            pvalues = p, fdr.level = fdr.level, significant = (qvalue <= 
                fdr.level), lambda = lambda)
    }
    else {
        retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, 
            pvalues = p, lambda = lambda)
    }
    class(retval) <- "qvalue"
    return(retval)
}

getColor <- function(subType, col = c("red","blue","green","magenta1","black","yellow","orangered","brown"), 
          cNames = c("1", "2", "3", "4","5","6","7","8")){
    names(col) <- cNames
    
    return(col[subType])
}

loader <- function(){require(gplots, quietly = TRUE)}
