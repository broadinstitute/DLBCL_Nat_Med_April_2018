# <R> <libdir>Topgenes_v1.R main -s<libdir> -m<expfile> -u<selectedgenes> -o<outputprefix>

options(warn=-1)

main=function(...) {

	args <- list(...)
	for(i in 1:length(args)) {
		flag <- substring(args[[i]], 0, 2)
		value <- substring(args[[i]], 3)
		if(flag=='-m') {
			expfile=value
		}
		else if(flag=='-o'){
			outputprefix=value
		}
		else if(flag=='-s'){
			libdir=value
		}
		else if(flag=='-u') {
			selectedgenes=value
		}

	}

	exp <- read.data(expfile)
	
	# Remove up to 2 leading non-numeric columns.
	if (all(is.na(as.numeric(exp[,1])))) {
		exp=exp[, -1,drop=F]
		if (all(is.na(as.numeric(exp[,1])))) {
			exp=exp[, -1,drop=F]
		}
	}
	
	exp[exp == -Inf] <- NA
	temp <- apply(exp,2,as.numeric)
	row.names(temp) <- row.names(exp)
	temp <- sweep(temp,1,apply(temp,1,median,na.rm=T))
	exp <- temp
	#temp <- scale(t(temp), center = TRUE, scale = FALSE)
	#exp <- t(temp)

	colnames(exp) <- gsub("\\.","-",colnames(exp))
	
	na.number  <- apply(exp,2,function(x) sum(is.na(x)))
	t <- which(na.number > 0.8*dim(exp)[1])

	if(length(t) > 0 ) {exp <- exp[,-t]}
	
	if (selectedgenes == "ALL") {
		convert2gct(exp,outputprefix,libdir)
	} else {
		if (selectedgenes < 1) {selectedgenes <- as.numeric(selectedgenes)* nrow(exp)}
		topexp <- topgenes(as.data.frame(exp),selectedgenes)
		convert2gct(topexp,outputprefix,libdir)
	}
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


convert2gct <- function(exp,outputprefix,libdir) {
	#gct.rows <- nrow(exp)
	#gct.cols <- ncol(exp)
	## for the expression file ##
       if (any(is.na(exp))) {

       #install impute package
		#sink(stdout(), type="message")
		if(!require("impute",lib.loc=libdir,quietly = TRUE)){
			#sink(stderr(), type="message")
			#sink(stdout(), type="message")
			junk <- capture.output(install.packages(paste(libdir,"impute_1.20.0.tar.gz",sep=""),lib=libdir, quiet = TRUE))
			#sink(stderr(), type="message")
			#sink(stdout(), type="message")
			require("impute",lib.loc=libdir, quietly =TRUE)
			#sink(stderr(), type="message")
		}
		
		na.number  <- apply(exp,2,function(x) sum(is.na(x)))
		t <- which(na.number > 0.8*dim(exp)[1])
		if(length(t) > 0 ) {exp <- exp[,-t]}

		temp <- impute.knn(as.matrix(exp))
		exp <-cbind("GeneName"=as.character(rownames(exp)),"DESCRIPTION"=as.character(rownames(exp)),temp$data)
		
		t <- which(is.na(exp[,1]))
		if (length(t) > 0 ) {exp <- exp[-t,]}
			# write out gct file
		        gct.rows <- nrow(exp)
			gct.cols <- ncol(exp)-2
			line1<-"#1.2"
			line2<-paste(gct.rows,gct.cols,sep="\t")
			write.table(line1,file=paste(outputprefix,".expclu.gct",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
			write.table(line2,file=paste(outputprefix,".expclu.gct",sep=""),sep="\t",quote=F,row.names=F,col.names=F,append=T)
			write.table(exp,file=paste(outputprefix,".expclu.gct",sep=""),sep="\t",quote=F,row.names=F,col.names=T,append=T)

	} else {
		exp <-cbind("GeneName"=as.character(rownames(exp)),"DESCRIPTION"=as.character(rownames(exp)),exp)
		t <- which(is.na(exp[,1]))
		if (length(t) > 0 ) {exp <- exp[-t,]}
		gct.rows <- nrow(exp)
		gct.cols <- ncol(exp)-2
		# write out gct file
		line1<-"#1.2"
		line2<-paste(gct.rows,gct.cols,sep="\t")
		write.table(line1,file=paste(outputprefix,".expclu.gct",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
		write.table(line2,file=paste(outputprefix,".expclu.gct",sep=""),sep="\t",quote=F,row.names=F,col.names=F,append=T)
		write.table(exp,file=paste(outputprefix,".expclu.gct",sep=""),sep="\t",quote=F,row.names=F,col.names=T,append=T)
	}
	
}


topgenes <- function(unifiedscale,selectedgenes) {
	count=0
	mads=sort(apply(unifiedscale, 1, mad,na.rm=T), decreasing=TRUE, index.return=TRUE)
	count=1
	if (selectedgenes == "Auto") {
		selectedgenes <- findpeaks(mads$x)
	} else {
		#num <- as.numeric(selectedgenes)
		#print(num)
		#t <- which(mads$x >= quantile(mads$x,num))
		#selectedgenes <- length(t)
		selectedgenes <- as.numeric(selectedgenes)
	}
	dataset<-matrix(nrow=0, ncol=ncol(unifiedscale))
	while (count <= selectedgenes){
		tmp<-unifiedscale[mads$ix[count],]
		dataset<-rbind(dataset, tmp )
		count<-count+1
	}
	return(dataset)
}

findpeaks <- function(tmad){
	den <- density(tmad)
	id <- which(den$y == max(den$y))
	X <- den$x
	Y <- den$y
	slope <- sapply((id+1):(length(X)-1),function(i){
		tslope <- (Y[i+1]-Y[i])/(X[i+1]-X[i])
		return(tslope)
	})
	t <- which(slope >= 0)
	if(length(t) == 0) {
		t <- which(tmad >= quantile(tmad,0.9))
		tlen <- length(t)
	} else {
		peak <- X[id+t[1]-1]
		tid <- which(tmad >= X[id+t[1]-1])
		tlen <- length(tid)
	}
	return(tlen)
}
