
#Command line calling script
# <R> <libdir>consensusReport_v5.R writeReport -o<expdata> -v<kclus> -s<bestclu> -u<allcluster> -w<markers> -p<cormatrix> -q<markersP> -r<heatmap> -t<heatmapall> -a<file.gif.2> -b<file.gif.3> -c<file.gif.4> -d<file.gif.5> -e<file.gif.6> -f<file.gif.7> -g<file.gif.8> -h<nozzle.path>
#sink(stdout(), type="message")
#require(CNTools, quietly = TRUE) || stop("Package CNTools not available")
#sink(stderr(), type="message")


#
#
writeMethylationReport <- function(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images) {
    #require(Nozzle.R1, lib.loc=nozzle.path,quietly = TRUE)
    print(kclus)
    print(expdata)
    #check <-readLines(bestclu, n=1)
    tgenenumber <- readLines(expdata,n=2)[2]
    genenumber <- unlist(strsplit(tgenenumber,"\t"))[1]
    bestclufile <- read.delim(bestclu,header=T,sep="\t",as.is = TRUE,skip=1)
    nnclu <- unique(bestclufile[,2])
    nclu <- length(nnclu)
    samplenum <- dim(bestclufile)[1]
    #tsample <-sapply(1:length(nnclu),function(i) {return(length(which(bestclufile[,2] == i)))})
    #clustable <- cbind("cluster"=c(1:length(nnclu)),"No.samples"=tsample)
    almarkers <- read.delim(markers,header=T,sep="\t",as.is = TRUE,skip=1)
    markersig <- read.delim(markersP,header=T,sep="\t",as.is = TRUE)
    markergenenumber <- nrow(markersig)
	allfile <- images
    consensusm <- allfile[nclu-1]
    allmems <- read.delim(allcluster,sep="\t",header=T,as.is = TRUE,stringsAsFactors=FALSE)
	#colnames(allmems) <- paste("K=" , seq(2 + length(which(images == "")), 8), sep = "")
       colnames(allmems) <- paste("K=" , seq(2, 8), sep = "")
    allmems <- cbind("SampleName"=as.character(rownames(allmems)),allmems)
    allmems <- allmems[order(allmems[,2]),]
    
    #Start HTML   

    fhReport <- newReport("Clustering of Methylation: consensus NMF")

# --- References ---

fullCitation <- newCitation( authors="Brunet, J.P., Tamayo, P., Golub, T.R. & Mesirov, J.P.", title="Metagenes and molecular pattern discovery using matrix factorization", publication="Proc Natl Acad Sci U S A", issue="12", number="101", pages="4164-9", year="2004", url="http://www.ncbi.nlm.nih.gov/pubmed?term=Metagenes%20and%20molecular%20pattern%20discovery%20using%20matrix%20factorization" );

webCitation <- newCitation( title="Broad Genepattern: NMFConsensus", url="http://genepattern.broadinstitute.org/gp/pages/index.jsf" );

fullCitation2 <- newCitation( authors="Rousseeuw, P.J.", title="Silhouettes: A graphical aid to the interpretation and validation of cluster analysis.", publication="J. Comput. Appl. Math.", issue="20", pages="53-65", year="1987" );

webCitation2 <- newCitation( title="R silhouette package", url="http://stat.ethz.ch/R-manual/R-patched/library/cluster/html/silhouette.html" );


	# --- Results ---
	figconsensu <- newFigure(basename(consensusm), fileHighRes = basename(consensusm),
		"The consensus matrix after clustering shows ", asParameter(nclu)," robust clusters with limited overlap between clusters." );
	figcormatrix <- newFigure(basename(cormatrix), fileHighRes = basename(cormatrix),
		"The correlation matrix also shows ",asParameter(nclu)," robust clusters.");
	figsilhouette <- newFigure(basename(kclus), fileHighRes = basename(kclus),
		"Silhouette width was calculated and the average silhouette width for all samples within one cluster was shown below according to different clusters (left panel). The robust cluster was pointed out by blue symbol (left panel) and the silhouette width of each sample in robust cluster was shown on right panel." );
# Disabled as the figure is currently incorrect
#	figheatmap <- newFigure(basename(heatmap), 
#		fileHighRes = basename(heatmap),
#		"Tumours were separated into ", asParameter(nclu), " clusters on the 
#		basis of gene methylation with ", asParameter(samplenum), " samples and ",
#		asParameter(markergenenumber), " marker genes. The color bar of the row 
#		indicates the marker genes for the corresponding cluster.");
	figheatmapall <- newFigure(basename(heatmapall), 
		fileHighRes = basename(heatmapall),
		"The gene methylation heatmap with a standard hierarchical clustering 
		for ", asParameter(samplenum), " samples and ",
		asParameter(genenumber), " most variable genes.");


bestclufile <- bestclufile[1:10,]
tabbestclus <- newTable( bestclufile, file=basename(bestclu),
				"List of samples with ",nclu," subtypes and silhouette width." );
allmems <- as.matrix(allmems[1:10,])
#write.table(allmems,file="allmems.txt",row.names=F,quote=F)
taballclus <- newTable( allmems, file=basename(allcluster),
				"List of samples belonging to each cluster in different k clusters." );

almarkers <- almarkers[1:10,]
tabmarker <- newTable( almarkers, file=basename(markers),
				"List of marker genes with p<= 0.05 (The positive value of column ",asParameter("difference")," means gene is upregulated in this subtype and vice versa)." );
			
	fhReport <- addToResults(fhReport,
		addTo(newSubSection("Silhouette width of each sample in robust cluster"), figsilhouette), 
# Disabled as the figure is currently incorrect
#		addTo(newSubSection("Gene methylation patterns of molecular subtype"), figheatmap, figheatmapall),
		addTo(newSubSection("Gene methylation patterns of molecular subtype"), figheatmapall),
		addTo(newSubSection("Consensus and correlation matrix"), figconsensu, figcormatrix),
		addTo(newSubSection("Samples assignment with silhouette width"), tabbestclus, taballclus),
		addTo(newSubSection("Marker genes of each subtype"), 
			newParagraph("Samples most representative of the clusters, hereby 
			called ", asParameter("core samples"), " were identified based on 
			positive silhouette width, indicating higher similarity to their 
			own class than to any other class member. Core samples were used to
			select differentially expressed marker genes for each subtype by 
			comparing the subclass versus the other subclasses, using Student's
			t-test."), tabmarker)); 


# --- Overview ---

fhReport <- addToIntroduction( fhReport,
				newParagraph( "This pipeline calculates clusters based on a consensus non-negative matrix factorization (NMF) clustering method ", asReference( fullCitation ),",", asReference( webCitation ),". This pipeline has the following features:"),newList( isNumbered=TRUE, newParagraph("Convert input data set to non-negativity matrix by column rank normalization."),newParagraph( "Classify samples into consensus clusters."),newParagraph( "Determine differentially expressed marker genes for each subtype." ) ) ); 

fhReport <- addToSummary( fhReport,
				newParagraph( "The ", asParameter(genenumber) ," most variable methylated genes were selected based on variation. The variation cutoff are set for each tumor type empirically by fitting a bimodal distriution. For genes with multiple methylation probes, we chose the most variable one to represent the gene.    Consensus NMF clustering of ", asParameter(samplenum) ," samples and ", asParameter(genenumber)," genes identified ", asParameter(nclu) ," subtypes with the stability of the clustering increasing for k = 2 to k = 8 and the average silhouette width calculation for selecting the robust clusters."));



# --- Methods ---

method1 <- addTo( newSubSection( "CNMF Method" ),
				newParagraph( "Non-negative matrix factorization (NMF) is an unsupervised learning algorithm
				that has been shown to identify molecular patterns when applied to gene expression data ", asReference( fullCitation ),",", asReference( webCitation ),".
				Rather than separating gene clusters based on distance computation, NMF detects contextdependent
				patterns of gene expression in complex biological systems." ) );

method2 <- addTo( newSubSection( "Silhouette Width" ),
				newParagraph("Silhouette width is defined as the ratio of average distance of each sample to samples in the same cluster to the smallest distance to samples not in the same cluster. If silhouette width is close to 1, it means that sample is well clustered. If silhouette width is close to -1, it means that sample is misclassified ", asReference( fullCitation2),"."));


fhReport <- addToMethods( fhReport, method1, method2)

# Removed as the median's are not longer used in the report
#fhReport <- addToInput( fhReport, 
#				newParagraph( asStrong( "Median-integrated mRNA expression data set" ),
#					" The median based integrated expression data set was assembled using column-centered Level 3 data generated from Affymetrix HT-HG-U133A GeneChips, Affymetrix Human Exon 1.0 ST GeneChips, and custom designed Agilent 244k feature Gene Expression Microarrays. This data set included every gene and every samples that has been profiled on one of these platform. If a gene was only assayed on one platform, this measurement was used. If the gene was assayed on two platforms, the average of the two measurements was used; if the gene was assayed on all platforms the median measurement was used." ),
#				newParameterList( asParameter( "gene expression file" ),
#					asFilename( expdata ) ) );				



fhReport <- addToReferences( fhReport, fullCitation, webCitation, fullCitation2,webCitation2 );



writeReport( fhReport);

}


