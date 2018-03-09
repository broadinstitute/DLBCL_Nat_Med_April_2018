
#Command line calling script
# <R> <libdir>consensusReport_v5.R writeReport -o<expdata> -v<kclus> -s<bestclu> -u<allcluster> -w<markers> -p<cormatrix> -q<markersP> -r<heatmap> -t<heatmapall> -a<file.gif.2> -b<file.gif.3> -c<file.gif.4> -d<file.gif.5> -e<file.gif.6> -f<file.gif.7> -g<file.gif.8> -h<nozzle.path>
#sink(stdout(), type="message")
#require(CNTools, quietly = TRUE) || stop("Package CNTools not available")
#sink(stderr(), type="message")


#
#
writehclumiRseqReport <- function(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images){
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
    #allmems <- cbind("SampleName"=as.character(rownames(allmems)),allmems)
    allmems <- allmems[order(allmems[,2]),]
    
    #Start HTML   
    fhReport <- newReport("Clustering of miRseq mature expression: consensus hierarchical")

# --- References ---

	fullCitation <- newCitation( authors="Bmonti, S., Tamayo, P., Mesirov, J. & Golub, T.R.", title="Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data", publication="Machine Learning",number="52", pages="91-118", year="2003");
	
	webCitation <- newCitation( title="Broad Genepattern: ConsensusClustering", url="http://genepattern.broadinstitute.org/gp/pages/index.jsf" );
	
	fullCitation2 <- newCitation( authors="Rousseeuw, P.J.", title="Silhouettes: A graphical aid to the interpretation and validation of cluster analysis.", publication="J. Comput. Appl. Math.", issue="20", pages="53-65", year="1987" );
	
	webCitation2 <- newCitation( title="R silhouette package", url="http://stat.ethz.ch/R-manual/R-patched/library/cluster/html/silhouette.html" );
	
	webCitation3 <- newCitation( title="RSEM", url="http://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html" );
	

# --- Results ---

figconsensu <- newFigure( basename(consensusm),
			"The consensus matrix after clustering shows ", asParameter(nclu)," robust clusters with limited overlap between clusters." );

figcormatrix <- newFigure( basename(cormatrix),"The correlation matrix also shows ",asParameter(nclu)," robust clusters.");
figsilhouette <- newFigure( basename(kclus),
			"Silhouette width was calculated and the average silhouette width for all samples within one cluster was shown below according to different clusters (left panel). The robust cluster was pointed out by blue symbol (left panel) and the silhouette width of each sample in robust cluster was shown on right panel." );
# Removed as the figure is currently incorrect
	figheatmap <- newFigure(basename(heatmap),
		"Tumours were separated into ", asParameter(nclu), " clusters on the 
		basis of miR expression with ", asParameter(samplenum)," samples and ",
		asParameter(markergenenumber), " marker miRs. The color bar of the row 
		indicates the marker miRs for the corresponding cluster.");

figheatmapall <- newFigure( basename(heatmapall),
			"The miR expression heatmap with a standard hierarchical clustering for ",asParameter(samplenum)," samples and ",asParameter(genenumber)," most variable miRs.");


bestclufile <- bestclufile[1:10,]
tabbestclus <- newTable( bestclufile, file=basename(bestclu),
				"List of samples with ",nclu," subtypes and silhouette width." );
allmems <- as.matrix(allmems[1:10,])
#write.table(allmems,file="allmems.txt",row.names=F,quote=F)
taballclus <- newTable( allmems, file=basename(allcluster),
				"List of samples belonging to each cluster in different k clusters." );

almarkers <- almarkers[1:10,]
tabmarker <- newTable( almarkers, file=basename(markers),
				"List of marker miRs with p<= 0.05 (The positive value of column ",asParameter("difference")," means miR is upregulated in this subtype and vice versa)." );
			
	fhReport <- addToResults(fhReport,
		addTo(newSubSection("Silhouette width of each sample in robust cluster"), figsilhouette), 
# Removed as the figure is currently incorrect
		addTo(newSubSection("miR expression patterns of molecular subtype"), figheatmap, figheatmapall),
#		addTo(newSubSection("miR expression patterns of molecular subtype"), figheatmapall),
		addTo(newSubSection("Consensus and correlation matrix"), figconsensu, figcormatrix),
		addTo(newSubSection("Samples assignment with silhouette width"), tabbestclus, taballclus),
		addTo(newSubSection("Marker miRs of each subtype"), 
			newParagraph("Samples most representative of the clusters, hereby
			called ", asParameter("core samples")," were identified based on 
			positive silhouette width, indicating higher similarity to their own 
			class than to any other class member. Core samples were used to 
			select differentially expressed marker miRs for each subtype by 
			comparing the subclass versus the other subclasses, using Student's
			t-test."), tabmarker)); 


# --- Overview ---

fhReport <- addToIntroduction( fhReport,
		newParagraph( "This pipeline calculates clusters based on consensus hierarchical clustering with agglomerative average linkage", asReference( fullCitation ),",", asReference( webCitation ),". This pipeline has the following features:"), newList( isNumbered=TRUE,newParagraph( "Classify samples into consensus clusters."),newParagraph( "Determine differentially expressed marker miRs for each subtype." ) ) ); 

fhReport <- addToSummary( fhReport,
		newParagraph( "We filtered the data to ", asParameter(genenumber) ," most variable miRs. Consensus average linkage hierarchical clustering of ", asParameter(samplenum) ," samples and ", asParameter(genenumber)," miRs identified ", asParameter(nclu) ," subtypes with the stability of the clustering increasing for k = 2 to k = 8 and the average silhouette width calculation for selecting the robust clusters.") );


# --- Methods ---

# Prepare the CNMF method section
method1 <- addTo( newSubSection( "Consensus Hierarchical Clustering" ),
		newParagraph( "Consensus Hierarchical clustering is a resampling-based clustering. It provides for a method to represent the consensus across multiple runs of a clustering algorithm and to assess the stability of the discovered clusters. To this end, perturbations of the original data are simulated by resampling techniques ",asReference( fullCitation ),",", asReference( webCitation ),"." ) );

method2 <- addTo( newSubSection( "Silhouette Width" ),
		newParagraph("Silhouette width is defined as the ratio of average distance of each sample to samples in the same cluster to the smallest distance to samples not in the same cluster. If silhouette width is close to 1, it means that sample is well clustered. If silhouette width is close to -1, it means that sample is misclassified ", asReference( fullCitation2),asReference( webCitation2 ),"."));

Input <- addTo(
				newParagraph("miRseq (MIMATs) of RPM value (reads per million reads aligned to miRBase mature) with log2 transformed was as the input data for the clustering"));
				#newParagraph("miRseq of raw counts was used to calculate differentially expressed markers for each cluster", asReference(webCitation3),".")));

# Add the methods to the report
fhReport <- addToMethods( fhReport, method1, method2)

#Add the Input to the report
fhReport <- addToInput( fhReport,Input)

# Report citations
fhReport <- addToReferences( fhReport, fullCitation, webCitation, fullCitation2,webCitation2,webCitation3);



writeReport( fhReport);

}


