# Filename: 
#  Authors:
#
#  Purpose:

# Command line calling script:
# <R> <libdir>consensusReport_v5.R writeReport -o<expdata> -v<kclus> -s<bestclu> -u<allcluster> -w<markers> -p<cormatrix> -q<markersP> -r<heatmap> -t<heatmapall> -a<file.gif.2> -b<file.gif.3> -c<file.gif.4> -d<file.gif.5> -e<file.gif.6> -f<file.gif.7> -g<file.gif.8> -h<nozzle.path>

#
#
writehclumRNAseqReport <- function(expdata,kclus,bestclu,allcluster,markers,cormatrix,markersP,heatmap,heatmapall, images){
    # Load the relevent information tha tis needed by the report
	tgenenumber <- readLines(expdata,n=2)[2]
    genenumber <- unlist(strsplit(tgenenumber,"\t"))[1]
    bestclufile <- read.delim(bestclu, header = TRUE, sep="\t", as.is = TRUE, skip = 1)
    nnclu <- unique(bestclufile[,2])
    nclu <- length(nnclu)
    samplenum <- dim(bestclufile)[1]
    almarkers <- read.delim(markers, header = TRUE, sep = "\t", as.is = TRUE, skip = 1)
    markersig <- read.delim(markersP, header = TRUE, sep = "\t",as.is = TRUE)
    markergenenumber <- nrow(markersig)

	# Prepare the reference to the images and set the names of the members
	allfile <- images
    consensusm <- allfile[nclu-1]
    allmems <- read.delim(allcluster, sep = "\t", header=TRUE, as.is = TRUE, stringsAsFactors = FALSE)
	
	# This is a bit of a hack, but since we normally only calculate up to eight, 
	# remove any higher or lower names that are not needed
	#colnames(allmems) <- paste("K=" , seq(2 + length(which(images == "")), 8), sep = "")
    #allmems <- cbind("SampleName" = as.character(rownames(allmems)), allmems)
    allmems <- allmems[order(allmems[,2]),]
    
	# Prepare the references for inclusion
	fullCitation <- newCitation( authors="Bmonti, S., Tamayo, P., Mesirov, J. & Golub, T.R.", title="Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data", publication="Machine Learning",number="52", pages="91-118", year="2003");
	
	webCitation <- newCitation( title="Broad Genepattern: ConsensusClustering", url="http://genepattern.broadinstitute.org/gp/pages/index.jsf" );
	
	fullCitation2 <- newCitation( authors="Rousseeuw, P.J.", title="Silhouettes: A graphical aid to the interpretation and validation of cluster analysis.", publication="J. Comput. Appl. Math.", issue="20", pages="53-65", year="1987" );
	
	webCitation2 <- newCitation( title="R silhouette package", url="http://stat.ethz.ch/R-manual/R-patched/library/cluster/html/silhouette.html" );
	
	webCitation3 <- newCitation( title="RSEM", url="http://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html" );


	# Start a new report
    fhReport <- newReport("Clustering of mRNAseq gene expression: consensus hierarchical")

	# Report overview section
	fhReport <- addToIntroduction(fhReport,
			newParagraph( "This pipeline calculates clusters based on consensus hierarchical clustering with agglomerative average linkage", 
					asReference( fullCitation ),",", asReference( webCitation ),". This pipeline has the following features:"), 
			newList( isNumbered=TRUE,
					newParagraph( "Classify samples into consensus clusters."),
					newParagraph( "Determine differentially expressed marker genes for each subtype." ) ) ); 

	fhReport <- addToSummary(fhReport,
			newParagraph( "The ", asParameter(genenumber) ," most variable genes were selected. Consensus average linkage hierarchical clustering of ", 
						asParameter(samplenum) ," samples and ", asParameter(genenumber)," genes identified ", asParameter(nclu) ," 
						subtypes with the stability of the clustering increasing for k = 2 to k = 8 and the average silhouette width calculation for selecting the robust clusters.") );
	

	# Prepare the table of the best clusters
	bestclufile <- bestclufile[1:10,]
	tabbestclus <- newTable(bestclufile, file = basename(bestclu), 
		"List of samples with ", nclu, " subtypes and silhouette width.");
	# Prepare the table of the list of samples and their clusters
	allmems <- as.matrix(allmems[1:10,])
	taballclus <- newTable(allmems, file = basename(allcluster), 
		"List of samples belonging to each cluster in different k clusters.");
	# Prepare the table of gene markers
	almarkers <- almarkers[1:10,]
	tabmarker <- newTable(almarkers, file = basename(markers),
		"List of marker genes with p <= 0.05 (The positive value of column ", 
		asParameter("difference"), " means gene is upregulated in this subtype and vice versa).");
	
	# Prepare the figures for inclusion in the results
	figconsensu <- newFigure(basename(consensusm), fileHighRes = basename(consensusm),
		"The consensus matrix after clustering shows ", asParameter(nclu),
		" clusters with limited overlap between clusters." );
	figcormatrix <- newFigure(basename(cormatrix), fileHighRes = basename(cormatrix),
		"The correlation matrix also shows ",asParameter(nclu)," clusters.");
	figsilhouette <- newFigure(basename(kclus), fileHighRes = basename(kclus),
		"The silhouette width was calculated for each sample and each value of ", 
		asParameter("k"), ". The left panel shows the average silhouette width 
		across all samples for each tested ", asParameter("k"), " (left panel).
		The right panels shows assignments of clusters to samples and the 
		silhouette width of each sample for the most robust clustering.");
# Removed as the figure is currently flawed
	figheatmap <- newFigure(basename(heatmap), fileHighRes = basename(heatmap),
		"Samples were separated into ", asParameter(nclu), " clusters. Shown are ", 
		asParameter(samplenum), " samples and ", asParameter(markergenenumber),
		" marker genes. The color bar of the row indicates the marker genes for
		the corresponding cluster.");
	figheatmapall <- newFigure(basename(heatmapall), fileHighRes = basename(heatmapall),
		"Heatmap with a standard hierarchical clustering for ",
		asParameter(samplenum), " samples and the ", asParameter(genenumber), " most variable genes.");

	# Render the main results
	fhReport <- addToResults(fhReport,
# Removed as the figure is currently flawed
		addTo(newSubSection("Gene expression patterns of molecular subtypes"), figheatmap, figheatmapall),
		#addTo(newSubSection("Gene expression patterns of molecular subtypes"), figheatmapall),
		addTo(newSubSection("Consensus and correlation matrix"), figconsensu, figcormatrix),
		addTo(newSubSection("Silhouette widths"), figsilhouette), 
		addTo(newSubSection("Samples assignment with silhouette width"), tabbestclus, taballclus),
		addTo(newSubSection("Marker genes of each subtype"), 
			newParagraph("Samples most representative of the clusters, hereby 
				called ", asParameter("core samples"), " were identified based 
				on positive silhouette width, indicating higher similarity to 
				their own class than to any other class member. Core samples 
				were used to select differentially expressed marker genes for 
				each subtype by comparing the subclass versus the other 
				subclasses, using Student's t-test."), tabmarker)); 
	
	# Prepare the CNMF method section
method1 <- addTo( newSubSection( "Consensus Hierarchical Clustering" ),
		newParagraph( "Consensus Hierarchical clustering is a resampling-based clustering. It provides for a method to represent the consensus across multiple runs of a clustering algorithm and to assess the stability of the discovered clusters. To this end, perturbations of the original data are simulated by resampling techniques ",asReference( fullCitation ),",", asReference( webCitation ),"." ) );

method2 <- addTo( newSubSection( "Silhouette Width" ),
		newParagraph("Silhouette width is defined as the ratio of average distance of each sample to samples in the same cluster to the smallest distance to samples not in the same cluster. If silhouette width is close to 1, it means that sample is well clustered. If silhouette width is close to -1, it means that sample is misclassified ", asReference( fullCitation2),asReference( webCitation2 ),"."));


	Input <- addTo(
				newParagraph("mRNAseq of normalized RSEM/RPKM value with log2 transformed was as the input RNAseq data for the clustering."),
				newParagraph("RSEM", asReference(webCitation3)," is used to estimate gene and transcript abundances and  these values are normalized to a fixed upper quaritile value of 1000 for gene and 300 for transcript level estimates."),
				newParagraph("RPKM for a given GeneX is calculated by:  (raw read counts * 10^9) / (total reads * length of GeneX).  Total reads is the lane yield after removing poor quality reads and the length of GeneX is defined as the median length of all transcripts associated with GeneX."));
					
				#newParagraph("mRNAseq of raw counts was used to calculate differentially expressed markers for each cluster", asReference(webCitation3),".")));
	
	# Add the methods to the report
	fhReport <- addToMethods( fhReport, method1, method2)
	
	#Add the Input to the report
	fhReport <- addToInput( fhReport,Input)

	# Report citations
	fhReport <- addToReferences( fhReport, fullCitation, webCitation, fullCitation2,webCitation2,webCitation3);

	# Render the report
	writeReport( fhReport);
}


