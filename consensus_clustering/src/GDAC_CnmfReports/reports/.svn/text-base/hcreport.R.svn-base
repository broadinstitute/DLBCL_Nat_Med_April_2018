# Filename: 
#  Authors:
#
#  Purpose:

# Command line calling script:
# <R> <libdir>consensusReport_v5.R writeReport -o<expdata> -v<kclus> -s<bestclu> -u<allcluster> -w<markers> -p<cormatrix> -q<markersP> -r<heatmap> -t<heatmapall> -a<file.gif.2> -b<file.gif.3> -c<file.gif.4> -d<file.gif.5> -e<file.gif.6> -f<file.gif.7> -g<file.gif.8> -h<nozzle.path>

#
#
writehcluReport <- function(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall,images,Inputdes,title,fname,FName){
    # Load the relevent information tha tis needed by the report
  #FName <- simpleCap(fname)
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
	fullCitation <- newCitation( authors="Monti, S., Tamayo, P., Mesirov, J. & Golub, T.R.", title="Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data", publication="Machine Learning",number="52", pages="91-118", year="2003");
	
	webCitation <- newCitation(authors="Wilkerson M and Waltman P",title="ConsensusClusterPlus: ConsensusClusterPlus", url="http://www.bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html" );
	
	fullCitation2 <- newCitation( authors="Rousseeuw, P.J.", title="Silhouettes: A graphical aid to the interpretation and validation of cluster analysis.", publication="J. Comput. Appl. Math.", issue="20", pages="53-65", year="1987" );
	
	fullCitation3 <- newCitation(
		authors = "Brunet, J.P., Tamayo, P., Golub, T.R. & Mesirov, J.P.", 
		title="Metagenes and molecular pattern discovery using matrix factorization", 
		publication="Proc Natl Acad Sci U S A", issue="12", number="101", pages="4164-9", year="2004", 
		url="http://www.ncbi.nlm.nih.gov/pubmed?term=Metagenes%20and%20molecular%20pattern%20discovery%20using%20matrix%20factorization" );
	
	webCitation2 <- newCitation( title="R silhouette package", url="http://stat.ethz.ch/R-manual/R-patched/library/cluster/html/silhouette.html" );
	
	webCitation3 <- newCitation( title="RSEM", url="http://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html" );


	# Start a new report
    #fhReport <- newReport(asParameter(title))
    #fhReport <- newReport("Clustering of mRNAseq gene expression: consensus hierarchical")
     fhReport <- newReport(title)	
	# Report overview section
	fhReport <- addToIntroduction(fhReport,
			newParagraph( "This pipeline calculates clusters based on consensus hierarchical clustering with agglomerative ward linkage", 
					asReference( fullCitation ),",", asReference( webCitation ),". This pipeline has the following features:"), 
			newList( isNumbered=TRUE,
					newParagraph( "Classify samples into consensus clusters."),
					newParagraph( "Determine differentially expressed marker ",fname,"s for each subtype." ) ) ); 

	fhReport <- addToSummary(fhReport,
			newParagraph( "Median absolute deviation (MAD) was used to select ", asParameter(genenumber) ," most variable ",fname,"s. Consensus ward linkage hierarchical clustering of ", 
						asParameter(samplenum) ," samples and ", asParameter(genenumber)," ",fname,"s identified ", asParameter(nclu) ," 
						subtypes with the stability of the clustering increasing for k = 2 to k = 10.") );
	

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
		"List of marker ",fname,"s with p <= 0.05 (The positive value of column ", 
		asParameter("difference"), " means ",fname," is upregulated in this subtype and vice versa).");
	
	# Prepare the figures for inclusion in the results
	figconsensu <- newFigure(basename(consensusm), fileHighRes = basename(consensusm),
		"The consensus matrix after clustering shows ", asParameter(nclu),
		" clusters with limited overlap between clusters." );
	#figcormatrix <- newFigure(basename(cormatrix), fileHighRes = basename(cormatrix),
	#	"The correlation matrix also shows ",asParameter(nclu)," clusters.");
	figsilhouette <- newFigure(basename(kclus), fileHighRes = basename(kclus),
		"The silhouette width was calculated for each sample and each value of ", 
		asParameter("k"), ". The left upper panel shows the average silhouette width 
		across all samples for each tested ", asParameter("k"), " (left upper panel). 
		The left lower panel shows the Cophenetic Correlation Coefficients for each tested ", asParameter("k"), ".   
		The right panel shows assignments of clusters to samples and the ", asParameter("k"), "
	        silhouette width of each sample for the most robust clustering.");
# Removed as the figure is currently flawed
	figheatmap <- newFigure(basename(heatmap), fileHighRes = basename(heatmap),
		"Samples were separated into ", asParameter(nclu), " clusters. Shown are ", 
		asParameter(samplenum), " samples and ", asParameter(markergenenumber),
		" marker ",fname,"s. The color bar of the row indicates the marker ",fname,"s for
		the corresponding cluster.");
	figheatmapall <- newFigure(basename(heatmapall), fileHighRes = basename(heatmapall),
		"Heatmap with a standard hierarchical clustering for ",
		asParameter(samplenum), " samples and the ", asParameter(genenumber), " most variable ",fname,"s.");

	# Render the main results
	fhReport <- addToResults(fhReport,
# Removed as the figure is currently flawed
		addTo(newSubSection(FName," expression patterns of molecular subtypes"), figheatmap, figheatmapall),
		#addTo(newSubSection("Gene expression patterns of molecular subtypes"), figheatmapall),
		#addTo(newSubSection("Consensus and correlation matrix"), figconsensu, figcormatrix),
		addTo(newSubSection("Silhouette widths, Cophenetic Correlation Coefficients and Consensus matrix"), figsilhouette,figconsensu), 
		addTo(newSubSection("Samples assignment with silhouette width"), tabbestclus, taballclus),
		addTo(newSubSection("Marker ",fname,"s of each subtype"), 
			newParagraph("Samples most representative of the clusters, hereby 
				called ", asParameter("core samples"), " were identified based 
				on positive silhouette width, indicating higher similarity to 
				their own class than to any other class member. Core samples 
				were used to select differentially expressed marker ",fname,"s for 
				each subtype by comparing the subclass versus the other 
				subclasses, using Student's t-test."), tabmarker)); 
	
	# Prepare the CHC method section
method1 <- addTo( newSubSection( "Consensus Hierarchical Clustering" ),
		newParagraph( "Consensus Hierarchical clustering is a resampling-based clustering. It provides for a method to represent the consensus across multiple runs of a clustering algorithm and to assess the stability of the discovered clusters. To this end, perturbations of the original data are simulated by resampling techniques. In our analysis, the R version of ConsensusClusterPlus(v1.18.0) was used ",asReference( fullCitation ),",", asReference( webCitation ),"." ) );

method2 <- addTo( newSubSection( "Silhouette Width" ),
		newParagraph("Silhouette width is defined as the ratio of average distance of each sample to samples in the same cluster to the smallest distance to samples not in the same cluster. If silhouette width is close to 1, it means that sample is well clustered. If silhouette width is close to -1, it means that sample is misclassified ", asReference( fullCitation2),asReference( webCitation2 ),"."));

# Prepare the cophenetic coeefficent method section
	copheneticMethod <- addTo(
		newSubSection("Cophenetic Correlation Coefficient and How to select the best cluster"),
		newParagraph("The cophenetic correlation coefficient is computed as the Pearson correlation of two distance matrices:"),
		newList( isNumbered=TRUE, newParagraph("Distance between samples induced by the consensus matrix."),
		         newParagraph("Distance between samples induced by the linkage used in reordering the consensus matrix.")),
		newParagraph("The cophenetic correlation coefficients and average silhouette values are used to 
                  determine the k with the most robust clusterings. From the plot of cophenetic correlation versus k, 
                  we select modes and the the point preceding the greatest decrease in cophenetic correlation coefficient, 
                 and from these choose the k with the highest average silhouette value."));
    

Input <- addTo(newParagraph(asParameter(Inputdes)),
            newParameterList( asParameter( "Input file" ),
                                                 asFilename( expdata )))
					
				#newParagraph("mRNAseq of raw counts was used to calculate differentially expressed markers for each cluster", asReference(webCitation3),".")));
	
	# Add the methods to the report
	fhReport <- addToMethods( fhReport, method1, method2,copheneticMethod)
	
	#Add the Input to the report
	fhReport <- addToInput( fhReport,Input)

	# Report citations
	fhReport <- addToReferences( fhReport, fullCitation, webCitation, fullCitation2,fullCitation3,webCitation2,webCitation3);

	# Render the report
	writeReport( fhReport);
}


#simpleCap <- function(x) {
#  s <- strsplit(x, " ")[[1]]
#  paste(toupper(substring(s, 1,1)), substring(s, 2),
#        sep="", collapse=" ")
#}
