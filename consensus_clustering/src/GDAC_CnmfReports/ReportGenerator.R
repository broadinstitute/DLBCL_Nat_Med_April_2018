# Filename: ReportGenerator.R 
#  Authors: R. Zupko II
#
#  Purpose: This file checks the paramters and selects the correct report to be generated.

# Command line calling script:
# <R> <libdir>ReportGenerator.R main -o<expdata> -v<kclus> -s<bestclu> -u<allcluster> -w<markers> -p<cormatrix> -q<markersP> -r<heatmap> -t<heatmapall> -a<file.gif.2> -b<file.gif.3> -c<file.gif.4> -d<file.gif.5> -e<file.gif.6> -f<file.gif.7> -g<file.gif.8> -h<nozzle.path> -y<libdir> -z<report> 

# Main entry point for the report generator, parses the parameters provided and 
# then generates the correct report. If the correct report is not supplied the
# R stop function is invoked.
main <- function(...) {
    # Create a vector to store the images we are given, for now assume seven
    images <- c(rep("", times = 7))
    # Parse the parameters provided
    args <- list(...)
    for(i in 1:length(args)) {
        flag <- substring(args[[i]], 0, 2)
        value <- substring(args[[i]], 3, nchar(args[[i]]))
        if (flag == '-a') {
            images[1] <- value
        } else if(flag == '-b') {
            images[2] <- value
        } else if(flag == '-c') {
            images[3] <- value
        } else if(flag == '-d') {
            images[4] <- value
        } else if(flag == '-e') {
            images[5] <- value
        } else if(flag == '-f') {
            images[6] <- value
        } else if(flag == '-g') {
            images[7] <- value
        } else if(flag == '-j') {
          images[8] <- value
        } else if(flag == '-k') {
          images[9] <- value
        } else if(flag == '-h') {
            nozzle.path = value
            require(Nozzle.R1, lib.loc = nozzle.path, quietly = TRUE)
        } else if(flag == "-o") {
            expdata <- value
        } else if(flag == "-q") {
            markersP <- value
        } else if(flag == "-r") {
            heatmap <- value
        } else if(flag == "-s") {
            bestclu <- value
        } else if(flag == "-t") {
            heatmapall <- value
        } else if(flag == "-u") {
            allcluster <- value
        } else if(flag == "-v") {
            kclus <- value
        } else if(flag == "-w") {
            markers <- value
		    } else if (flag == "-y") {
			    libdir <- value
		    } else if (flag == "-z") {
			    reportType <- value
		    }
    }
  
    Inputmrna <- "The median based integrated expression data set was assembled using column-centered Level 3 data generated from Affymetrix HT-HG-U133A GeneChips, Affymetrix Human Exon 1.0 ST GeneChips, and custom designed Agilent 244k feature Gene Expression Microarrays. This data set included every gene and every samples that has been profiled on one of these platform. If a gene was only assayed on one platform, this measurement was used. If the gene was assayed on two platforms, the average of the two measurements was used; if the gene was assayed on all platforms the median measurement was used."
    Inputmrnaseq <- "mRNAseq of normalized RSEM/RPKM value with log2 transformed was as the input RNAseq data for the clustering. RSEM is used to estimate gene and transcript abundances and  these values are normalized to a fixed upper quaritile value of 1000 for gene and 300 for transcript level estimates. RPKM for a given GeneX is calculated by:  (raw read counts * 10^9) / (total reads * length of GeneX). Total reads is the lane yield after removing poor quality reads and the length of GeneX is defined as the median length of all transcripts associated with GeneX." 
    Inputmirna <- "The input file is miR array data from Agilent H-miR_8x15K/H-miR_8x15Kv2 plstforms."
    InputmiRseq <- "miRseq (at precursor level) of RPM value (reads per million reads aligned to miRBase precursor) with log2 transformed was as the input data for the clustering."
    InputmiRseqmature <- "miRseq (MIMATs) of RPM value (reads per million reads aligned to miRBase mature) with log2 transformed was as the input data for the clustering."
    Inputme <-"For a given gene, we select the probe with the maximum standard deviation across all beta values. Then we discard any probes with a standard deviation below a specified cutoff. The default cutoff is .2, but it can be tuned based on the desired output file size."
    InputCN <- "Copy number data file (the peak region with absolute value was used) = All Lesions File actual copy number part (all_lesions.conf_##.txt, where ## is the confidence level). The all lesions file is from GISTIC pipeline and summarizes the results from the GISTIC run. It contains data about the significant regions of amplification and deletion as well as which samples are amplified or deleted in each of these regions. The identified regions are listed down the first column, and the samples are listed across the first row, starting in column 10."
    InputCNT <- "Copy number data file (the threshold part was used) = All Lesions File actual copy number part (all_lesions.conf_##.txt, where ## is the confidence level). The all lesions file is from GISTIC pipeline and summarizes the results from the GISTIC run. It contains data about the significant regions of amplification and deletion as well as which samples are amplified or deleted in each of these regions. The identified regions are listed down the first column, and the samples are listed across the first row, starting in column 10."
    Inputrppa <- "The RPPA Level 3 data was used as the input for clustering; protein measurements corrected by median centering across antibodies."
    
    
    # Render the report
	if (reportType == "mrna") {
		source(paste(libdir, "/reports/cNMFreport.R", sep = ""))
		Input <- Inputmrna
		title <- "Clustering of mRNA expression: consensus NMF"
		fname <- "gene"
		FName <- "Gene"
    writecNMFReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	} else if (reportType == "mRNAseq") {
    source(paste(libdir, "/reports/cNMFreport.R", sep = ""))
	  Input <- Inputmrnaseq
    title <- "Clustering of mRNAseq gene expression: consensus NMF"
	  fname <- "gene"
	  FName <- "Gene"
	  writecNMFReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	} else if (reportType == "mirna") {
	  source(paste(libdir, "/reports/cNMFreport.R", sep = ""))
	  Input <- Inputmirna
	  title <- "Clustering of miR expression: consensus NMF"
	  fname <- "miR"
	  FName <- "miR"
	  writecNMFReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	} else if (reportType == "miRseq") {
	  source(paste(libdir, "/reports/cNMFreport.R", sep = ""))
	  Input <- InputmiRseq
    title <- "Clustering of miRseq precursor expression: consensus NMF"
	  fname <- "miR"
	  FName <- "miR"
	  writecNMFReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	} else if (reportType == "miRseqMature") {
	  source(paste(libdir, "/reports/cNMFreport.R", sep = ""))
	  Input <- InputmiRseqmature
    title <- "Clustering of miRseq mature expression: consensus NMF"
	  fname <- "miR"
	  FName <- "miR"
	  writecNMFReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	}else if (reportType == "methylation") {
	  source(paste(libdir, "/reports/cNMFreport.R", sep = ""))
		Input <- Inputme
    title <- "Clustering of Methylation: consensus NMF"
		fname <- "gene"
		FName <- "Gene"
		writecNMFReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	} else if (reportType == "nmfCN") {
	  source(paste(libdir, "/reports/cNMFreport.R", sep = ""))
	  Input <- InputCN 
    title <- "Clustering of copy number data by focal peak region with absolute value: consensus NMF"
	  fname <- "gene"
	  FName <- "Gene"
	  writecNMFReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	} else if (reportType == "nmfCNT") {
	  source(paste(libdir, "/reports/cNMFreport.R", sep = ""))
	  Input <- InputCNT
    title <- "Clustering of copy number data by peak region with threshold value: consensus NMF"
	  fname <- "gene"
	  FName <- "Gene"
	  writecNMFReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	}	else if (reportType == "nmfRppa") {
	  source(paste(libdir, "/reports/cNMFreport.R", sep = ""))
	  Input <- Inputrppa
    title <- "Clustering of RPPA data: consensus NMF"
	  fname <- "protein"
	  FName <- "Protein"
	  writecNMFReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	} else if (reportType == "hclumRNA") {
		source(paste(libdir, "/reports/hcreport.R", sep = ""))
    Input <- Inputmrna
    title <- "Clustering of mRNA expression: consensus NMF"
		title <- "Clustering of mRNA expression: consensus hierarchical"
		fname <- "gene"
		FName <- "Gene"
		writehcluReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	} else if (reportType == "hclumRNAseq") {
	  source(paste(libdir, "/reports/hcreport.R", sep = ""))
	  Input <- Inputmrnaseq
	  title <- "Clustering of mRNAseq gene expression: consensus hierarchical"
	  fname <- "gene"
	  FName <- "Gene"
	  writehcluReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	} else if (reportType == "hclumiRNA") {
	  source(paste(libdir, "/reports/hcreport.R", sep = ""))
	  Input <- Inputmirna
	  title <- "Clustering of miR expression: consensus hierarchical"
	  fname <- "miR"
	  FName <- "miR"
	  writehcluReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	} else if (reportType == "hclumiRseq") {
	  source(paste(libdir, "/reports/hcreport.R", sep = ""))
		Input <- InputmiRseq
		title <- "Clustering of miRseq precursor expression: consensus hierarchical"
		fname <- "miR"
		FName <- "miR"
		writehcluReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)	
	} else if (reportType == "hclumiRseqMature") {
	  source(paste(libdir, "/reports/hcreport.R", sep = ""))
		Input <- InputmiRseqmature
		title <- "Clustering of miRseq mature expression: consensus hierarchical"
		fname <- "miR"
		FName <- "miR"
		writehcluReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	} else if (reportType == "hclRppa") {
	  source(paste(libdir, "/reports/hcreport.R", sep = ""))
		Input <- Inputrppa
		title <- "Clustering of RPPA data: consensus hierarchical"
		fname <- "protein"
		FName <- "Protein"
		writehcluReport(expdata,kclus,bestclu,allcluster,markers,markersP,heatmap,heatmapall, images,Input,title,fname,FName)
	} else {
		stop(paste("The report type supplied was not understood: ", reportType, sep = "")) 
	} 
}

