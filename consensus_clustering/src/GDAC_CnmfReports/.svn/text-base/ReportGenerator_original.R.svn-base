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
        } else if(flag == '-h') {
            nozzle.path = value
            require(Nozzle.R1, lib.loc = nozzle.path, quietly = TRUE)
        } else if(flag == "-o") {
            expdata <- value
        } else if(flag == "-p") {
            cormatrix <- value
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
    # Render the report
	if (reportType == "mrna") {
		source(paste(libdir, "/reports/CnmfMrnaReport.R", sep = ""))
		writeMrnaReport(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images)
	} else if (reportType == "mirna") {
		source(paste(libdir, "/reports/CnmfMirnaReport.R", sep = ""))
		writeMirnaReport(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images)
	} else if (reportType == "methylation") {
		source(paste(libdir, "/reports/CnmfMethylationReport.R", sep = ""))
		writeMethylationReport(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images)
	} else if (reportType == "mRNAseq") {
		source(paste(libdir, "/reports/CnmfMrnaseqReport.R", sep = ""))
		writemRNAseqReport(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images)
	} else if (reportType == "hclumRNAseq") {
		source(paste(libdir, "/reports/hclusMrnaseq.R", sep = ""))
		writehclumRNAseqReport(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images)
	} else if (reportType == "hclumiRseq") {
		source(paste(libdir, "/reports/hclusMiRseq.R", sep = ""))
		writehclumiRseqReport(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images)	
	} else if (reportType == "hclumiRseqMature") {
		source(paste(libdir, "/reports/hclusMiRseqMature.R", sep = ""))
		writehclumiRseqReport(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images)	
	} else if (reportType == "miRseq") {
		source(paste(libdir, "/reports/CnmfmiRseqReport.R", sep = ""))
		writemiRseqReport(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images)
	} else if (reportType == "miRseqMature") {
		source(paste(libdir, "/reports/CnmfMirnaMatureReport.R", sep = ""))
		writemiRseqReport(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images)
	} else if (reportType == "nmfCN") {
		source(paste(libdir, "/reports/CnmfCNReport.R", sep = ""))
		writenmfCNReport(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images)
	} else if (reportType == "nmfCNT") {
		source(paste(libdir, "/reports/CnmfCNReport_threshold.R", sep = ""))
		writenmfCNReport(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images)
	}	else if (reportType == "nmfRppa") {
		source(paste(libdir, "/reports/CnmfRPPAReport.R", sep = ""))
		writenmfRPPAReport(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images)
	} else if (reportType == "hclRppa") {
		source(paste(libdir, "/reports/hcluRPPAReport.R", sep = ""))
		writehcluRPPAReport(expdata, kclus, bestclu, allcluster, markers, cormatrix, markersP, heatmap, heatmapall, images)
	} else {
		stop(paste("The report type supplied was not understood: ", reportType, sep = "")) 
	} 
}

