#
#  GenePattern Methodology for:
#
#  Metagenes and Molecular Pattern Discovery using Matrix Factorization
#  Jean-Philippe Brunet, Pablo Tamayo, Todd. R. Golub, and Jill P. Mesirov
# 
#  Authors:  Pablo Tamayo (tamayo@genome.wi.mit.edu), R. Zupko II (rzupko@broadinstitute.org)
#
#  Based on the original matlab version written by Jean-Philippe Brunet (brunet@broad.mit.edu) and
#  with additional contributions from: Ted Liefeld (liefeld@broad.mit.edu)   
#  Date:  November 27, 2003
#
#  Last change March 3, 2005: modifications to make the output more readable.
#
#  Execute from an R console window with this command:
#  source("<this file>", echo = TRUE)
#  E.g. someoutput <- mynmf2(input.ds="c:\\nmf\\all_aml.res",k.init=2,k.final=5,num.clusterings=20,maxniter=500) 
#
#  For details on the method see:
#
#  Proc. Natl. Acad. Sci. USA 2004 101: 4164-4169
#  http://www.broad.mit.edu/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=89
#
#  Input parameters
#
#   input.ds
#                       input gene expression dataset in GCT or RES format
#   k.init
#                       initial value of k
#   k.final
#                       final value of k
#   num.clusterings
#                       number of NMF clusterings to build consensus matrix
#   maxniter
#                       maximum number of NMF iterations
#   error.function
#                       NMF error function: "divergence" of "euclidean"
#   rseed
#                       random number generator seed
#   directory
#                       file directory where to store the result files
#   stopconv
#                       how many no change checks are needed to stop NMF iterations (convergence)
#   stopfreq
#                       frequency (NMF iterations) of "no change" checks 
#   non.interactive.run 
#                       flag controlling if the plots are produced interatively (Rgui and saved) or only saved in files
#   doc.string
#                       prefix to be added to the output files
#
#  Output files are (prefix with "doc.string")
#
#   params.txt 
#                       run parameters and time of execution
#   membership.gct		
#			membership results for samples at all values of K
#   cophenetic.txt 
#			cophenetic values for each K
#   cophenetic.plot.jpeg
#			plot of cophenetic for each value of K		
#   consensus.k#.gct (for each value of K)
#			consensus matrix for k=#
#   consensus.plot.k#.jpeg (for each value of K)
#			plot of consensus matrix for k=#
#   graphs.k#.jpeg (for each value of K)

# <R> <libdir>CNMF_gp_v2.R main -s<libdir> -m<expfile> -u<k.int> -v<k.final> -o<outputprefix>

options(device = "png")
options(warn = -1)

# Main entry point for the program
main <- function(...) {
	#read the input files
	args <- list(...)
	for(i in 1:length(args)) {
		flag <- substring(args[[i]], 0, 2)
		value <- substring(args[[i]], 3)
		if(flag == '-m') {
			expfile = value
		}
		else if(flag == '-o'){
			outputprefix = value
		}
		if(flag == '-s'){
			libdir = value
		}
		else if(flag == '-u') {
			k.int = value
		}
		else if(flag == '-v') {
			k.final = value
		}
	}
	# Source the files that we need
	source(paste(libdir, "/CNMF_Preprocess.R", sep = ""))
	source(paste(libdir, "/CNMF_Data_Processing.R", sep = ""))
	# Load the data set
	print('loading_data_set\n')
	MSIG.Preprocess.Dataset(input.ds = expfile, output.ds = paste(outputprefix, ".normalized.gct", sep = ""), normalization = NULL)
	norfile <- list.files(getwd(), pattern = "\\.normalized.gct", full.names = T)
	# Run the main algorithm
	outputprefix
	consensusNMF(input.ds = norfile[1], k.init = k.int, k.final = k.final, num.clusterings = 20, maxniter = 1500, error.function = "divergence", rseed = 12345678, stopconv = 40, stopfreq = 10, non.interactive.run = F, doc.string = outputprefix, libdir = libdir)
}

# Load the relevent libraries for the module.
loader <- function() {
	require(doMC)
	require(jit)
}

#
#
consensusNMF <- function(input.ds, k.init, k.final, num.clusterings, maxniter, error.function, rseed=123456789,stopconv = 40, stopfreq = 10, non.interactive.run = F, doc.string = "", libdir = "", ...) {
	# Check to make sure the error function is supported, the only one in this 
	# module is divergence
	if (error.function != "divergence") {
		stop(paste("Un-supported error function: ", error.function, sep = ""))
	}
	# Save input parameters
	directory = ""
	filename <- paste(directory, doc.string, ".params.txt", sep="", collapse="")
	
	writeParameters(filename, input.ds, k.init, k.final, num.clusterings, maxniter, error.function, rseed, directory, stopconv, stopfreq, non.interactive.run, doc.string)
	# Prepare the parameters for the run
	k.init <- as.integer(k.init)
	k.final <- as.integer(k.final)
	num.clusterings <- as.integer(num.clusterings)
	n.iter <- as.integer(maxniter)
	if (!is.na(rseed)){
    	seed <- as.integer(rseed)
	}
	D <- CNMF.read.dataset(input.ds)
	col.names <- names(D)
	A <- data.matrix(D)
	# Threshold negative values to small quantity 
	A[A < 0] <- .Machine$double.eps
	cols <- ncol(A)
	num.k <- k.final - k.init + 1
	rho <- vector(mode = "numeric", length = num.k)
	k.vector <- vector(mode = "numeric", length = num.k)
	k.index <- 1
	connect.matrix.ordered <- array(0, c(num.k, cols, cols))
	all.membership <- c()
	# Load the required libraries 
    suppressPackageStartupMessages(loader())
	# Prepare for parallel code
	registerDoMC()
	print(paste("INFO: Cores: ", getDoParWorkers()))
	print(paste("INFO: DoMC Version: ", getDoParVersion()))
	# Run for the required number of ranks
	for (k in k.init:k.final) {
		# Loop over the clusterings and calculate the divergence
		assign <- foreach(i = 1:num.clusterings, .combine = cbind) %dopar% {
			# Find the nonnegative matrix factorization
		   	NMF.out <- NMF.div(V = A, k = k, maxniter = n.iter, seed = seed + i, stopconv = stopconv, stopfreq = stopfreq)
			# Find the membership
			members <- c()
			jit(1)
			for (j in 1:cols) { 
	    		members <- c(members, order(NMF.out$H[,j], decreasing = T)[1])
			}
			return(members)
	    }
		assign <- matrix(assign, nrow = num.clusterings, byrow = TRUE)
		# Compute consensus matrix
	    connect.matrix <- matrix(0, nrow = cols, ncol = cols)
	    jit(1)
		for (i in 1:num.clusterings) {
			for (j in 1:cols) {
   		     	for (p in 1:cols) {
					if (j == p) {
						connect.matrix[j, p] <- connect.matrix[j, p] + 1
					} else if (assign[i, j] == assign[i, p]) {
                    	connect.matrix[j, p] <- connect.matrix[j, p] + 1
            	  	}
        	   	}
    		}
	    }
	    connect.matrix <- connect.matrix / num.clusterings
		# Compute the distance matrix
	    dist.matrix <- as.dist(1 - connect.matrix)
	    HC <- hclust(dist.matrix, method="average")
		# Update the rho for the index
	    dist.coph <- cophenetic(HC)
	    k.vector[k.index] <- k
	    rho[k.index] <- signif((cor(dist.matrix, dist.coph)), digits = 4)
		# Order the matrix
		jit(1)
	    for (i in 1:cols) {
	    	for (j in 1:cols) {
        		connect.matrix.ordered[k.index, i, j] <- connect.matrix[HC$order[i], HC$order[j]]
    		}
	  	}
		# Compute consensus clustering membership and update the k index
		all.membership <- cbind(all.membership, cutree(HC, k = k))
    	k.index <- k.index + 1
	}
	# Save consensus matrices in one file
	require(Cairo, quietly = TRUE)
	source(paste(libdir, "/CNMF_Plot.R", sep = ""))
	# Plot the all clusters as a single image
	CairoPNG(filename = paste(directory, doc.string, ".", "consensus.all.k.plot.png", sep=""))
	nf <- layout(matrix(c(1:8), 2, 4, byrow=T), c(1, 1), c(1, 1, 1, 1), TRUE)
	for (k in 1:num.k) { 
    	CNMF.matrix.abs.plot(connect.matrix.ordered[k,,], log = F, main = paste("k=", k.vector[k]), sub = paste("Cophenetic coef.=", rho[k]), ylab = "samples", xlab ="samples")
  	}
	y.range <- c(1 - 2*(1 - min(rho)), 1)
  	plot(k.vector, rho, main ="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, xlab = "k", ylab="Cophenetic correlation", type = "n")
	lines(k.vector, rho, type = "l", col = "black")
  	points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")
	dev.off()	
	# Plot each image individually
	for (k in 1:num.k) { 
		CairoPNG(filename = paste(directory, doc.string, ".", "consensus.plot.k", k.vector[k], ".png", sep=""), width = 600, height = 600)
    	CNMF.matrix.abs.plot(connect.matrix.ordered[k,,], log = F, main = paste("k=", k.vector[k]), ylab = "samples", xlab ="samples")
		dev.off()
  	}
	# Plot the cophenetic coefficient
	CairoPNG(filename = paste(directory, doc.string, ".cophenetic.coefficient.png", sep = ""), width = 600, height = 800)
	y.range <- c(1 - 2*(1 - min(rho)), 1)
	plot(k.vector, rho, main ="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, xlab = "k", ylab="Cophenetic correlation", type = "n")
	lines(k.vector, rho, type = "l", col = "black")
	points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")
    dev.off()
	# Prepare the results to be saved
	resultsmembership <- data.frame(all.membership)
	rownames(resultsmembership) <- col.names
	colnames(resultsmembership) <- c("membership", paste("membership.", 1:(ncol(resultsmembership) - 1), sep = ""))
	# Write the membership matrix
	filename <- paste(directory, doc.string, ".", "membership", ".txt", sep="", collapse="")
	rownames(resultsmembership) <- gsub("\\.","-", rownames(resultsmembership))
	write.table(resultsmembership, file = filename, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", col.names = TRUE, row.names = TRUE)
	# Write the cophenetic coefficient for down stream use
	filename <- paste(directory, doc.string, ".cophenetic.coefficient.txt", sep = "")
	rho <- data.frame(rho, row.names = k.init:k.final)
	write.table(rho, file = filename, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", col.names = TRUE, row.names = TRUE)
}

#
#
NMF.div <- function(V, k, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {
	# Precalcuate the correct size
	M <- ncol(V)
	size <- nrow(V) * M 
    # Randomize and initialize W and H with random numbers
	set.seed(seed)
    W <- matrix(runif(nrow(V) * k), nrow = nrow(V), ncol = k)
    H <- matrix(runif(k * M), nrow = k, ncol = M)
    # Prepare the rest of the variables
    new.membership <- vector(mode = "numeric", length = M)
	old.membership <- vector(mode = "numeric", length = M)
    no.change.count <- 0
	# Run until the max iterator or there has been convergence
	jit(1)
	for (t in 1:maxniter) {
		# Calculate the updated W and H		
        H <- H * (t(W) %*% (V / (W %*% H))) + .Machine$double.eps
        norm <- apply(W, MARGIN = 2, FUN = sum)
		jit(1)
        for (i in 1:k) {
        	H[i,] <- H[i,] / norm[i]
        }
        W <- W * ((V / (W %*% H)) %*% t(H)) + .Machine$double.eps
        norm <- apply(H, MARGIN = 1, FUN = sum)
        jit(1)
		for (i in 1:k) {
        	W[,i] <- W[,i] / norm[i]
        }
        # Check to see if there has been convergence 
		if (t %% stopfreq == 0) {
        	jit(1)
			for (j in 1:M) {
            	new.membership[j] <- order(H[,j], decreasing=T)[1]
            }
			new.membership <- as.vector(new.membership)
            if (sum(new.membership == old.membership) == M) {
            	no.change.count <- no.change.count + 1
            	if (no.change.count == stopconv) { 
					break
				}
            } else {
            	no.change.count <- 0
            }
            old.membership <- new.membership
        }
	}
    return(list(H = H))
}

# Write the parameters for the module to the file indicated.
#
writeParameters <- function(filename, input.ds, k.init, k.final, num.clusterings, maxniter, error.function, rseed, directory, stopconv, stopfreq, non.interactive.run, doc.string) {
	input.ds
	time.string <- as.character(as.POSIXlt(Sys.time(),"GMT"))
	write(paste("Run of NMF on ", time.string), file=filename)
	write(paste("input.ds =", input.ds, sep=" "), file=filename, append=T) 
	write(paste("k.init = ", k.init, sep=" "), file=filename, append=T) 
	write(paste("k.final =", k.final, sep=" "), file=filename, append=T) 
	write(paste("num.clusterings =", num.clusterings, sep=" "), file=filename, append=T) 
	write(paste("maxniter =", maxniter, sep=" "), file=filename, append=T) 
	write(paste("error.function =", error.function, sep=" "), file=filename, append=T) 
	write(paste("rseed =", rseed, sep=" "), file=filename, append=T) 
	write(paste("directory =", directory, sep=" "), file=filename, append=T) 
	write(paste("stopconv =", stopconv, sep=" "), file=filename, append=T) 
	write(paste("stopfreq =", stopfreq, sep=" "), file=filename, append=T)
	write(paste("non.interctive.run =", non.interactive.run, sep=" "), file=filename, append=T) 
	write(paste("doc.string =", doc.string, sep=" "), file=filename, append=T) 
}





