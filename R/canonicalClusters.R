# Functions for generating canonical clusters
# Optimised for minimal output file size, not speed, as this should only be used
# to generate new canonical cluster files and not on a regular basis
# USAGE:
# Requires:
#	[a] GAPI gtu file (already contained in data.frame)
#	[b] Path of matching GAPI intu file
#	[c] Directory for saving cluster files (*.RData)
# 	clusterData <- readRawGtu(<FILENAME>)
# 	lapply(clusterData, convertLines, gapiINTU=<FILENAME>, outDIR=<OUTDIR>)
# N.B. Large memory requirements
#-----------------------------------------------------------------------------------------------------

convertLines <- function(gapiGTU, gapiINTU, outDIR=getwd(), splitIDs=TRUE, returnData=FALSE, returnInfo=TRUE){
	if (returnData == returnInfo){
		stop("Cannot return cluster results and info - choose one option.")
	}
	snp <- names(eval(as.list(sys.call(-1))[[2]]))[as.numeric(gsub("[^0-9]", "", deparse(substitute(gapiGTU))))]
	mySnp <- unlist(strsplit(snp,"X"))[length(unlist(strsplit(snp,"X")))]
	GTU <- gapiGTU
	INTall <- readSnpIntu(mySnp, gapiINTU)
	INT <- INTall$intu
	if (is.null(dim(INT))){
		print(paste("WARNING: No intensity data available for SNP",mySnp,"in",gapiINTU,sep=" "))
		gsCluster <- clusterSnp(x=NA, mySnp, DIR=outDIR)
	}else{
		# Split intensity data into A and B allele sets
		index <- rep(c(0,1),(length(row.names(INT))/2))
		use <- index==0
		A <- INT[use,]
		use <- index==1
		B <- INT[use,]
		namesA <- names(A)
		namesB <- names(B)
		if (splitIDs){	# Take last two elements of ID split on underscore
			namesA <- sapply(names(A),idDouble)
			namesB <- sapply(names(B),idDouble)
		}
		if (sum(namesA==namesB)==length(row.names(INT))/2){
			# Calculate angle (theta) from signal intensities
			t <- atan2(as.numeric(as.character(B)),as.numeric(as.character(A))) / (pi/2) # Make sure values range from 0 to 1
			# Calculate ray (R) from signal intensities
			# r <- sqrt(as.numeric(as.character(A))^2 + as.numeric(as.character(B))^2)	# This is the Euclididean distance i.e. the hypotenuse of the triangle for conversion of Cartesian to polar BUT Illumina DOESN'T use this
			r <- as.numeric(as.character(A)) + as.numeric(as.character(B))	# This is the Manhattan distance, which is what Illumina use
			names(r) <- namesA
			names(t) <- namesA
			if (sum(row.names(GTU) == namesA) == length(namesA)){
				results <- cbind(data.frame(r,t), GTU)
				if (returnData){
					gsCluster <- clusterSnp(results, mySnp, DIR=outDIR, returnCluster=TRUE)
					output <- list(raw=results,cluster=gsCluster)
					return(output)
				}else if (returnInfo){
					gsCluster <- clusterSnp(results, mySnp, DIR=outDIR, returnInfo=TRUE)
				}else{
					gsCluster <- clusterSnp(results, mySnp, DIR=outDIR)
				}
			}else{
				stop(paste("Different samples in intensity and genotype files for SNP",mySnp,sep=" "))
			}
		}else{
			stop(paste("IDs for A allele data and B allele data do not match or have different order for SNP", mySnp,sep=" "))
		}
	}	
}

clusterSnp <- function(x, mySnp, DIR=getwd(), confThresh = 0.99, clusterStat = "median", saveCluster=TRUE, returnCluster=FALSE, returnInfo=TRUE){
	# Check
	if (clusterStat!="median" && clusterStat != "mean"){
		stop("Invalid summary statistic for clusters.  Select either 'median' or 'mean'.")
	}
	if (confThresh < 0 | confThresh > 1){
		stop("Invalid value for confidence threshold.  Please supply value between 0 and 1.")
	}
	if (returnCluster == returnInfo){
		stop("Cannot return cluster results and info - choose one option.")
	}
	# Body
	if (is.null(dim(x))){
		summary <- 0
		model <- 0
		clusterData <- list(model=model,summary=summary)
		if (saveCluster){
			save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""))
		}
		if (returnCluster){
			return(clusterData)
		}
		if (returnInfo){
			return(cbind(SNP=mySnp,FLAG="N"))
		}
	}else{
		newX <- x[which( as.numeric(as.character(x$conf)) >= confThresh & as.character(x$call) != "NN"),]
		summary <- tapply(newX$t,as.character(newX$call),clusterStat)
		# Check whether the results are sensible
		if (length(summary)==3){
			if ( summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))] > max(summary[which(substr(names(summary),1,1)==substr(names(summary),2,2))]) || summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))] < min(summary[which(substr(names(summary),1,1)==substr(names(summary),2,2))]) ){
				# Theta value for heterozygote is smaller than AA homozygote or larger than BB homozygote so we don't want to call from this
				summary <- 0
				model <- 0
				clusterData <- list(model=model,summary=summary)
				if (saveCluster){
					save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""))
				}
				if (returnCluster){
					return(clusterData)
				}
				if (returnInfo){
					return(cbind(SNP=mySnp,FLAG="F"))
				}
			}else{
				temp <- lm(r ~ ns(t, knots=summary[2], Boundary.knots=c(summary[1],summary[3])), data=newX, model=FALSE)
				model <- temp[match(c("coefficients","rank","qr","terms"),names(temp))]
				class(model) <- class(temp)
				summary <- rbind(tapply(newX$t,as.character(newX$call),clusterStat),table(as.character(newX$call)))
				clusterData <- list(model=model,summary=summary)
				if (saveCluster){
					save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""), compression_level=9)
				}
				if (returnCluster){
					return(clusterData)
				}
				if (returnInfo){
					return(cbind(SNP=mySnp,FLAG="P"))
				}
			}
		}else if (length(summary)==2){
			temp <- lm(r ~ ns(t, Boundary.knots=c(summary[1],summary[2])), data=newX, model=FALSE)
			model <- temp[match(c("coefficients","rank","qr","terms"),names(temp))]
			class(model) <- class(temp)
			summary <- rbind(tapply(newX$t,as.character(newX$call),clusterStat),table(as.character(newX$call)))
			clusterData <- list(model=model,summary=summary)
			if (saveCluster){
				save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""), compression_level=9)
			}
			if (returnCluster){
				return(clusterData)
			}
			if (returnInfo){
				return(cbind(SNP=mySnp,FLAG="P"))
			}
		}else if (length(summary)==1){
			temp <- lm(r ~ ns(t, Boundary.knots=c(min(t),max(t))), data=newX, model=FALSE)
			model <- temp[match(c("coefficients","rank","qr","terms"),names(temp))]
			class(model) <- class(temp)
			summary <- rbind(tapply(newX$t,as.character(newX$call),clusterStat),table(as.character(newX$call)))
			clusterData <- list(model=model,summary=summary)
			if (saveCluster){
				save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""), compression_level=9)
			}
			if (returnCluster){
				return(clusterData)
			}
			if (returnInfo){
				return(cbind(SNP=mySnp,FLAG="P"))
			}
		}else{
			print(paste("WARNING: Number of genotypes should be between 1 and 3! Cannot generate canonical clusters for ",mySnp,".",sep=""))
			model <- 0
			summary <- 0
			clusterData <- list(model=model,summary=summary)
			if (saveCluster){
				save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""))
			}
			if (returnCluster){
				return(clusterData)
			}
			if (returnInfo){
				return(cbind(SNP=mySnp,FLAG="N"))
			}
		}
	}
}
