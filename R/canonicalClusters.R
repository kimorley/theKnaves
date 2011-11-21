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
	alleles <- unlist(INTall$map$Alleles)
	if (is.na(sum(as.numeric(unlist(INT))))){
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
					gsCluster <- clusterSnp(results, mySnp, alleles, DIR=outDIR, returnCluster=TRUE)
					output <- list(raw=results,cluster=gsCluster)
					return(output)
				}else if (returnInfo){
					gsCluster <- clusterSnp(results, mySnp, alleles, DIR=outDIR, returnInfo=TRUE)
				}else{
					gsCluster <- clusterSnp(results, mySnp, alleles, DIR=outDIR)
				}
			}else{
				stop(paste("Different samples in intensity and genotype files for SNP",mySnp,sep=" "))
			}
		}else{
			stop(paste("IDs for A allele data and B allele data do not match or have different order for SNP", mySnp,sep=" "))
		}
	}	
}

clusterSnp <- function(data, mySnp, alleles, DIR=getwd(), confThresh = 0.99, clusterStat = "median", saveCluster=TRUE, returnCluster=FALSE, returnInfo=TRUE){
	# Check
	if (missing(data)){
		stop("Must supply 'data' argument.")
	}
	if (missing(mySnp)){
		stop("Must supply 'mySnp' argument.")
	}
	if (missing(alleles)){
		stop("Must supply 'alleles' argument.")
	}
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
	if (is.null(data)){
		clusterData <- list(model=NA, summary=NA, alleles=alleles)
		if (saveCluster){ save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep="")) }
		if (returnCluster){ return(clusterData) }
		if (returnInfo){ return(cbind(SNP=mySnp, FLAG="N", N_SAMPLES=0)) }
	}else{
		cc <- data[which( as.numeric(as.character(data$conf)) >= confThresh & as.character(data$call) != "NN"),]	# Complete cases with high confidence
		subsetGenotypes <- function(genotype, data){
			xx <- subset(data, call==genotype)
			if (length(xx[,1]) >= 5){	# Remove genotypes with less than five samples in the cluster
				return(xx)
			}
		}
		ccList <- lapply(unique(as.character(cc$call)), subsetGenotypes, data=cc)
		names(ccList) <- unique(as.character(cc$call))
		cc <- do.call(rbind,ccList)
		if (length(cc[,1])==0){
			print(paste("WARNING: No data for",mySnp))
			clusterData <- list(summary=NA, model=NA, alleles=alleles)
			if (saveCluster){ save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""), compression_level=9) }
			if (returnCluster){ return(clusterData) }
			if (returnInfo){ return(cbind(SNP=mySnp, FLAG="P", N_SAMPLES=0)) }
		}else{
			summary <- data.frame(rbind(tapply(cc$t,as.character(cc$call),clusterStat),table(as.character(cc$call))))
			if (length(summary)==3){	# All three genotypes have adequate data
				aa <- ccList[[which(names(ccList)==paste(substr(alleles,1,1),substr(alleles,1,1),sep=""))]]
				ab <- ccList[[which(names(ccList)==paste(substr(alleles,1,1),substr(alleles,2,2),sep=""))]]
				bb <- ccList[[which(names(ccList)==paste(substr(alleles,2,2),substr(alleles,2,2),sep=""))]]
				if ( median(ab$t) < median(aa$t) || median(ab$t) > median(bb$t) ){
					# Theta value for heterozygote is smaller than AA homozygote or larger than BB homozygote so we don't want to call from this
					clusterData <- list(model=NA, summary=NA, alleles=alleles)
					if (saveCluster){ save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""), compression_level=9) }
					if (returnCluster){ return(clusterData) }
					if (returnInfo){ return(cbind(SNP=mySnp, FLAG="P", N_SAMPLES=0)) }
				}else{
					# Ellipses for BAF calculation
					aaBounds <- predict(ellipsoidhull(as.matrix(cbind(aa$t,aa$r))))
					abBounds <- predict(ellipsoidhull(as.matrix(cbind(ab$t,ab$r))))
					bbBounds <- predict(ellipsoidhull(as.matrix(cbind(bb$t,bb$r))))
					# Spline linear model for LRR calculation
					temp <- lm(r ~ ns(t, knots=median(ab$t), Boundary.knots=c(median(aa$t),median(bb$t))), data=cc, model=FALSE)
					model <- temp[match(c("coefficients","rank","qr","terms"),names(temp))]
					class(model) <- class(temp)
					# Combine info and save/return
					clusterData <- list(aaBounds=aaBounds, abBounds=abBounds, bbBounds=bbBounds, summary=summary, model=model, alleles=alleles)
					if (saveCluster){ save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""), compression_level=9) }
					if (returnCluster){ return(clusterData) }
					if (returnInfo){ return(cbind(SNP=mySnp, FLAG="P", N_SAMPLES=sum(summary[2,]))) }
				}
			}else if (length(summary)==2){	# Only two genotypes have adequate data
				a <- ccList[[which(names(ccList)==names(summary)[1])]]
				b <- ccList[[which(names(ccList)==names(summary)[2])]]
				# Ellipses for BAF calculation
				aBounds <- predict(ellipsoidhull(as.matrix(cbind(a$t,a$r))))
				bBounds <- predict(ellipsoidhull(as.matrix(cbind(b$t,b$r))))
				# Spline linear model for LRR calculation
				temp <- lm(r ~ ns(t, Boundary.knots=c(median(a$t),median(b$t))), data=cc, model=FALSE)
				model <- temp[match(c("coefficients","rank","qr","terms"),names(temp))]
				class(model) <- class(temp)
				clusterData <- list(aBounds=aBounds, bBounds=bBounds, summary=summary, model=model, alleles=alleles)
				# Combine info and save/return
				if (saveCluster){ save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""), compression_level=9) }
				if (returnCluster){ return(clusterData) }
				if (returnInfo){ return(cbind(SNP=mySnp, FLAG="P", N_SAMPLES=sum(summary[2,]))) }
			}else if (length(summary)==1){	# Only one genotype has adequate data
				a <- ccList[[which(names(ccList)==names(summary)[1])]]
				# Ellipses for BAF calculation
				aBounds <- predict(ellipsoidhull(as.matrix(cbind(a$t,a$r))))
				# Spline linear model for LRR calculation
				temp <- lm(r ~ ns(t, Boundary.knots=c(min(a$t),max(a$t))), data=cc, model=FALSE)
				model <- temp[match(c("coefficients","rank","qr","terms"),names(temp))]
				class(model) <- class(temp)
				# Combine info and save/return
				clusterData <- list(aBounds=aBounds, summary=summary, model=model, alleles=alleles)
				if (saveCluster){ save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""), compression_level=9) }
				if (returnCluster){ return(clusterData) }
				if (returnInfo){ return(cbind(SNP=mySnp, FLAG="P", N_SAMPLES=sum(summary[2,]))) }
			}else{
				print(paste("WARNING: No data for",mySnp))
				clusterData <- list(summary=NA, model=NA, alleles=alleles)
				if (saveCluster){ save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""), compression_level=9) }
				if (returnCluster){ return(clusterData) }
				if (returnInfo){ return(cbind(SNP=mySnp, FLAG="P", N_SAMPLES=0)) }
			}
		}
	}
}
