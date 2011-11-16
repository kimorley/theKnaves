# Calculation of LRR and BAF
# This requires:
#	[a] GAPI gtu file (already contained in data.frame)
#	[b] Path of matching GAPI intu file
#	[c] Location of cluster files (plus current chromosome being processed)
# This function returns the data as a list, with or without the input data, of SNP-objects.
# To save individual files use saveLRRBAF command.
#-----------------------------------------------------------------------------------------------------

callLRRBAF <- function(gapiGTU, gapiINTU, clusterPATH, CHR, minClusterSize=1, keepOriginal=FALSE, splitIDs=TRUE){
	# Setup and execution
	snp <- names(eval(as.list(sys.call(-1))[[2]]))[as.numeric(gsub("[^0-9]", "", deparse(substitute(gapiGtu))))]
	mySnp <- unlist(strsplit(snp,"X"))[length(unlist(strsplit(snp,"X")))]
	GTU <- gapiGTU
	INT <- readSnpIntu(mySnp, gapiINTU)
	INT <- INT$intu
	if (is.na(INT)){
		callData <- calculations(INT, mySnp, clusterPATH, minClusterSize)	# This will return NA in this instance
		return(callData)
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
				results <- cbind(data.frame(r,t), GTU)	# This is the input data - R, theta, genotype call, confidence
				callData <- calculations(results, mySnp, clusterPATH, CHR, minClusterSize)	# This is the LRR and BAF called from the canonical clusters
				return(callData)
			}else{
				stop("Different samples in intensity and genotype files.")
			}
		}else{
			stop("IDs for A allele data and B allele data do not match or have different order.")
		}
	}
}
	
calculations <- function(target, SNP, clusterPATH, CHR, minClusterSize, keepOriginal=FALSE){
	if (is.na(target)){
		return(NA)
	}else{
		# Load canonical cluster for SNP
		if (file.exists(paste(clusterPATH,"/chr",CHR,"/clusterFile-",SNP,".RData",sep=""))){
			loadCluster <- try(load(paste(clusterPATH,"/chr",CHR,"/clusterFile-",SNP,".RData",sep="")))
			if (inherits(loadCluster, 'try-error')){
				print(paste("Canonical cluster file for SNP ",SNP," does not exist in specified location.",sep=""))
				return(clusterData <- list(model=0,summary=0))
			}
		}else{
			print(paste("Canonical cluster file for SNP ",SNP," does not exist in specified location.",sep=""))
			return(clusterData <- list(model=0,summary=0))
		}
		# Call LRR and BAF for samples
		if (sum(clusterData$summary)==0){
			target <- cbind(target,rPred="NA",lrr="NA",baf="NA")
		}else if (min(clusterData$summary[2,]) < minClusterSize){	# If we do not have enough samples, return nothing for this SNP
			target <- cbind(target,rPred="NA",lrr="NA",baf="NA")
		}else{
			target <- cbind(target,rPred=predict(clusterData$model,target))
			use <- as.character(target$call) %in% names(clusterData$summary[1,])	# Which ones have genotypes seen in the canonical clusters
			target$rPred <- as.numeric(ifelse(use==T,target$rPred,"NA"))	# Set to missing samples with genotypes not seen in canonical clusters (will catch NN as these excluded from clusters)
			target <- cbind(target,lrr=log2(target$r/target$rPred))		# Generate log2R			
			if (length(clusterData$summary[1,])==3){
				# Check order of clusters
				if (clusterData$summary[1,which(substr(names(clusterData$summary[1,]),1,1)!=substr(names(clusterData$summary[1,]),2,2))] == clusterData$summary[1,2]  && clusterData$summary[1,1] == min(clusterData$summary[1,])){
					# Calculate BAF
					baf <- target$t
					baf <- ifelse(baf <= clusterData$summary[1,1], 0, baf)
					baf <- ifelse(baf >= clusterData$summary[1,3], 1, baf)
					baf <- ifelse(baf != 0 & baf <= clusterData$summary[1,2], 0.5*((baf - clusterData$summary[1,1])/(clusterData$summary[1,2] - clusterData$summary[1,1])), baf)
					baf <- ifelse(baf != 1 & baf >= clusterData$summary[1,2], 0.5+0.5*((baf - clusterData$summary[1,2])/(clusterData$summary[1,3] - clusterData$summary[1,2])), baf)
					baf <- ifelse(baf > 1, 1, baf)
					baf <- ifelse(baf < 0, 0, baf)
					target <- cbind(target,baf=baf)
				}else{
					target <- cbind(target,baf=NA)
				}
			}else if (length(clusterData$summary[1,])==2){
				# Check clusters
				if (sum(substr(names(clusterData$summary[1,]),1,1)!=substr(names(clusterData$summary[1,]),2,2))==0){ # No heterozygote so can only calculate BAF==0/1
					# Calculate BAF
					baf <- target$t
					baf <- ifelse(baf <= clusterData$summary[1,1], 0, baf)
					baf <- ifelse(baf >= clusterData$summary[1,2], 1, baf)
					baf <- ifelse(baf != 0 || baf != 1, NA, baf)
					baf <- ifelse(baf > 1, 1, baf)
					baf <- ifelse(baf < 0, 0, baf)
					target <- cbind(target,baf=baf)						
				}else if (clusterData$summary[1,which(substr(names(clusterData$summary[1,]),1,1)!=substr(names(clusterData$summary[1,]),2,2))] == clusterData$summary[1,2]){ # Only AA and AB so cannot calculate BAF > 0.5
					# Calculate BAF
					baf <- target$t
					baf <- ifelse(baf <= clusterData$summary[1,1], 0, baf)
					baf <- ifelse(baf != 0 & baf <= clusterData$summary[1,2], 0.5*((baf - clusterData$summary[1,1])/(clusterData$summary[1,2] - clusterData$summary[1,1])), baf)
					baf <- ifelse(baf > clusterData$summary[1,2], NA, baf)
					baf <- ifelse(baf > 1, 1, baf)
					baf <- ifelse(baf < 0, 0, baf)
					target <- cbind(target,baf=baf)
				}else if (clusterData$summary[1,which(substr(names(clusterData$summary[1,]),1,1)!=substr(names(clusterData$summary[1,]),2,2))] == clusterData$summary[1,1]){ # Only AB and BB so cannot calculate BAF < 0.5
					baf <- target$t
					baf <- ifelse(baf >= clusterData$summary[1,2], 1, baf)
					baf <- ifelse(baf != 1 & baf >= clusterData$summary[1,2], 0.5+0.5*((baf - clusterData$summary[1,1])/(clusterData$summary[1,2] - clusterData$summary[1,1])), baf)
					baf <- ifelse(baf < clusterData$summary[1,1], NA, baf)
					baf <- ifelse(baf > 1, 1, baf)
					baf <- ifelse(baf < 0, 0, baf)
					target <- cbind(target,baf=baf)
				}else{
					target <- cbind(target,baf=NA)
				}
			}else{
				target <- cbind(target,baf=NA)
			}
		}
		if (keepOriginal){
			return(target)
		}else{
			return(subset(target,select=c(lrr,baf)))
		}
	}
}

