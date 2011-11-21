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
	snp <- names(eval(as.list(sys.call(-1))[[2]]))[as.numeric(gsub("[^0-9]", "", deparse(substitute(gapiGTU))))]
	mySnp <- unlist(strsplit(snp,"X"))[length(unlist(strsplit(snp,"X")))]
	GTU <- gapiGTU
	INT <- readSnpIntu(mySnp, gapiINTU)
	INT <- INT$intu
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
	
calculations <- function(target, SNP, clusterPATH, CHR, minClusterSize, keepOriginal=FALSE){
	if (sum(is.na(target))==1){
		return(NA)
	}else{
		# Load canonical cluster for SNP
		if (file.exists(paste(clusterPATH,"/chr",CHR,"/clusterFile-",SNP,".RData",sep=""))){
			loadCluster <- try(load(paste(clusterPATH,"/chr",CHR,"/clusterFile-",SNP,".RData",sep="")))
			if (inherits(loadCluster, 'try-error')){
				print(paste("Canonical cluster file for SNP ",SNP," does not exist in specified location.",sep=""))
				return(clusterData <- list(model=NULL,summary=NULL))
			}
		}else{
			print(paste("Canonical cluster file for SNP ",SNP," does not exist in specified location.",sep=""))
			return(clusterData <- list(model=NULL,summary=NULL))
		}
		# Call LRR and BAF for samples
		if (sum(is.na(clusterData$summary))==1){
			target <- cbind(target,rPred="NA",lrr="NA",baf="NA")
		}else if (min(clusterData$summary[2,]) < minClusterSize){	# If we do not have enough samples, return nothing for this SNP
			target <- cbind(target,rPred="NA",lrr="NA",baf="NA")
		}else{
			# LRR calculations
			target <- cbind(target,rPred=predict(clusterData$model,target))
			use <- as.character(target$call) %in% names(data.frame(clusterData$summary))	# Which ones have genotypes seen in the canonical clusters
			target$rPred <- as.numeric(ifelse(use==T,target$rPred,"NA"))			# Set to missing samples with genotypes not seen in canonical clusters (will catch NN as these excluded from clusters)
			target <- cbind(target,lrr=log2(target$r/target$rPred))					# Generate log2R	
			# BAF calculations
			targetTest <- function(target){
				if (is.na(target[3])){	# If no genotype call made, return NA
					return(NA)
				}else{
					if (length(clusterData$summary[1,])==3){	# If there were three genotype groups in the canonical clusters
						if ( point.in.polygon(target[2], target[1], clusterData$aaBounds[,1], clusterData$aaBounds[,2]) ){	# Point within AA cluster
							if (target[3]==paste(substr(clusterData$alleles,1,1),substr(clusterData$alleles,1,1),sep="")){	# Genotype call same as canonical cluster
								return(0)
							}else{
								return(NA)
							}
						}else if ( point.in.polygon(target[2], target[1], clusterData$abBounds[,1], clusterData$abBounds[,2]) ){ # Point within AB cluster
							if (target[3]==paste(substr(clusterData$alleles,1,1),substr(clusterData$alleles,2,2),sep="")){	# Genotype call same as canonical cluster
								return(0.5)
							}else{
								return(NA)
							}
						}else if ( point.in.polygon(target[2], target[1], clusterData$bbBounds[,1], clusterData$bbBounds[,2]) ){ # Point within BB cluster
							if ( target[3]==paste(substr(clusterData$alleles,2,2),substr(clusterData$alleles,2,2),sep="") ){	# Genotype call same as canonical cluster
								return(1)
							}else{
								return(NA)
							}
						}else{
							theta <- as.numeric(clusterData$summary[1,])
							if (target[3]==paste(substr(clusterData$alleles,1,1),substr(clusterData$alleles,1,1),sep="") && target[2] < 0){	# AA with negative theta
								return(0)
							}else if (target[3]==paste(substr(clusterData$alleles,2,2),substr(clusterData$alleles,2,2),sep="") && target[2] > 1){ # BB with large positive theta
								return(1)
							}else if (as.numeric(target[2]) < theta[1]){	# Positive theta value in space between AA and AB
								value <- 0.5*((as.numeric(target[2]) - theta[1])/(theta[2] - theta[1]))
								return( ifelse(value < 0, 0, value) )	# Return interpolated value
							}else if (as.numeric(target[2]) > theta[2]){ # Positive theta value in space between AB and BB
								value <- 0.5+0.5*((as.numeric(target[2]) - theta[2])/(theta[3] -theta[2]))
								return( ifelse(value > 1, 1, value) )
							}else{
								return(NA)
							} 
						} 
					}else if (length(clusterData$summary[1,])==2){	# If there were two genotype groups in the canonical clusters
						if (clusterData$alleles %in% names(clusterData$summary)){	# XX and XY
							XX <- subset(clusterData$summary, select=names(clusterData$summary) != clusterData$alleles)
							XY <- subset(clusterData$summary, select=names(clusterData$summary) == clusterData$alleles)
						 	if ( point.in.polygon(target[2], target[1], clusterData$bBounds[,1], clusterData$bBounds[,2]) ){ # Point within XY cluster
								if (target[3]==names(XY)){	# Genotype call same as canonical cluster
									return(0.5)
								}else{
									return(NA)
								}
							}else if (XX[1,1] < XY[1,1]){	# XX BAF is 0
								if ( point.in.polygon(target[2], target[1], clusterData$aBounds[,1], clusterData$aBounds[,2]) ){	# Point within XX cluster
									if (target[3]==names(XX)){	# Genotype call same as canonical cluster
										return(0)
									}else{
										return(NA)
									}
								}else{
									theta <- as.numeric(clusterData$summary[1,])
									if ((target[3] %in% names(clusterData$summary))==FALSE){	# Ensure we don't call BAF for genotype not seen in canonical clusters
										return(NA)
									}else if (target[3]==names(XX) && target[2] < 0){	# AA with negative theta
										return(0)
									}else if (as.numeric(target[2]) < theta[1]){	# Positive theta value in space between AA and AB
										value <- 0.5*((as.numeric(target[2]) - theta[1])/(theta[2] - theta[1]))
										return( ifelse(value < 0, 0, value) )	# Return interpolated value
									}else{
										return(NA)
									} 
								}
							}else if (XX[1,1] >= XY[1,1]){	# XX BAF is 1
								if ( point.in.polygon(target[2], target[1], clusterData$aBounds[,1], clusterData$aBounds[,2]) ){	# Point within XX cluster
									if (target[3]==names(XX)){	# Genotype call same as canonical cluster
										return(1)
									}else{
										return(NA)
									}
								}else{
									theta <- as.numeric(clusterData$summary[1,])
									if ((target[3] %in% names(clusterData$summary))==FALSE){	# Ensure we don't call BAF for genotype not seen in canonical clusters
										return(NA)
									}else if (target[3]==names(XX) && target[2] > 1){ # BB with large positive theta
										return(1)
									}else if (as.numeric(target[2]) > theta[2]){ # Positive theta value in space between AB and BB
										value <- 0.5+0.5*((as.numeric(target[2]) - theta[2])/(theta[3] -theta[2]))
										return( ifelse(value > 1, 1, value) )
									}else{
										return(NA)
									} 
									
								}
							}else{	# Fail
								return(NA)
							}
						}else{ # XX & YY
							XX <- clusterData$summary[,1]
							YY <- clusterData$summary[,2]
							if (XX[1,1] < YY[1,1]){
								if ( point.in.polygon(target[2], target[1], clusterData$aBounds[,1], clusterData$aBounds[,2]) ){	# Point within XX cluster
									if (target[3]==names(XX)){	# Genotype call same as canonical cluster
										return(0)
									}else{
										return(NA)
									}
								}else if ( point.in.polygon(target[2], target[1], clusterData$bBounds[,1], clusterData$bBounds[,2]) ){ # Point within YY cluster
									if (target[3]==names(YY)){	# Genotype call same as canonical cluster
										return(1)
									}else{
										return(NA)
									}
								}else{
									theta <- as.numeric(clusterData$summary[1,])
									if (target[3]==clusterData$alleles){	# Ensure we don't call BAF for genotype not seen in canonical clusters
										return(NA)
									}else if (target[3]==names(XX) && target[2] < 0){	# XX with negative theta
										return(0)
									}else if (target[3]==names(YY) && target[2] > 1){ # BB with large positive theta
										return(1)
									}else{		# Cannot interpolate as no AB cluster
										return(NA)
									}
								}  
							}else if (YY[1,1] >= XX[1,1]){
								if ( point.in.polygon(target[2], target[1], clusterData$aBounds[,1], clusterData$aBounds[,2]) ){	# Point within XX cluster
									if (target[3]==names(XX)){	# Genotype call same as canonical cluster
										return(1)
									}else{
										return(NA)
									}
								}else if ( point.in.polygon(target[2], target[1], clusterData$bBounds[,1], clusterData$bBounds[,2]) ){ # Point within YY cluster
									if (target[3]==names(YY)){	# Genotype call same as canonical cluster
										return(0)
									}else{
										return(NA)
									}
								}else{
									theta <- as.numeric(clusterData$summary[1,])
									if (target[3]==clusterData$alleles){	# Ensure we don't call BAF for genotype not seen in canonical clusters
										return(NA)
									}else if (target[3]==names(YY) && target[2] < 0){	# XX with negative theta
										return(0)
									}else if (target[3]==names(XX) && target[2] > 1){ # BB with large positive theta
										return(1)
									}else{		# Cannot interpolate as no AB cluster
										return(NA)
									}
								}
							}else{	# Fail
								return(NA)
							}
						}
					}else if (length(clusterData$summary[1,])==1){
						if ( point.in.polygon(target[2], target[1], clusterData$aBounds[,1], clusterData$aBounds[,2]) ){	# Point within XX cluster
							if (target[3]==names(clusterData$summary)){	# Genotype call same as canonical cluster
								if (clusterData$alleles==names(clusterData$summary)){
									return(0.5)
								}else if (names(clusterData$summary)==paste(substr(clusterData$alleles,1,1),substr(clusterData$alleles,1,1),sep="")){
									return(0)
								}else{
									return(1)
								}
							}else{
								return(NA)
							}
						}else{
							return(NA)
						}
					}
				}
			}
			baf <- apply(target, 1, targetTest)
			target <- cbind(target,baf)
		}
		if (keepOriginal){
			return(target)
		}else{
			return(subset(target,select=c(lrr,baf)))
		}
	}
}

