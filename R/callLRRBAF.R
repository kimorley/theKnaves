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
			callData <- calculations(results, mySnp, clusterPATH, CHR, minClusterSize, keepOriginal)	# This is the LRR and BAF called from the canonical clusters
			return(callData)
		}else{
			stop("Different samples in intensity and genotype files.")
		}
	}else{
		stop("IDs for A allele data and B allele data do not match or have different order.")
	}
}
	
calculations <- function(target, SNP, clusterPATH, CHR, minClusterSize, keepOriginal=FALSE){
	print(SNP)
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
			target <- cbind(target,rPred=NA,lrr=NA,baf=NA)
		}else if (min(clusterData$summary[2,]) < minClusterSize){	# If we do not have enough samples, return nothing for this SNP
			target <- cbind(target,rPred=NA,lrr=NA,baf=NA)
		}else{
			# LRR calculations
			target <- cbind(target,rPred=predict(clusterData$model,target))
			use <- as.character(target$call) %in% names(data.frame(clusterData$summary))	# Which ones have genotypes seen in the canonical clusters
			target$rPred <- as.numeric(ifelse(use==T,target$rPred,NA))			# Set to missing samples with genotypes not seen in canonical clusters (will catch NN as these excluded from clusters)
			target <- cbind(target,lrr=log2(target$r/target$rPred))					# Generate log2R	
			# BAF calculations
			bafTest <- function(target){
				if (is.na(target[3])){	# If no genotype call made, return NA
					return(NA)
				}else if (unlist(target[3]) %in% names(clusterData$summary)==FALSE){	# Have we seen this genotype in the canonical clusters?	
					return(NA)
				}else if (!is.na(sum(clusterData[[which(names(clusterData)==target[3])]]))){	# If we did not generate an ellipse for the genotype seen in the target sample
					bounds <- clusterData[[which(names(clusterData)==target[3])]]
					if (point.in.polygon(target[2], target[1], bounds[,1], bounds[,2])){	# Is target sample point within ellipse?
						return(clusterData$summary[3,names(clusterData$summary)==target[3]])
					}else{
						if (clusterData$summary[3,names(clusterData$summary)==target[3]]==0 && as.numeric(target[2]) < as.numeric(clusterData$summary[1,names(clusterData$summary)==target[3]])){ # AA with negative theta
							return(0)
						}else if (clusterData$summary[3,names(clusterData$summary)==target[3]]==1 && as.numeric(target[2]) >= as.numeric(clusterData$summary[1,names(clusterData$summary)==target[3]])){ # BB with larget positive theta
							return(1)
						}else{
							theta <- clusterData$summary[1,]
							names(theta) <- names(clusterData$summary)
							if (length(clusterData$summary)==3){
								if (as.numeric(target[2]) < theta[which(names(theta)==clusterData$alleles)]){	# Positive theta value in space between AA and AB
									value <- 0.5*((as.numeric(target[2]) - theta[which(theta==min(theta))])/(theta[which(names(theta)==clusterData$alleles)] - theta[which(theta==min(theta))]))
									return( ifelse(value < 0, 0, value) )	# Return interpolated value
								}else if (as.numeric(target[2]) >= theta[which(names(theta)==clusterData$alleles)]){ # Positive theta value in space between AB and BB
									value <- 0.5+0.5*((as.numeric(target[2]) - theta[which(names(theta)==clusterData$alleles)])/(theta[which(theta==max(theta))] - theta[which(names(theta)==clusterData$alleles)]))
									return( ifelse(value > 1, 1, value) )
								}else{
									return(NA)
								} 
							}else if (length(clusterData$summary)==2){
								if (theta[which(names(theta)==clusterData$alleles)] == max(theta)){	# Can only evaluate between AA and AB
									if (target[2] <= theta[which(names(theta)==clusterData$alleles)]){	# Point falls in valid space
										value <- 0.5*( ( as.numeric(target[2]) - theta[which(names(theta)!=clusterData$alleles)] ) / ( theta[which(names(theta)==clusterData$alleles)] - theta[which(names(theta)!=clusterData$alleles)] ) )
										return( ifelse(value < 0, 0, value) )	# Return interpolated value	
									}else{
										return(NA)
									}
								}else if (theta[which(names(theta)==clusterData$alleles)] == min(theta)){ # Can only evaluate between AB and BB
									if (target[2] >= theta[which(names(theta)==clusterData$alleles)]){
										value <- 0.5 + 0.5*( ( as.numeric(target[2]) - theta[which(names(theta)==clusterData$alleles)] )/( theta[which(names(theta)!=clusterData$alleles)] - theta[which(names(theta)==clusterData$alleles)] ) )
										return( ifelse(value > 1, 1, value) )
									}else{
										return(NA)
									}
								}else{
									return(NA)
								}
							}else{	# If cluster only has one cloud, can't interpolate
								return(NA)
							}
						}
					}
				}else{	# Genotype doesn't exist in canonical cluster so cannot return value
					return(NA)
				}
			}
			baf <- apply(target, 1, bafTest)
			target <- cbind(target,baf=unlist(baf))
			if (keepOriginal){
				return(target)
			}else{
				return(subset(target,select=c(lrr,baf)))
			}
		}
	}
}

