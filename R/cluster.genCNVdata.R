cluster.genCNVdata <- function(gapiGTU, gapiINTU, clusterPATH, CHR, minClusterSize=1, keepOriginal=FALSE, splitIDs=TRUE){
	# Setup and execution
	snp <- names(eval(as.list(sys.call(-1))[[2]]))[as.numeric(gsub("[^0-9]", "", deparse(substitute(gapiGTU))))]
	SNP <- unlist(strsplit(snp,"X"))[length(unlist(strsplit(snp,"X")))]
	GTU <- gapiGTU
	INT <- readSnpIntu(SNP, gapiINTU)
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
		# r <- sqrt(as.numeric(as.character(A))^2 + as.numeric(as.character(B))^2)	# This is the Euclidean distance i.e. the hypotenuse of the triangle for conversion of Cartesian to polar BUT Illumina DOESN'T use this
		r <- as.numeric(as.character(A)) + as.numeric(as.character(B))	# This is the Manhattan distance, which is what Illumina use
		names(r) <- namesA
		names(t) <- namesA
		if (sum(row.names(GTU) == namesA) == length(namesA)){
			target <- cbind(data.frame(r,t), GTU)	# This is the input data - R, theta, genotype call, confidence
			if (sum(is.na(target))==1){
				return(NA)
			}else{
				# Load canonical clusters for SNP
				if (file.exists(paste(clusterPATH,"/",CHR,"/clusterFile-",SNP,".RData",sep=""))){
					loadCluster <- try(load(paste(clusterPATH,"/",CHR,"/clusterFile-",SNP,".RData",sep="")))
					if (inherits(loadCluster, 'try-error')){
						print(paste("Canonical cluster file for SNP ",SNP," not found at ",clusterPATH,"/",CHR,"/clusterFile-",SNP,".RData",sep=""))
						return(cData <- list(m=NULL,s=NULL))
					}
				}else{
					print(paste("Canonical cluster file for SNP ",SNP," not found at ",clusterPATH,"/",CHR,"/clusterFile-",SNP,".RData",sep=""))
					return(cData <- list(m=NULL,s=NULL))
				}	
				# Call LRR and BAF for samples		
				if (sum(is.na(cData$s))==1){
					target <- cbind(target,rPred=NA,lrr=NA,baf=NA)
				}else if (min(cData$s[row.names(cData$s)=="n",]) < minClusterSize){	# If we do not have enough samples, return nothing for this SNP
					target <- cbind(target,rPred=NA,lrr=NA,baf=NA)
				}else{
					# LRR calculations
					target <- cbind(target,rPred=predict(cData$m,target))
					use <- as.character(target$call) %in% names(data.frame(cData$s))	# Which ones have genotypes seen in the canonical clusters
					use <- ifelse(use==F & as.numeric(as.character(target$t)) <= max(cData$s[row.names(cData$ss)=="max",]) & as.numeric(as.character(target$t)) >= min(cData$s[row.names(cData$ss)=="min",]), TRUE, use)
					target$rPred <- as.numeric(ifelse(use==T,target$rPred,NA))								# Set to missing samples with no genotype and theta values outside those seen in canonical samples
					target <- cbind(target,lrr=log2(target$r/target$rPred))									# Generate log2R	
					# BAF calculations
					baf <- function(target){
						theta <- cData$s[row.names(cData$s)=="median",]
						names(theta) <- names(cData$s)
						if (length(cData$s)==3){
							if ( target[2] < cData$s[row.names(cData$s)=="median",cData$s[1,]==0] ){	# Less than AA median
								return(0)
							}else if ( target[2] > cData$s[row.names(cData$s)=="median",cData$s[1,]==1] ){	# Greater than BB median
								return(1)
							}else if (as.numeric(target[2]) < theta[which(names(theta)==cData$a)]){	# Positive theta value in space between AA and AB
								value <- 0.5*((as.numeric(target[2]) - theta[which(theta==min(theta))])/(theta[which(names(theta)==cData$a)] - theta[which(theta==min(theta))]))
								return( ifelse(value < 0, 0, value) )	# Return interpolated value
							}else if (as.numeric(target[2]) >= theta[which(names(theta)==cData$a)]){ # Positive theta value in space between AB and BB
								value <- 0.5+0.5*((as.numeric(target[2]) - theta[which(names(theta)==cData$a)])/(theta[which(theta==max(theta))] - theta[which(names(theta)==cData$a)]))
								return( ifelse(value > 1, 1, value) )
							}else{
								return(NA)
							} 
						}else if (length(cData$s)==2){
							if (cData$a %in% names(cData$s)==FALSE){	# Heterozygote not observed in canonical clusters
								if ( target[2] < cData$s[row.names(cData$s)=="median",cData$s[1,]==0] ){	# Less than AA median
									return(0)
								}else if ( target[2] > cData$s[row.names(cData$s)=="median",cData$s[1,]==1] ){	# Greater than BB median
									return(1)
								}else{
									return(NA)
								}
							}else{
								if (theta[which(names(theta)==cData$a)] == max(theta)){	# Can only evaluate between AA and AB
									if ( target[2] < cData$s[row.names(cData$s)=="median",cData$s[1,]==0] ){	# Less than AA median
										return(0)
									}else if (target[2] <= theta[which(names(theta)==cData$a)]){	# Point falls in valid space
										value <- 0.5*( ( as.numeric(target[2]) - theta[which(names(theta)!=cData$a)] ) / ( theta[which(names(theta)==cData$a)] - theta[which(names(theta)!=cData$a)] ) )
										return( ifelse(value < 0, 0, value) )	# Return interpolated value	
									}else{
										return(NA)
									}
								}else if (theta[which(names(theta)==cData$a)] == min(theta)){ # Can only evaluate between AB and BB
									if ( target[2] > cData$s[row.names(cData$s)=="median",cData$s[1,]==1] ){	# Greater than BB median
										return(1)
									}else if (target[2] >= theta[which(names(theta)==cData$a)]){
										value <- 0.5 + 0.5*( ( as.numeric(target[2]) - theta[which(names(theta)==cData$a)] )/( theta[which(names(theta)!=cData$a)] - theta[which(names(theta)==cData$a)] ) )
										return( ifelse(value > 1, 1, value) )
									}else{
										return(NA)
									}
								}else{
									return(NA)
								}
							}						
						}else{	# If cluster only has one cloud, can't interpolate
							if (target[3] == names(cData$s)){
								if ((cData$s[1,] == 0) && (target[2] < cData$s[row.names(cData$s)=="median",cData$s[1,]==0])){
									return(0)
								}else if ((cData$s[1,] == 1) && (target[2] > cData$s[row.names(cData$s)=="median",cData$s[1,]==1])){
									return(1)
								}else{
									return(NA)
								}
							}else{
								return(NA)
							}
						}
					}
					output <- apply(target, 1, baf)
					target <- cbind(target,baf=unlist(output))
					if (keepOriginal){
						return(target)
					}else{
						return(subset(target,select=c(lrr,baf)))
					}
				}
			}
			return(target)
		}else{
			stop("Different samples in intensity and genotype files.")
		}
	}else{
		stop("IDs for A allele data and B allele data do not match or have different order.")
	}
}
