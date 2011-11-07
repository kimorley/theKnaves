# Calculation of LRR and BAF
# Function reads in list of dataframes (one for each SNP) and then generates canonical clusters from 
# reference data and calculates LRR and BAF.
# This function returns the data as a list, with or without the input data, of SNP-objects.
# To save individual files use saveLRRBAF command.
calcLRRBAF <- function(data, clusterINTU, clusterGTU, minClusterSize=5, keepOriginal=FALSE){
	# Checks
	if (missing(data)){
		stop("Must supply 'data' argument")
	}
	if (missing(clusterINTU)){
		stop("Must supply location of intensity data for canonical samples.")
	}
	if (missing(clusterGTU)){
		stop("Must supply location of genotype call and confidence data for canonical samples.")
	}
	# Setup and execution
	for (i in 1:length(data)){	# Loop through each SNP
		target <- subset(data[[i]])
		# Generate canonical cluster for SNP
		canonicalData <- convertRawSnp(SNP=names(data)[i],intuFile=clusterINTU,gtuFile=clusterGTU)
		canonicalCluster <- generateClusters(canonicalData)
		clusterData <- canonicalCluster[[1]]
		# Call LRR and BAF for samples
		if (min(clusterData$summary[2,]) < minClusterSize){	# If we do not have enough samples, return nothing for this SNP
			target <- cbind(target,rPred="NA",lrr="NA",baf="NA")
		}else{
			target <- cbind(target,rPred=predict(clusterData$model,target))
			use <- as.character(target$call) %in% names(clusterData$summary[1,])	# Which ones have genotypes seen in the canonical clusters
			target$rPred <- as.numeric(ifelse(use==T,target$rPred,"NA"))	# Set to missing samples with genotypes not seen in canonical clusters (will catch NN as these excluded from clusters)
			target <- cbind(target,lrr=log2(target$r/target$rPred))				
			if (length(clusterData$summary[1,])==3){
				# Check order of clusters
				if (clusterData$summary[1,which(substr(names(clusterData$summary[1,]),1,1)!=substr(names(clusterData$summary[1,]),2,2))] == clusterData$summary[1,2]  && clusterData$summary[1,1] == min(clusterData$summary[1,])){
					# Calculate BAF
					baf <- target$t
					baf <- ifelse(baf <= clusterData$summary[1,1], 0, baf)
					baf <- ifelse(baf >= clusterData$summary[1,3], 1, baf)
					baf <- ifelse(baf != 0 & baf <= clusterData$summary[1,2], 0.5*((baf - clusterData$summary[1,1])/(clusterData$summary[1,2] - clusterData$summary[1,1])), baf)
					baf <- ifelse(baf != 1 & baf >= clusterData$summary[1,2], 0.5+0.5*((baf - clusterData$summary[1,2])/(clusterData$summary[1,3] - clusterData$summary[1,2])), baf)
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
					target <- cbind(target,baf=baf)						
				}else if (clusterData$summary[1,which(substr(names(clusterData$summary[1,]),1,1)!=substr(names(clusterData$summary[1,]),2,2))] == clusterData$summary[1,2]){ # Only AA and AB so cannot calculate BAF > 0.5
					# Calculate BAF
					baf <- target$t
					baf <- ifelse(baf <= clusterData$summary[1,1], 0, baf)
					baf <- ifelse(baf != 0 & baf <= clusterData$summary[1,2], 0.5*((baf - clusterData$summary[1,1])/(clusterData$summary[1,2] - clusterData$summary[1,1])), baf)
					baf <- ifelse(baf > clusterData$summary[1,2], NA, baf)
					target <- cbind(target,baf=baf)
				}else if (clusterData$summary[1,which(substr(names(clusterData$summary[1,]),1,1)!=substr(names(clusterData$summary[1,]),2,2))] == clusterData$summary[1,1]){ # Only AB and BB so cannot calculate BAF < 0.5
					baf <- target$t
					baf <- ifelse(baf >= clusterData$summary[1,2], 1, baf)
					baf <- ifelse(baf != 1 & baf >= clusterData$summary[1,2], 0.5+0.5*((baf - clusterData$summary[1,1])/(clusterData$summary[1,2] - clusterData$summary[1,1])), baf)
					baf <- ifelse(baf < clusterData$summary[1,1], NA, baf)
					target <- cbind(target,baf=baf)
				}else{
					target <- cbind(target,baf=NA)
				}
			}else{
				target <- cbind(target,baf=NA)
			}
		}
		if (keepOriginal){
			data[[i]] <- target
		}else{
			data[[i]] <- subset(target,select=c(lrr,baf))
		}
	}
	return(data)
}		
