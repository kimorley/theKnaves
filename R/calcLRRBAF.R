calcLRRBAF <- function(data, cluster, minClusterSize=5){
	# Checks
	if (missing(data)){
		stop("Must supply 'data' argument")
	}
	if (missing(cluster)){
		stop("Must supply 'cluster' argument")
	}
	# Setup and execution
	if (sum(names(data) %in% names(cluster))==length(cluster)){
		for (i in 1:length(data)){
			target <- subset(data[[i]])
			ref <- cluster[[which(names(cluster) == names(data)[i])]]	# Canonical cluster data for this SNP
			if (sum(ref$summary)==0){
				warning(paste("No data in reference clusters; cannot generate cluster for ",names(data[i]),sep=""))	
			}else if (min(ref$summary[2,]) < minClusterSize){
				target <- cbind(target,rPred="NA",baf="NA")
			}else{
				target <- cbind(target,rPred=predict(ref$model,target))
				use <- as.character(target$call) %in% names(ref$summary[1,])	# Which ones have genotypes seen in the canonical clusters
				target$rPred <- as.numeric(ifelse(use==T,target$rPred,"NA"))	# Set to missing samples with genotypes not seen in canonical clusters (will catch NN as these excluded from clusters)
				target <- cbind(target,lrr=log2(target$r/target$rPred))				
				if (length(ref$summary[1,])==3){
					# Check order of clusters
					if (ref$summary[1,which(substr(names(ref$summary[1,]),1,1)!=substr(names(ref$summary[1,]),2,2))] == ref$summary[1,2]  && ref$summary[1,1] == min(ref$summary[1,])){
						# Calculate BAF
						baf <- target$t
						baf <- ifelse(baf <= ref$summary[1,1], 0, baf)
						baf <- ifelse(baf >= ref$summary[1,3], 1, baf)
						baf <- ifelse(baf != 0 & baf <= ref$summary[1,2], 0.5*((baf - ref$summary[1,1])/(ref$summary[1,2] - ref$summary[1,1])), baf)
						baf <- ifelse(baf != 1 & baf >= ref$summary[1,2], 0.5+0.5*((baf - ref$summary[1,2])/(ref$summary[1,3] - ref$summary[1,2])), baf)
						target <- cbind(target,baf=baf)
					}else{
						target <- cbind(target,baf=NA)
					}
				}else if (length(ref$summary[1,])==2){
					# Check clusters
					if (sum(substr(names(ref$summary[1,]),1,1)!=substr(names(ref$summary[1,]),2,2))==0){ # No heterozygote so can only calculate BAF==0/1
						# Calculate BAF
						baf <- target$t
						baf <- ifelse(baf <= ref$summary[1,1], 0, baf)
						baf <- ifelse(baf >= ref$summary[1,2], 1, baf)
						baf <- ifelse(baf != 0 || baf != 1, NA, baf)
						target <- cbind(target,baf=baf)						
					}else if (ref$summary[1,which(substr(names(ref$summary[1,]),1,1)!=substr(names(ref$summary[1,]),2,2))] == ref$summary[1,2]){ # Only AA and AB so cannot calculate BAF > 0.5
						# Calculate BAF
						baf <- target$t
						baf <- ifelse(baf <= ref$summary[1,1], 0, baf)
						baf <- ifelse(baf != 0 & baf <= ref$summary[1,2], 0.5*((baf - ref$summary[1,1])/(ref$summary[1,2] - ref$summary[1,1])), baf)
						baf <- ifelse(baf > ref$summary[1,2], NA, baf)
						target <- cbind(target,baf=baf)
					}else if (ref$summary[1,which(substr(names(ref$summary[1,]),1,1)!=substr(names(ref$summary[1,]),2,2))] == ref$summary[1,1]){ # Only AB and BB so cannot calculate BAF < 0.5
						baf <- target$t
						baf <- ifelse(baf >= ref$summary[1,2], 1, baf)
						baf <- ifelse(baf != 1 & baf >= ref$summary[1,2], 0.5+0.5*((baf - ref$summary[1,1])/(ref$summary[1,2] - ref$summary[1,1])), baf)
						baf <- ifelse(baf < ref$summary[1,1], NA, baf)
						target <- cbind(target,baf=baf)
					}else{
						target <- cbind(target,baf=NA)
					}
				}else{
					target <- cbind(target,baf=NA)
				}
			}
			data[[i]] <- target
		}	
	}else{
		stop("Input data is different probe-set to canonical clusters.")
	}
	return(data)
}		
