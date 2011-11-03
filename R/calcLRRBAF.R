calcLRRBAF <- function(data, CHR, clusterDIR=getwd(), minClusterSize=5){
	# Checks
	if (missing(data)){
		stop("Must supply 'data' argument")
	}
	if (missing(CHR)){
		stop("Must supply 'CHR' argument")
	}
	# Setup and execution
	for (i in 1:length(data)){
		target <- subset(data[[i]])
		# Check whether canonical cluster file exists for this SNP and if so, load the model and summary list
		if (file.exists(paste(clusterDIR,"/clusterFile-",names(data)[i],".gzip",sep=""))){
			loadCluster <- try(load(paste(clusterDIR,"/clusterFile-",names(data)[i],".gzip",sep="")))
			if (inherits(loadCluster, 'try-error')){
				warning("Canonical cluster file for SNP ",names(data[i])," does not exist in specified location.")
				summary <- 0
				model <- 0
				clusterData <- list(model=model,summary=summary)
			}
		}else{
			warning("Canonical cluster file for SNP ",names(data[i])," does not exist in specified location.")
			summary <- 0
			model <- 0
			clusterData <- list(model=model,summary=summary)
		}
		# Call LRR and BAF for samples
		if (sum(clusterData$summary)==0){
			warning(paste("No data in reference clusters; cannot generate cluster for ",names(data)[i],sep=""))	
			target <- cbind(target,rPred="NA",baf="NA")
		}else if (min(clusterData$summary[2,]) < minClusterSize){
			target <- cbind(target,rPred="NA",baf="NA")
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
		data[[i]] <- target
	}
	return(data)
}		
