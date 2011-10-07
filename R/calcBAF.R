calcBAF <- function(data, cluster="scottishFounders"){
	stop("This function not active yet!")
	# We use linear interpolation to determine BAF but this is truncated at 0/1 so we are not using a full linear model
	# If samples have a genotype call that is not seen in the canonical cluster sample, they are given a missing value, not extrapolated
	# Check
	# Functions
	if (sum(names(data) %in% names(scottishFounders$model))==length(scottishFounders$model)){
		
		for (i in 1:length(data)){
			newX <- subset(data[[i]], select=c(theta))					# Theta values for target samples
			sum <- scottishFounders$sum[[which(names(scottishFounders$sum) == names(data)[i])]]	# Canonical cluster summary for this SNP
			if (sum==0){									# If we did not generate a canonical cluster, return NA for all samples
				baf <- rep("NA",length(newX[,1]))							
				data[[i]] <- cbind(data[[i]], baf=baf)
			}else{
				calculate <- function(x, summary){
					if (length(summary)==3){	
						if (x >=max(summary[which(substr(names(summary),1,1)==substr(names(summary),2,2))]) ){
							# If theta greater than the BB mean, BAF = 1
							return(1)
						}else if (x <=min(summary[which(substr(names(summary),1,1)==substr(names(summary),2,2))])){
							# If theta less than the AA mean, BAF = 0
							return(0)
						}else if (x <= summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))]){
							# Interpolate between AA and AB
							return( 0.5*(x-min(summary[which(substr(names(summary),1,1)==substr(names(summary),2,2))]))/(summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))]-min(summary[which(substr(names(summary),1,1)==substr(names(summary),2,2))])) )						
						}else if (x >= summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))]){
							# Interpolate between AB and BB
							return( 0.5+0.5(x-summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))])/(max(summary[which(substr(names(summary),1,1)==substr(names(summary),2,2))]) - summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))]) )
						}else{
							return("NA")
						}
					}else if (length(summary)==2){
						# Check if test genotype = cluster genotype
						# If so interpolate between genotypes, otherwise missing
						if (length(summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))])==0){ # No heterozygote
							# Can score BAF for AA and BB with theta <= or >= mean/median of cluster
							# All others missing
						}else{	# Missing one homozygote
							# Can score BAF outside homozygote group and between homozygote and AB mean/median
							# All others missing
						}
					}else{
						return("NA")
					}
				
				}
				?baf <- mapply(calculate,newX,sum)				
				use <- as.character(data[[i]]$call) %in% names(scottishFounders$sum[[which(names(scottishFounders$model) == names(data)[i])]])	# Which ones have genotypes seen in the canonical clusters
				outlrr <- ifelse(use=T,lrr,"NA")						# Set to missing those samples with genotypes not seen in canonical clusters (this will also catch NN as these were excluded from clusters)
				data[[i]] <- cbind(data[[i]], baf=baf)
				
			}
		}
		
	}else{
		stop("Input data is different probe-set to canonical clusters.")
	}
	return(data)
}
