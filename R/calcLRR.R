calcLRR <- function(data, cluster="scottishFounders"){
	stop("This function not active yet!")
	# We use linear interpolation/extrapolation to determine R-expected
	# If samples have a genotype call that is not seen in the canonical cluster sample, they are given a missing value, not extrapolated
	# Check
	# Functions
	if (sum(names(data) %in% names(scottishFounders$model))==length(scottishFounders$model)){
		for (i in 1:length(data)){
			newX <- subset(data[[i]], select=c(theta))				# Theta values for target samples
			model <- scottishFounders$model[[which(names(scottishFounders$model) == names(data)[i])]]	# Canonical cluster model for this SNP
			sum <- scottishFounders$sum[[which(names(scottishFounders$sum) == names(data)[i])]]	# Canonical cluster summary for this SNP
			if (model==0 || sum==0){
				lrr <- rep("NA",length(newX[,1]))	# If we did not generate a cluster, return NA for all samples
				data[[i]] <- cbind(data[[i]], lrr=lrr)
			}else{
				newY <- predict(model, newX, na.action = na.exclude)	# Predicted R values from model
				lrr <- log2(data[[i]]$R/newY)							# Interpolated/extrapolated value for all samples
				use <- as.character(data[[i]]$call) %in% names(scottishFounders$sum[[which(names(scottishFounders$model) == names(data)[i])]])	# Which ones have genotypes seen in the canonical clusters
				outlrr <- ifelse(use=T,lrr,"NA")						# Set to missing those samples with genotypes not seen in canonical clusters (this will also catch NN as these were excluded from clusters)
				data[[i]] <- cbind(data[[i]], lrr=lrr)
			}
		}
	}else{
		stop("Input data is different probe-set to canonical clusters.")
	}
	return(data)
}
