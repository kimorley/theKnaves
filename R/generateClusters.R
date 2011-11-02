generateClusters <- function(data, confThresh = 0.99, clusterStat = "median"){
	# Check
	if (missing(data)){
		stop("Must supply 'data' argument")
	}
	if (!is.list(data)){
		stop("Input data must be in list format")
	}
	if (clusterStat=="median"){
		print("Generating median-based clusters")
	}else if (clusterStat=="mean"){
		print("Generating mean-based clusters")
	}else{
		stop("Invalid summary statistic for clusters.  Select either 'median' or 'mean'.")
	}
	if (confThresh < 0 | confThresh > 1){
		stop("Invalid value for confidence threshold.  Please supply value between 0 and 1.")
	}
	# Functions
	cluster <- function(x){
		newX <- x[which( as.numeric(as.character(x$conf)) >= confThresh & as.character(x$call) != "NN"),]
		if (length(newX[,1])<=length(as.numeric(x$conf))/2){	# This is arbitary; basically just saying if we have to remove 50% of sample the cluster won't be used.  Need to evaluate this.
			warning(gettextf("Not enough data to generate canonical clusters"), domain = NA)
			summary <- 0
			return(summary)
		}else{
			summary <- tapply(newX$t,as.character(newX$call),clusterStat)
			# Check whether the results are sensible
			if (length(summary)==3){
				if ( summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))] > max(summary[which(substr(names(summary),1,1)==substr(names(summary),2,2))]) || summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))] < min(summary[which(substr(names(summary),1,1)==substr(names(summary),2,2))]) ){
					# Theta value for heterozygote is smaller than AA homozygote or larger than BB homozygote so we don't want to call from this
					summary <- 0
					model <- 0
					clusterData <- list(model=model,summary=summary)
					return(clusterData)
				}else{
					model <- lm(r ~ ns(t, knots=summary[2], Boundary.knots=c(summary[1],summary[3])), data=newX)
					summary <- rbind(tapply(newX$t,as.character(newX$call),clusterStat),table(as.character(newX$call)))
					clusterData <- list(model=model,summary=summary)
					return(clusterData)
				}
			}else if (length(summary)==2){
				model <- lm(r ~ ns(t, Boundary.knots=c(summary[1],summary[2])), data=newX)
				summary <- rbind(tapply(newX$t,as.character(newX$call),clusterStat),table(as.character(newX$call)))
				clusterData <- list(model=model,summary=summary)
				return(clusterData)
			}else if (length(summary)==1){
				model <- lm(r ~ ns(t, Boundary.knots=c(min(t),max(t))), data=newX)
				summary <- rbind(tapply(newX$t,as.character(newX$call),clusterStat),table(as.character(newX$call)))
				clusterData <- list(model=model,summary=summary)
				return(clusterData)
			}else{
				warning("Number of genotypes should be between 1 and 3! Cannot generate canonical clusters.")
			}
		}
	}
	# Apply
	output <- lapply(data,cluster)
	return(output)
}