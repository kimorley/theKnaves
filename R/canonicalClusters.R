canonicalClusters <- function(data, thresh = 0.9, stat = "median"){
	stop("This function not active yet!")
	# Check
	if (missing(data)){
		stop("Must supply 'data' argument")
	}
	if (!is.list(data)){
		stop("Input data must be in list format")
	}
	if (stat=="median"){
		print("Generating median-based clusters")
	}else if (stat=="mean"){
		print("Generating mean-based clusters")
	}else{
		stop("Invalid summary statistic for clusters.  Select either 'median' or 'mean'.")
	}
	if (thresh < 0 | thresh > 1){
		stop("Invalid value for confidence threshold.  Please supply value between 0 and 1.")
	}
	# Functions
	clusterSummary <- function(x){
		newX <- x[which( as.numeric(as.character(x$conf)) >= THRESH & as.character(dataSet$call) != "NN"),]
		if (length(newX[,1])<=length(as.numeric(x$conf))/2){	# This is arbitary; basically just saying if we have to remove 50% of sample the cluster won't be used.  Need to evaluate this.
			warning(gettextf("Not enough data to generate canonical clusters"), domain = NA)
			summary <- 0
			return(summary)
		}else{
			modelData <- data.frame( theta=as.numeric(as.character(newX$theta)), call=as.character(newX$call) )
			summary <- tapply(modelData$theta,modelData$call,stat)
			# Check whether the results are sensible
			if (length(summary)==3){
				if ( summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))] > max(summary[which(substr(names(summary),1,1)==substr(names(summary),2,2))]) || summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))] < min(summary[which(substr(names(summary),1,1)==substr(names(summary),2,2))]) ){
					# Theta value for heterozygote is smaller than AA homozygote or larger than BB homozygote so we don't want to call from this
					summary <- 0
				}
			}
			return(summary)
		}
	}
	clusterModel <- function(x){
		newX <- x[which( as.numeric(as.character(x$conf)) >= THRESH & as.character(dataSet$call) != "NN"),]
		if (length(newX[,1])<=length(as.numeric(x$conf))/2){	# This is arbitary; basically just saying if we have to remove 50% of sample the cluster won't be used.  Need to evaluate this.
			warning(gettextf("Not enough data to generate canonical clusters"), domain = NA)
			model <- 0
			return(model)
		}else{
			modelData <- data.frame( R=as.numeric(as.character(newX$R)), theta=as.numeric(as.character(newX$theta)) )
			model <- lm(R ~ theta, data=modelData)
			return(model)
		}
	}
	# Apply
	baf <- lapply(data,clusterSummary)
	lrr <- lapply(data,clusterModel)
	CC <- list(model=lrr,sum=baf)
	return(CC)
}