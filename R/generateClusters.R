# Generate canonical clusters given a data file containing R and theta (run convertRawFiles or convertRawSnp first)
generateClusters <- function(data, DIR=getwd(), confThresh = 0.99, clusterStat = "median", saveCluster=FALSE){
	# Check
	if (missing(data)){
		stop("Must supply 'data' argument")
	}
	if (!is.list(data)){
		stop("Input data must be in list format")
	}
	if (clusterStat!="median" && clusterStat != "mean"){
		stop("Invalid summary statistic for clusters.  Select either 'median' or 'mean'.")
	}
	if (confThresh < 0 | confThresh > 1){
		stop("Invalid value for confidence threshold.  Please supply value between 0 and 1.")
	}
	# Functions
	cluster <- function(x){
		mySnp <- names(eval(as.list(sys.call(-1))[[2]]))[as.numeric(gsub("[^0-9]", "", deparse(substitute(x))))]
		newX <- x[which( as.numeric(as.character(x$conf)) >= confThresh & as.character(x$call) != "NN"),]
		if (length(newX[,1])<=length(as.numeric(x$conf))/2){	# This is arbitary; basically just saying if we have to remove 50% of sample the cluster won't be used.  Need to evaluate this.
			warning(gettextf("Not enough data to generate canonical clusters"), domain = NA)
			summary <- 0
			model <- 0
			clusterData <- list(model=model,summary=summary)
			if (saveCluster){
				save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""))
			}
			return(clusterData)
		}else{
			summary <- tapply(newX$t,as.character(newX$call),clusterStat)
			# Check whether the results are sensible
			if (length(summary)==3){
				if ( summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))] > max(summary[which(substr(names(summary),1,1)==substr(names(summary),2,2))]) || summary[which(substr(names(summary),1,1)!=substr(names(summary),2,2))] < min(summary[which(substr(names(summary),1,1)==substr(names(summary),2,2))]) ){
					# Theta value for heterozygote is smaller than AA homozygote or larger than BB homozygote so we don't want to call from this
					summary <- 0
					model <- 0
					clusterData <- list(model=model,summary=summary)
					if (saveCluster){
						save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""))
					}
					return(clusterData)
				}else{
					modelTemp <- lm(r ~ ns(t, knots=summary[2], Boundary.knots=c(summary[1],summary[3])), data=newX, model=FALSE)
					model <- list()
					class(model) <- class(modelTemp)
					model$terms <- modelTemp$terms
					model$coefficients <- modelTemp$coefficients
					model$rank <- modelTemp$rank
					model$qr <- modelTemp$qr
					save(model,file=paste(mySnp,"-test.RData",sep=""),compression_level=9)
					summary <- rbind(tapply(newX$t,as.character(newX$call),clusterStat),table(as.character(newX$call)))
					clusterData <- list(model=model,summary=summary)
					if (saveCluster){
						save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""))
					}
					return(clusterData)
				}
			}else if (length(summary)==2){
				modelTemp <- lm(r ~ ns(t, Boundary.knots=c(summary[1],summary[2])), data=newX, model=FALSE)
				model <- list()
				class(model) <- class(modelTemp)
				model$coefficients <- modelTemp$coefficients
				model$terms <- modelTemp$terms
				model$rank <- modelTemp$rank
				model$qr <- modelTemp$qr
				summary <- rbind(tapply(newX$t,as.character(newX$call),clusterStat),table(as.character(newX$call)))
				clusterData <- list(model=model,summary=summary)
				if (saveCluster){
					save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""))
				}
				return(clusterData)
			}else if (length(summary)==1){
				modelTemp <- lm(r ~ ns(t, Boundary.knots=c(min(t),max(t))), data=newX, model=FALSE)
				model <- list()
				class(model) <- class(modelTemp)
				model$coefficients <- modelTemp$coefficients
				model$terms <- modelTemp$terms
				model$rank <- modelTemp$rank
				model$qr <- modelTemp$qr
				summary <- rbind(tapply(newX$t,as.character(newX$call),clusterStat),table(as.character(newX$call)))
				clusterData <- list(model=model,summary=summary)
				if (saveCluster){
					save(clusterData, file=paste(DIR,"/clusterFile-",mySnp,".RData",sep=""))
				}
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