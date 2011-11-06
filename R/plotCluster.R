plotCluster <- function(target,SNP,clusterDIR=getwd(),save=F,print=T){
	if (missing(target)){
		stop("Please supply target list object.")
	}
	if (missing(SNP)){
		stop("Supply SNP for which clusters are to be plotted.")
	}
	if (file.exists(paste(clusterDIR,"/clusterFile-",SNP,".gzip",sep=""))){
		loadCluster <- try(load(paste(clusterDIR,"/clusterFile-",SNP,".gzip",sep="")))
		if (inherits(loadCluster, 'try-error')){
			stop("Canonical cluster file for SNP ",SNP," does not exist in specified location.")
		}
	}else{
		stop("Canonical cluster file for SNP ",SNP," does not exist in specified location.")
	}
	data <- target[names(target)==SNP][[1]]	
	testData <- data.frame(t=seq(0,1,0.01))
	testData <- cbind(testData,r=predict(clusterData$model,testData))
	clusterPlot <- ggplot() + 
			layer(
					data = data, mapping = aes(x = t, y = r, colour = call),
					geom = "point"
			) + 
			layer(
					data = testData, mapping = aes(x = t, y = r),
					geom = "line"
			) +
			scale_colour_brewer(palette="Set1") + ylab("R") + xlab(expression(theta)) + opts(title = SNP)
	if (print){
		print(clusterPlot)
	}
	if (save){
		ggsave(file="clusterPlot.pdf",height=5,width=5)
	}
}