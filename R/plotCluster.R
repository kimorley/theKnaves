plotCluster <- function(CHR, SNP, targetINTU, targetGTU, clusterPATH, clusterDATA, outDIR=getwd(), save=F, print=T){
	# Checks
	if (missing(CHR)){
		stop("Supply chromosome on which SNP is located.")
	}
	if (missing(SNP)){
		stop("Supply SNP for which clusters are to be plotted.")
	}
	# Load model for SNP
	# Load canonical cluster for SNP
	if (file.exists(paste(clusterPATH,"/chr",CHR,"/clusterFile-",SNP,".RData",sep=""))){
		loadCluster <- try(load(paste(clusterPATH,"/chr",CHR,"/clusterFile-",SNP,".RData",sep="")))
		if (inherits(loadCluster, 'try-error')){
			print(paste("Canonical cluster file for SNP ",SNP," does not exist in specified location.",sep=""))
			return(clusterData <- list(model=NULL,summary=NULL))
		}
	}else{
		print(paste("Canonical cluster file for SNP ",SNP," does not exist in specified location.",sep=""))
		return(clusterData <- list(model=NULL,summary=NULL))
	}
	# Read target data
	target <- convertRawSnp(SNP,targetINTU, targetGTU)
	data <- target[[1]]
	# Read reference data
	ref <- convertRawSnp(SNP,paste(clusterDATA,"/intu/scot-founders-",CHR,".intu",sep=""),paste(clusterDATA,"/gtu/scot-founders-",CHR,".gtu",sep=""))
	reftemp <- ref[[1]]
	refdata <- subset(reftemp, call != "NN" & as.numeric(as.character(conf)) >= 0.99)
	# Plotting data	
	testData <- data.frame(t=seq(0,1,0.01))
	testData <- cbind(testData,r=predict(clusterData$model,testData))

	clusterPlot <- ggplot() + 
			layer(
					data = refdata, mapping = aes(x = t, y = r, colour = call), 
					alpha = 0.1,
					geom = "point"				
			) +
			layer(
					data = data, mapping = aes(x = t, y = r, colour = call),
					geom = "point"
			) + 
			layer(
					data = testData, mapping = aes(x = t, y = r),
					geom = "line"
			) +
			layer(
					data= het, mapping = aes(x = x, y = y),
					size = 0.5,
					linetype = 2,
					geom = "path"
			) +
			layer(
					data= hom1, mapping = aes(x = x, y = y),
					size = 0.5,
					linetype = 2,
					geom = "path"
			) +
			layer(
					data= hom2, mapping = aes(x = x, y = y),
					size = 0.5,
					linetype = 2,
					geom = "path"
			) +
			scale_colour_brewer(palette="Set1") + ylab("R") + xlab(expression(theta)) + opts(title = SNP)
	if (print){
		print(clusterPlot)
	}
	if (save){
		ggsave(file=paste(outDIR,"/clusterPlot-",SNP,".pdf",sep=""),height=8,width=8)
	}
}