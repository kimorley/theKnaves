processFinalReport <- function(ID,inDir=getwd(),outDir=getwd(),chrStart=1,chrStop=22,writeOut=FALSE){
	# Check
	if (length(ID)!=1){
		stop("ID should be single identifier.")
	}
	# Setup
	header <- c("CHR","Median","Mean","SD","MAD")
	qcBAF <- data.frame()
	qcLRR <- data.frame()
	qcLRRadj <- data.frame()
	qcMarker <- data.frame()
	qcData <- data.frame()
	# Execution
	for (CHR in chrStart:chrStop){
		data <- read.table(paste(inDir,"/",ID,"/","FR-",CHR,".txt",sep=""),h=T)
		data <- data[order(data[,4]), ]
		fr <- data[which(data$Position > 0),] # Remove unmapped markers
		qcBAF <- rbind(qcBAF,c(CHR,summaryBAF(fr$BAF)))
		qcLRR <- rbind(qcLRR,c(CHR,summaryLRR(fr$LRR)))
		fr <- cbind(fr,"LRRadj"=adjLRR(fr$LRR))
		qcLRRadj <- rbind(qcLRRadj,c(CHR,summaryLRR(fr$LRRadj)))
		qcMarker <- rbind(qcMarker,c(CHR,summaryMarkers(fr$Position)))
		qcData <- rbind(qcData,fr)
	}
	qcBAF <- rbind(qcBAF,c("Genome",summaryBAF(qcData$BAF)))
	qcLRR <- rbind(qcLRR,c("Genome",summaryLRR(qcData$LRR)))
	qcLRRadj <- rbind(qcLRRadj,c("Genome",summaryLRR(qcData$LRRadj)))
	qcMarker <- rbind(qcMarker,c("Genome",summaryMarkers(qcData$Position)))
	names(qcBAF) <- header
	names(qcLRR) <- header
	names(qcLRRadj) <- header
	names(qcMarker) <- c("CHR","Count","Median_InterMarkerDist")
	output <- list(BAF=qcBAF,LRR=qcLRR,LRRadj=qcLRRadj,markerInfo=qcMarker,data=qcData)
	if (writeOut){
		write.table(qcBAF,file=paste(outDir,"/",ID,"-BAF.txt",sep=""),quote=F,row.names=F)
		write.table(qcLRR,file=paste(outDir,"/",ID,"-LRR.txt",sep=""),quote=F,row.names=F)
		write.table(qcLRRadj,file=paste(outDir,"/",ID,"-LRRadj.txt",sep=""),quote=F,row.names=F)
		write.table(qcMarker,file=paste(outDir,"/",ID,"-Marker.txt",sep=""),quote=F,row.names=F)
	}
	return(output)
}