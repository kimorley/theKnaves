# Writes out the main input files for PennCNV with appropriate headers
writePennCNVdata <- function(data, CHR, ID, outDIR=getwd()){
	outData <- subset(data,select=c(SNP,BAF,waveLRR))
	names(outData) <- c("Name",paste(ID,"-",CHR,"A.B Allele Freq",sep=""),paste(ID,"-",CHR,"A.Log R Ratio",sep=""))
	write.table(outData,file=paste(outDIR,"/",ID,"-",CHR,".txt",sep=""),quote=F,row.names=F,sep = "\t")
}

writePennCNVlists <- function(trioIDs, CHR, outDIR=getwd()){
	# Write trio list
	list <- c(paste(outDIR,"/",trioIDs[3],"-",CHR,".txt",sep=""),paste(outDIR,"/",trioIDs[2],"-",CHR,".txt",sep=""),paste(outDIR,"/",trioIDs[1],"-",CHR,".txt",sep=""))
	write.table(t(list),file=paste(outDIR,"/","famList-",CHR,".txt",sep=""),quote=F,row.names=F,col.names=F,sep = "\t")
	# Write individual list
	list <- c(paste(outDIR,"/",trioIDs[3],"-",CHR,".txt",sep=""),paste(outDIR,"/",trioIDs[2],"-",CHR,".txt",sep=""),paste(outDIR,"/",trioIDs[1],"-",CHR,".txt",sep=""))
	write.table(list,file=paste(outDIR,"/","indList-",CHR,".txt",sep=""),quote=F,row.names=F,col.names=F,sep = "\t")
}

writePennCNVTrio <- function(trioData,outDir=getwd()){
	# Check
	if (!is.list(trioData)){
		stop("Input should be list of trio data lists generated using processTrio command.")
	}
	# Setup
	trioIDs <- c()
	# Execute
	# Write PennCNV input files
	for (i in 1:3){
		ID <- strsplit(as.character(unique(trioData[,i]$data$SampleID)),"_")[[1]][length(strsplit(as.character(unique(trioData[,i]$data$SampleID)),"_")[[1]])]
		trioIDs <- c(trioIDs,ID)
		outData <- subset(trioData[,i]$data,select=c(SNP,BAF,waveLRR))
		names(outData) <- c("Name",paste(ID,"A.B Allele Freq",sep=""),paste(ID,"A.Log R Ratio",sep=""))
		write.table(outData,file=paste(outDir,"/",ID,".txt",sep=""),quote=F,row.names=F,sep = "\t")
	}
	# Write trio list
	list <- c(paste(outDir,"/",trioIDs[1],".txt",sep=""),paste(outDir,"/",trioIDs[2],".txt",sep=""),paste(outDir,"/",trioIDs[3],".txt",sep=""))
	write.table(t(list),file=paste(outDir,"/","famList.txt",sep=""),quote=F,row.names=F,col.names=F,sep = "\t")
	# Write individual list
	list <- c(paste(outDir,"/",trioIDs[1],".txt",sep=""),paste(outDir,"/",trioIDs[2],".txt",sep=""),paste(outDir,"/",trioIDs[3],".txt",sep=""))
	write.table(list,file=paste(outDir,"/","indList.txt",sep=""),quote=F,row.names=F,col.names=F,sep = "\t")
}
