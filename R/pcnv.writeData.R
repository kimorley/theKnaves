# Writes out the main input files for PennCNV with appropriate headers
pcnv.writeData <- function(data, CHR, ID, outDIR=getwd()){
	ID <- as.character(ID)
	CHR <- as.numeric(CHR)
	if (sum(c("SNP","BAF","waveLRR")%in%names(data))!=3){
		stop("Did not find necessary variables: SNP, BAF, waveLRR")
	}
	outData <- subset(data,select=c(SNP,BAF,waveLRR))
	names(outData) <- c("Name",paste(ID,"-",CHR,"A.B Allele Freq",sep=""),paste(ID,"-",CHR,"A.Log R Ratio",sep=""))
	write.table(outData,file=paste(outDIR,"/",ID,"-",CHR,".txt",sep=""),quote=F,row.names=F,sep = "\t")
}


