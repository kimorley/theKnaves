# Writes out the main input files for PennCNV with appropriate headers
pcnv.writeData <- function(data, CHR, ID, outDIR){
	ID <- as.character(ID)
	CHR <- as.numeric(CHR)
	names(data) <- c("Name",paste(ID,"-",CHR,"A.B Allele Freq",sep=""),paste(ID,"-",CHR,"A.Log R Ratio",sep=""))
	write.table(data,file=paste(outDIR,"/",ID,"-",CHR,".txt",sep=""),quote=F,row.names=F,sep = "\t")
}


