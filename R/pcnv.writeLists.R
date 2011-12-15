# For set of trio ID, writes out individual file list and trio file lists for PennCNV
pcnv.writeLists <- function(trioIDs, CHR, outDIR=getwd()){
	print("Warning: this function assumes the trio IDs are in PLINK pedigree format (child, father, mother).")
	# Write trio list
	list <- c(paste(outDIR,"/",trioIDs[2],"-",CHR,".txt",sep=""),paste(outDIR,"/",trioIDs[3],"-",CHR,".txt",sep=""),paste(outDIR,"/",trioIDs[1],"-",CHR,".txt",sep=""))
	write.table(t(list),file=paste(outDIR,"/","famList-",CHR,".txt",sep=""),quote=F,row.names=F,col.names=F,sep = "\t")
	# Write individual list
	list <- c(paste(outDIR,"/",trioIDs[2],"-",CHR,".txt",sep=""),paste(outDIR,"/",trioIDs[3],"-",CHR,".txt",sep=""),paste(outDIR,"/",trioIDs[1],"-",CHR,".txt",sep=""))
	write.table(list,file=paste(outDIR,"/","indList-",CHR,".txt",sep=""),quote=F,row.names=F,col.names=F,sep = "\t")
}
