pcnv.trioFiles <- function(FAM, TRIOPATH, DATAPATH){
	TRIO <- read.table(paste(TRIOPATH,"/",FAM,"/trio.txt",sep=""),colClasses="character",h=F)
	names(TRIO) <- c("offspring","father","mother")
	filecount = 0
	for (ID in TRIO){
		chrcount = 0
		for (CHR in seq(21,22,1)){
			FILENAME <- paste(DATAPATH,"/",ID,"/adj-",CHR,".txt",sep="")
			if (file.exists(FILENAME)){
				data = read.table(FILENAME,h=T,colClasses=c("character","numeric","numeric"))
				names(data) <- c("Name",paste(ID,"-",CHR,"A.B Allele Freq",sep=""),paste(ID,"-",CHR,"A.Log R Ratio",sep=""))
				write.table(data,file=paste(TRIOPATH,"/",FAM,"/",ID,"-",CHR,".txt",sep=""), row.names = F, quote = F, sep = "\t")
				chrcount = chrcount+1
			}else{
				next
			}
			if (chrcount == 2){
				filecount = filecount+1
			}
		}
	}
	if (filecount == 3){
		# Write trio list
		list <- c(paste(TRIO[2],".txt",sep=""),paste(TRIO[3],".txt",sep=""),paste(TRIO[1],".txt",sep=""))
		write.table(t(list),file=paste(TRIOPATH,"/",FAM,"/famList-",CHR,".txt",sep=""),quote=F,row.names=F,col.names=F,sep = "\t")
		# Write individual list
		list <- c(paste(TRIO[2],".txt",sep=""),paste(TRIO[3],".txt",sep=""),paste(TRIO[1],".txt",sep=""))
		write.table(list,file=paste(TRIOPATH,"/",FAM,"/indList-",CHR,".txt",sep=""),quote=F,row.names=F,col.names=F,sep = "\t")
	}
	return(filecount)
}