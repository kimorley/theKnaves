# Take list of regions (from standard PennCNV output) and estimate confidence
pcnv.callConf <- function(callList, datafile, delFreq=0.01, dupFreq=0.01){
	# Function for extracting regions
	getConfidence <- function(callList){
		startSnp <- unlist(strsplit(as.character(callList[6]),"="))[2]
		endSnp <- unlist(strsplit(as.character(callList[7]),"="))[2]
		region <- cbind(callList[1],startSnp,endSnp,delFreq,dupFreq)
	}
	confRegions <- data.frame(t(apply(callList,1,getConfidence)))
	write.table(confRegions, file="confregions", row.names=F, quote=F, col.names=F, sep="\t")
	REGIONLIST <- "-candlist confregions"
	FILE <- datafile
	cmd <- paste("perl /nfs/team143/ddd/bin/penncnv/detect_cnv.pl -validate -hmm /nfs/team143/ddd/bin/penncnv/lib/hh550.hmm -pfb /lustre/scratch107/projects/ddd/users/km5/theKnaves/scotland-PFB.txt ",FILE,REGIONLIST," -out outfile -conf",sep=" ")
	system(cmd)
	if (file.exists("outfile") ){
		cnConf <- try(read.table("outfile",h=F,colClasses="character"), silent=T)
		if (inherits(cnConf, 'try-error')){
			cnConf <- 0 
		} 
	}else{
		cnConf <- 0
	}
	return(cnConf)
}