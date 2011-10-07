jointPennCNV <- function(chrStart=1,chrStop=22,outDir=getwd()){
	print("Warning: This function assumes that you have run the processTrio command first.")
	chromosomes <- seq(chrStart,chrStop,1)
	callCNV <- function(CHR){
		PFB <- paste(" -pfb ", "/lustre/scratch103/sanger/km5/ddd/snp-cnv/scotland/penncnv_in/pfb/scotland-PFB-",CHR,".txt ",sep="") 
		LIST <- paste(" --list ", outDir,"/","famList.txt ",sep="") 
		LOG <- paste(" -log ",outDir,"/",CHR,"-joint.log ",sep="") 
		OUTFILE <- paste(outDir,"/",CHR,".jointcnv",sep="")
		cmd <- paste("perl /nfs/team143/ddd/bin/penncnv/detect_cnv.pl -joint -hmm /nfs/team143/ddd/bin/penncnv/lib/hh550.hmm ",PFB,LIST,LOG," -out ",OUTFILE,sep="")
		system(cmd)
		if (file.exists(OUTFILE) ){
			chrCalls <- read.table(OUTFILE,h=F,colClasses="character")
		}else{
			chrCalls <- 0
		}
		return(chrCalls)
	}
	jointCalls <- lapply(chromosomes,callCNV)
	listLabel <- c()
	for (i in chromosomes){
		listLabel <- c(listLabel,paste("chr",i,sep=""))
	}
	names(jointCalls) <- listLabel
	return(jointCalls)
}