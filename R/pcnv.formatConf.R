# Processes output file from PennCNV -conf command
pcnv.formatConf <- function(cnCall){
	cnProbe <- as.numeric(strsplit(as.character(cnCall[2]),"=")[[1]][2])
	cnLength <- strsplit(strsplit(as.character(cnCall[3]),"=")[[1]][2],",")[[1]][1]
	if (length(strsplit(strsplit(as.character(cnCall[3]),"=")[[1]][2],",")[[1]]) > 1){
		for (i in 2:length(strsplit(strsplit(as.character(cnCall[3]),"=")[[1]][2],",")[[1]])){
			cnLength <- paste(cnLength,strsplit(strsplit(as.character(cnCall[3]),"=")[[1]][2],",")[[1]][i],sep="")
		}
	}
	cnLength <- as.numeric(cnLength)
	chr <- strsplit(strsplit(as.character(cnCall[1]),":")[[1]][1],"r")[[1]][2]
	start <- strsplit(strsplit(as.character(cnCall[1]),":")[[1]][2],"-")[[1]][1]
	stop <- strsplit(strsplit(as.character(cnCall[1]),":")[[1]][2],"-")[[1]][2]
	cn <- as.numeric(strsplit(as.character(cnCall[4]),"=")[[1]][2])
	startSnp <- strsplit(as.character(cnCall[6]),"=")[[1]][2] 
	stopSnp <- strsplit(as.character(cnCall[7]),"=")[[1]][2]
	conf <- strsplit(as.character(cnCall[8]),"=")[[1]][2]
	temp <- cbind(chr,start,stop,cn,startSnp,stopSnp,cnProbe,cnLength,conf)
	return(temp)
}
