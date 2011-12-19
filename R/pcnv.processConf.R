pcnv.processConf <- function(x){
	cnProbe <- as.numeric(strsplit(as.character(x[2]),"=")[[1]][2])
	cnLength <- strsplit(strsplit(as.character(x[3]),"=")[[1]][2],",")[[1]][1]
	if (length(strsplit(strsplit(as.character(x[3]),"=")[[1]][2],",")[[1]]) > 1){
		for (i in 2:length(strsplit(strsplit(as.character(x[3]),"=")[[1]][2],",")[[1]])){
			cnLength <- paste(cnLength,strsplit(strsplit(as.character(x[3]),"=")[[1]][2],",")[[1]][i],sep="")
		}
	}
	cnLength <- as.numeric(cnLength)
	chr <- strsplit(strsplit(as.character(x[1]),":")[[1]][1],"r")[[1]][2]
	start <- strsplit(strsplit(as.character(x[1]),":")[[1]][2],"-")[[1]][1]
	stop <- strsplit(strsplit(as.character(x[1]),":")[[1]][2],"-")[[1]][2]
	cn <- as.numeric(strsplit(as.character(x[4]),"=")[[1]][2])
	who <- unlist(strsplit(as.character(unlist(strsplit(as.character(x[5]),"\\."))[length(unlist(strsplit(as.character(x[5]),"/")))]),"-"))[1]
	trioid <- names(trioList)[which(who==trioList)]
	startSnp <- strsplit(as.character(x[6]),"=")[[1]][2] 
	stopSnp <- strsplit(as.character(x[7]),"=")[[1]][2]
	conf <- as.numeric(strsplit(as.character(x[8]),"=")[[1]][2])
	cnvid <- as.character(paste(chr,start,stop,sep=":"))
	temp <- cbind(chr,start,stop,cn,startSnp,stopSnp,cnProbe,cnLength,conf,trioid,cnvid)
	return(temp)
}
