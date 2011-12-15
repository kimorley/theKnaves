# Call chromosome-wide for single sample in PennCNV
# Prefix should be file name (assuming format prefixCHR.txt)
pcnv.callSingle <- function(prefix, START, STOP){
	chromosomes <- seq(START,STOP,1)
	callCNV <- function(CHR){
		FILENAME <- paste(prefix,CHR,".txt",sep="")
		if (file.exists(FILENAME) ){
#			data <- try(read.table(FILENAME,h=T,colClasses=c("character","numeric","numeric")))
			data <- try(read.table(FILENAME,h=T,))
			data <- subset(data, select=c(SNPName,BAF,LRR))
			if (inherits(data, 'try-error')){
				return(NULL)
			}else{
				names(data) <- c("SNP","BAF","waveLRR")
				pcnv.writeData(data, CHR=CHR, ID="adj")
				FILE <- paste("adj-",CHR,".txt",sep="")
				cmd <- paste("perl /nfs/team143/ddd/bin/penncnv/detect_cnv.pl -test -hmm /nfs/team143/ddd/bin/penncnv/lib/hh550.hmm -pfb /lustre/scratch107/projects/ddd/users/km5/theKnaves/scotland-PFB.txt ",FILE,"-log logfile -out outfile",sep=" ")
				system(cmd)
				if (file.exists("outfile") ){
					cnCall <- try(read.table("outfile",h=F,colClasses="character"), silent=T)
					if (inherits(cnCall, 'try-error')){
						cnCall <- 0 
					} 
				}else{
					cnCall <- 0
				}
				if (sum(length(cnCall)) != 7){
					file.remove("outfile","logfile")
					return(NULL)
				}else{
					cnConf <- pcnv.callConf(cnCall, FILE)
					report <- data.frame(t(apply(cnCall,1,pcnv.formatConf)))
					report <- cbind(report,paste(unlist(strsplit(date()," "))[4],unlist(strsplit(date()," "))[2],unlist(strsplit(date()," "))[6],"-",unlist(strsplit(date()," "))[5],sep=""),packageDescription("theKnaves",field="Version"))
					names(report) <- c("chr","start_bp","stop_bp","cn","start_snp","stop_snp","cnv_probes","cnv_length","cnv_conf","date_processed","pipe_version")
					return(report)
					file.remove("outfile","confregions","logfile")
				}
			}
		}else{
			return(NULL)
		}
	}
	singleCalls <- lapply(chromosomes,callCNV)
	listLabel <- c()
	for (i in chromosomes){
		listLabel <- c(listLabel,paste("chr",i,sep=""))
	}
	names(singleCalls) <- listLabel
	singleCalls <- singleCalls[!sapply(singleCalls, is.null)]
	singleCalls <- do.call(rbind, singleCalls)
	return(singleCalls)
}
