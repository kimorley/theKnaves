# Usage: jointCalls <- lapply(chromosomes, pcnv.callJoint, DATAPATH = DATAPATH, TRIO = trioList)
pcnv.callJoint <- function(CHR, DATAPATH, TRIO){
	print(CHR)
	filecount = 0
	for (ID in TRIO){
		FILENAME <- paste(DATAPATH,"/",ID,"/adj-",CHR,".txt",sep="")
		if (file.exists(FILENAME)){
			data = read.table(FILENAME,h=T,colClasses=c("character","numeric","numeric"))
			names(data) <- c("Name",paste(ID,"-",CHR,"A.B Allele Freq",sep=""),paste(ID,"-",CHR,"A.Log R Ratio",sep=""))
			write.table(data,file=paste(ID,".txt",sep=""), row.names = F, quote = F, sep = "\t")
			filecount = filecount+1
		}else{
			next
		}
	}
	if (filecount==3){
		# Write trio list
		list <- c(paste(TRIO[2],".txt",sep=""),paste(TRIO[3],".txt",sep=""),paste(TRIO[1],".txt",sep=""))
		write.table(t(list),file=paste("famList-",CHR,".txt",sep=""),quote=F,row.names=F,col.names=F,sep = "\t")
		# Write individual list
		list <- c(paste(TRIO[2],".txt",sep=""),paste(TRIO[3],".txt",sep=""),paste(TRIO[1],".txt",sep=""))
		write.table(list,file=paste("indList-",CHR,".txt",sep=""),quote=F,row.names=F,col.names=F,sep = "\t")
		# Run joint calling in PennCNV
		LIST <- paste(" --list ", getwd(),"/","famList-",CHR,".txt ",sep="") 
		cmd <- paste("perl /nfs/ddd0/penncnv/detect_cnv.pl -joint -hmm /nfs/ddd0/penncnv/lib/hh550.hmm -pfb /lustre/scratch107/projects/ddd/users/km5/theKnaves/scotland-PFB.txt",LIST,"-out out1",sep=" ")
		system(cmd)
		if (file.exists("out1") ){
			cnCall <- try(read.table("out1",h=F,colClasses="character"), silent=T)
			if (inherits(cnCall, 'try-error')){
				cnCall <- 0 
			} 
		}else{
			cnCall <- 0
		}
		#file.remove("outfile")
		if (sum(length(cnCall)) != 9){
			#file.remove(paste(TRIO[1],".txt",sep=""),paste(TRIO[2],".txt",sep=""),paste(TRIO[3],".txt",sep=""))
			return(NULL)
		}else{
			report <- data.frame(t(apply(cnCall,1,pcnv.formatJointCall)))
			names(report) <- c("chr","start_bp","stop_bp","cn","start_snp","stop_snp","cnv_probes","cnv_length","id")
			findDeNovo <- function(x){
				all <- pcnv.checkInherit(x)
				if (sum(length(all))==0){
					return(NULL)
				}else{
					getConfidence <- function(x){
						label <- paste("chr",x[1],":",x[2],"-",x[3],sep="")
						region <- cbind(label,x[5],x[6],0.01,0.01)
					}
					confRegions <- data.frame(t(apply(all,1,getConfidence)))
					write.table(confRegions, file="confregions", row.names=F, quote=F, col.names=F, sep="\t")
					REGIONLIST <- "-candlist confregions"
					LIST <- paste(" --list ", getwd(),"/","indList-",CHR,".txt ",sep="") 
					cmd <- paste("perl /nfs/ddd0/penncnv/detect_cnv.pl -validate -hmm /nfs/ddd0/penncnv/lib/hh550.hmm -pfb /lustre/scratch107/projects/ddd/users/km5/theKnaves/scotland-PFB.txt ",LIST,REGIONLIST," -out out2 -conf",sep=" ")
					system(cmd)
					if (file.exists("out2") ){
						cnConf <- try(read.table("out2",h=F,colClasses="character"), silent=T)
						if (inherits(cnConf, 'try-error')){
							cnConf <- 0 
						} 
					}else{
						cnConf <- 0
					}
					pConf <- data.frame(t(apply(cnConf,1,pcnv.processConf)))
					o <- subset(pConf, X10=="offspring")
					p <- subset(pConf, X10=="father")
					m <- subset(pConf, X10=="mother")
					pcnv.combineInfo <- function(x){
						if (x[10] %in% o[,11]){
							x <- c(x,as.numeric(as.character(o[which(o$X11==as.character(x[10])),9])))
						}else{
							x <- c(x,0)
						}
						if (x[10] %in% p[,11]){
							x <- c(x,as.numeric(as.character(p[which(p$X11==as.character(x[10])),9])))
						}else{
							x <- c(x,0)
						}
						if (x[10] %in% m[,11]){
							x <- c(x,as.numeric(as.character(m[which(m$X11==as.character(x[10])),9])))
						}else{
							x <- c(x,0)
						}
						return(x)
					}
					summary <- data.frame(t(apply(all,1,pcnv.combineInfo)))
					summary <- cbind(summary,paste(unlist(strsplit(date()," "))[4],unlist(strsplit(date()," "))[2],unlist(strsplit(date()," "))[6],"-",unlist(strsplit(date()," "))[5],sep=""),packageDescription("theKnaves",field="Version"))
					names(summary)[12:16] <- c("off_conf","pat_conf","mat_conf","date_processed","pipe_version")
					file.remove("out1","out2","confregions",paste(TRIO[1],".txt",sep=""),paste(TRIO[2],".txt",sep=""),paste(TRIO[3],".txt",sep=""))
					return(summary)
				}
			}
			results <- findDeNovo(report)
		}
	}else{
		return(NULL)
	}
}



