# Function for retriving LL from PennCNV and calculating LL weights
# Usage: mapply(pcnv.acghRegion, CNVID=<>, ID=<>, TRIOPATH=<>,DATAPATH=<>,SIMPLIFY=F)

pcnv.acghRegion <- function(CNVID, ID, TRIOPATH, DATAPATH){
	cnvid <- as.character(CNVID)
	print(cnvid)
	FAM <- as.character(ID)
	CHR <- as.numeric(unlist(strsplit(as.character(cnvid),":"))[1])
	START <- as.numeric(unlist(strsplit(as.character(cnvid),":"))[2])
	STOP <- as.numeric(unlist(strsplit(as.character(cnvid),":"))[3])
	map <- read.table(paste("/lustre/scratch107/projects/ddd/users/km5/theKnaves/maps/chr",CHR,".map",sep=""),h=F,colClasses=c("numeric","character","numeric"))
	map <- map[order(map[,3]),]
	names(map) <- c("chr","SNP","position")
	currentRegion <- map[which(map$position >= START & map$position <= STOP),]
	STARTSNP <- ifelse(length(currentRegion[,1]) > 0, currentRegion$SNP[which(currentRegion$position==min(currentRegion$position))], NA)	
	ENDSNP <- ifelse(length(currentRegion[,1]) > 0, currentRegion$SNP[which(currentRegion$position==max(currentRegion$position))], NA)
	nullReport <- data.frame(cbind(as.character(FAM),as.character(cnvid),as.character(STARTSNP),as.character(ENDSNP),as.character(min(currentRegion$position)), as.character(max(currentRegion$position)),as.character(-1),as.character(max(currentRegion$position)-min(currentRegion$position)),as.character(NA),as.character(NA),as.character(NA),as.character(NA),as.character(NA),paste(unlist(strsplit(date()," "))[4],unlist(strsplit(date()," "))[2],unlist(strsplit(date()," "))[3],"-",unlist(strsplit(date()," "))[5],sep=""),packageDescription("theKnaves",field="Version")))
	names(nullReport) <- c("sanger_id","cnv_id","start_snp","stop_snp","bp_start","bp_stop","inc_snps","region_length","cns_father","cns_mother","cns_offspring","wlik","pDN","date_processed","pipeline_version")
	filecount = 0
	if (length(currentRegion[,1]) >= 3){
		if (file.exists(paste(TRIOPATH,"/",FAM,"/trio.txt",sep=""))){
			setwd(paste(TRIOPATH,"/",FAM,sep=""))
			trioList <- read.table("trio.txt",colClasses="character",h=F)
			names(trioList) <- c("offspring","father","mother")
			filecount=0
			for (ID in trioList){
				FILENAME <- paste(DATAPATH,"/",ID,"/adj-",CHR,".txt",sep="")
				if (file.exists(FILENAME) ){
					a <- try(read.table(FILENAME,h=T,colClasses=c("character","numeric","numeric")))
					names(a)[1] = "SNP"
					pcnv.writeData(a, CHR, ID)
					filecount = filecount+1
				}else{
					next
				}
			}
			if (filecount!=3){	# Cannot find all parents, call offspring only
				FILENAME <- paste(DATAPATH,"/",trioList[1],"/adj-",CHR,".txt",sep="")
				if (file.exists(FILENAME) ){
					a <- try(read.table(FILENAME,h=T,colClasses=c("character","numeric","numeric")))
					names(a)[1] = "SNP"
					pcnv.writeData(a, CHR, ID)
					# Write individual list
					list <- paste(getwd(),"/",trioList[1],"-",CHR,".txt",sep="")
					write.table(list,file=paste(getwd(),"/","indList-",CHR,".txt",sep=""),quote=F,row.names=F,col.names=F,sep = "\t")
					# Run PennCNV to get LL
					confRegions = data.frame(region=paste("chr",CHR,":",min(currentRegion$position),"-",max(currentRegion$position),sep=""),STARTSNP,ENDSNP,delfreq=0.001,dupfreq=0.001)
					write.table(confRegions, file="confregions", row.names=F, quote=F, col.names=F, sep="\t")
					REGIONLIST <- "-candlist confregions"		
					LIST <- paste(" --list indList-",CHR,".txt ",sep="") 
					LOG <- "acghRegion-loglik.log"
					OUTFILE <- "acghRegion.loglik"
					cmd <- paste("perl /nfs/ddd0/penncnv/detect_cnv.pl -validate -hmm /nfs/ddd0/penncnv/lib/hh550.hmm -pfb /lustre/scratch107/projects/ddd/users/km5/theKnaves/scotland-PFB.txt",LIST,REGIONLIST,"-valilog",LOG,"-out",OUTFILE,sep=" ")
					system(cmd)
					if (file.exists(LOG)){
						cnLL <- try(read.table(LOG,h=F,colClasses=c("character","character","numeric","numeric","numeric","numeric","numeric")))
						file.remove(LOG)
						llData <- cnLL[,3:7]
						if (sum(is.na(llData[,1])) == 0 && sum(llData) != 0){	# Avoids situation where too much missing data for PennCNV to make estimate
							trioCNstates <- expand.grid(kid = seq(1, 5, 1))
							calcLL <- function(trioCNstates){
								#cbind(trioCNstates[1],trioCNstates[2],trioCNstates[3],llData[1,trioCNstates[1]],llData[2,trioCNstates[2]],llData[3,trioCNstates[3]],(llData[1,trioCNstates[1]] + llData[2,trioCNstates[2]] + llData[3,trioCNstates[3]]))
								-2*(llData[1,trioCNstates[1]])
							}
							trioLL <- apply(trioCNstates,1,calcLL)
							tester <- cbind(trioCNstates,trioLL)
							tester <- cbind(tester,delta=tester$trioLL-min(tester$trioLL))
							weightTotal <- sum(exp(tester$delta*-1/2));
							tester <- cbind(tester,weight=((exp(tester$delta*-1/2))/weightTotal))
							state = subset(tester, weight==max(weight))
							report <- data.frame(cbind(as.character(FAM),as.character(cnvid),as.character(STARTSNP),as.character(ENDSNP),as.character(min(currentRegion$position)), as.character(max(currentRegion$position)), as.character(length(currentRegion[,1])),as.character(max(currentRegion$position)-min(currentRegion$position)),as.character(NA),as.character(NA),as.character(state[1]-1),as.character(state[4]),as.character(NA),paste(unlist(strsplit(date()," "))[4],unlist(strsplit(date()," "))[2],unlist(strsplit(date()," "))[3],"-",unlist(strsplit(date()," "))[5],sep=""),packageDescription("theKnaves",field="Version")))
							names(report) <- c("sanger_id","cnv_id","start_snp","stop_snp","bp_start","bp_stop","inc_snps","region_length","cns_father","cns_mother","cns_offspring","wlik","pDN","date_processed","pipeline_version")
							return(report)
						}else{
							return(nullReport)
						}
					}else{
						return(nullReport)
					}
				}else{
					return(nullReport)
				}			
			}else{
				# Write out trio file list in PennCNV format
				pcnv.writeLists(trioList, CHR)
				# Run PennCNV to get LL
				confRegions = data.frame(region=paste("chr",CHR,":",min(currentRegion$position),"-",max(currentRegion$position),sep=""),STARTSNP,ENDSNP,delfreq=0.001,dupfreq=0.001)
				write.table(confRegions, file="confregions", row.names=F, quote=F, col.names=F, sep="\t")
				REGIONLIST <- "-candlist confregions"		
				LIST <- paste(" --list indList-",CHR,".txt ",sep="") 
				LOG <- "acghRegion-loglik.log"
				OUTFILE <- "acghRegion.loglik"
				cmd <- paste("perl /nfs/ddd0/penncnv/detect_cnv.pl -validate -hmm /nfs/ddd0/penncnv/lib/hh550.hmm -pfb /lustre/scratch107/projects/ddd/users/km5/theKnaves/scotland-PFB.txt",LIST,REGIONLIST,"-valilog",LOG,"-out",OUTFILE,sep=" ")
				system(cmd)
				if (file.exists(LOG)){
					cnLL <- try(read.table(LOG,h=F,colClasses=c("character","character","numeric","numeric","numeric","numeric","numeric")))
					file.remove(LOG)
					llData <- cnLL[,3:7]
					if (sum(is.na(llData[,1])) == 0 && sum(llData) != 0){	# Avoids situation where too much missing data for PennCNV to make estimate
						trioCNstates <- expand.grid(dad = seq(1, 5, 1), mum = seq(1, 5, 1), kid = seq(1, 5, 1))
						calcLL <- function(trioCNstates){
							#cbind(trioCNstates[1],trioCNstates[2],trioCNstates[3],llData[1,trioCNstates[1]],llData[2,trioCNstates[2]],llData[3,trioCNstates[3]],(llData[1,trioCNstates[1]] + llData[2,trioCNstates[2]] + llData[3,trioCNstates[3]]))
							-2*(llData[1,trioCNstates[1]] + llData[2,trioCNstates[2]] + llData[3,trioCNstates[3]])
						}
						trioLL <- apply(trioCNstates,1,calcLL)
						tester <- cbind(trioCNstates,trioLL)
						tester <- cbind(tester,delta=tester$trioLL-min(tester$trioLL))
						weightTotal <- sum(exp(tester$delta*-1/2));
						tester <- cbind(tester,weight=((exp(tester$delta*-1/2))/weightTotal))
						state = subset(tester, weight==max(weight))
						pDN <- sum(tester$weight[!apply(tester,1,function(x) x[3] %in% c(x[1],x[2]))])
						report <- data.frame(cbind(as.character(FAM),as.character(cnvid),as.character(STARTSNP),as.character(ENDSNP),as.character(min(currentRegion$position)), as.character(max(currentRegion$position)), as.character(length(currentRegion[,1])),as.character(max(currentRegion$position)-min(currentRegion$position)),as.character(state[1]-1),as.character(state[2]-1),as.character(state[3]-1),as.character(state[6]),as.character(pDN),paste(unlist(strsplit(date()," "))[4],unlist(strsplit(date()," "))[2],unlist(strsplit(date()," "))[3],"-",unlist(strsplit(date()," "))[5],sep=""),packageDescription("theKnaves",field="Version")))
						names(report) <- c("sanger_id","cnv_id","start_snp","stop_snp","bp_start","bp_stop","inc_snps","region_length","cns_father","cns_mother","cns_offspring","wlik","pDN","date_processed","pipeline_version")
						return(report)
					}else{
						return(nullReport)
					}
				}else{
					return(nullReport)
				}
			}
		}else{
			return(nullReport)
		}
	}else{
		if (length(currentRegion[,1]) > 1){
			report <- data.frame(cbind(as.character(FAM),as.character(cnvid),as.character(STARTSNP),as.character(ENDSNP),as.character(min(currentRegion$position)), as.character(max(currentRegion$position)),as.character(length(currentRegion[,1])),as.character(max(currentRegion$position)-min(currentRegion$position)),as.character(NA),as.character(NA),as.character(NA),as.character(NA),as.character(NA),paste(unlist(strsplit(date()," "))[4],unlist(strsplit(date()," "))[2],unlist(strsplit(date()," "))[3],"-",unlist(strsplit(date()," "))[5],sep=""),packageDescription("theKnaves",field="Version")))
			names(report) <- c("sanger_id","cnv_id","start_snp","stop_snp","bp_start","bp_stop","inc_snps","region_length","cns_father","cns_mother","cns_offspring","wlik","pDN","date_processed","pipeline_version")
			return(report)
		}else{
			report <- data.frame(cbind(as.character(FAM),as.character(cnvid),as.character(STARTSNP),as.character(ENDSNP),as.character(NA), as.character(NA),as.character(length(currentRegion[,1])),as.character(NA),as.character(NA),as.character(NA),as.character(NA),as.character(NA),as.character(NA),paste(unlist(strsplit(date()," "))[4],unlist(strsplit(date()," "))[2],unlist(strsplit(date()," "))[3],"-",unlist(strsplit(date()," "))[5],sep=""),packageDescription("theKnaves",field="Version")))
			names(report) <- c("sanger_id","cnv_id","start_snp","stop_snp","bp_start","bp_stop","inc_snps","region_length","cns_father","cns_mother","cns_offspring","wlik","pDN","date_processed","pipeline_version")
			return(report)
		}
	}
}
