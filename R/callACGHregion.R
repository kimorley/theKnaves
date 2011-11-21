# Function for calling region identified by aCGH in SNP array data.



acgh <- read.table("aCGH-results.txt", as.is=TRUE, header=TRUE, sep="\t")
trioList <- scan("trio.txt",what=c("character","character","character"))
pfb <- read.table("scotland-PFB.txt",h=F,colClasses=c("character","numeric","numeric","numeric"))
pfb <- pfb[order(pfb[,2], pfb[,3]),]
names(pfb) <- c("snp","chr","position","pfb")

chromosomes <- sort(as.numeric(sapply(unique(acgh$chr),function(x){unlist(strsplit(x,"r"))[2]})))

for (i in sort(as.numeric(sapply(unique(acgh$chr),function(x){unlist(strsplit(x,"r"))[2]})))){
	# Read file for individual and QC
	for (i in trioList){
		adjData <- adjustLRR(dataFile=paste(i,"/raw-",CHR,".txt",sep=""),mapFile=paste("chr",CHR,".map",sep=""))
		writePennCNVdata(adjData, CHR, i)
	}
	# Write out trio file list in PennCNV format
	writePennCNVlists(trioList, CHR)
}

CHR <- 1
START <- 202259070
STOP <- 202272604
currentRegion <- subset(pfb, chr==CHR & position >= START & position <= STOP)

if (length(currentRegion[,1]) < 2){
	print(paste("Not enough probes to call region ",CHR,":","START","-","STOP",sep=""))
}else{
	# Call region
}

pennCnvLogLik <- function(trioIDs, cnvID, STARTSNP, ENDSNP, CHR, outDIR=getwd()){
	PFB <- paste(" -pfb ",outDIR,"/region-pfb.txt",sep="") 
	LIST <- paste(" --list ", outDIR,"/","indList-",CHR,".txt ",sep="") 
	LOG <- paste(" -valilog ", outDIR,"/",cnvID,"-loglik.log ",sep="") 
	REGION <- paste("-delfreq 0.001 -dupfreq 0.0001 -startsnp", STARTSNP, "-endsnp", ENDSNP, sep=" ")
	OUTFILE <- paste(outDIR,"/",cnvID,".loglik",sep="")
	cmd <- paste("perl /nfs/team143/ddd/bin/penncnv/detect_cnv.pl -validate -hmm /nfs/team143/ddd/bin/penncnv/lib/hh550.hmm",PFB,LIST,LOG,"-out",OUTFILE,REGION,sep=" ")
	system(cmd)
	if (file.exists(paste(outDIR,"/",cnvID,"-loglik.log",sep="")) ){
		cnCall <- try(read.table(paste(outDIR,"/",cnvID,"-loglik.log",sep=""),h=F,colClasses="character"))
		names(cnCall) <- c("STARTSNP","ENDSNP","CN0","CN1","CN2","CN3","CN4")
		row.names(cnCall) <- c(trioIDs[3],trioIDs[2],trioIDs[1])
		if (inherits(cnCall, 'try-error')){
			cnCall <- 0 
		} 
	}else{
		cnCall <- 0
	}
	return(cnCall)
}

test <- pennCnvLogLik(trioList, "test", currentRegion$snp[which(currentRegion$position==min(currentRegion$position))], currentRegion$snp[which(currentRegion$position==max(currentRegion$position))], CHR)


callACGHregion <- function(){
	stop("This function in need of SERIOUS development. Don't even attempt to use it!")
	# Input variables
	CHR <- 
	START <- 
	STOP <- 
	NPROBES <- 
	DAD
	MUM
	KID
	FID
	pfbFile
	lrrbafDir (DATAPATH)
	DIR
	
	# Read in aCGH and trio info
	acgh <- read.table("aCGH-results.txt", as.is=TRUE, header=TRUE, sep="\t")
	trioList <- scan("trio.txt",what=c("character","character","character"))
	
	# Copy trio files for relevant chromosome to working directory and generate trio list file and individual list file
	system("mkdir working")
	
	
	copyFile <- function(IID, DATAPATH){
		cmd <- paste("cp ",DATAPATH,"/",IID,"/FR-",CHR,".txt",sep="")
		system(cmd)
	}
	lapply(trioList, copyFile)
	# Write individual list
	list <- c(paste(outDir,"/",trioList[1],".txt",sep=""),paste(outDir,"/",trioList[2],".txt",sep=""),paste(outDir,"/",trioList[3],".txt",sep=""))
	write.table(list,file=paste(outDir,"/","indList.txt",sep=""),quote=F,row.names=F,col.names=F,sep = "\t")
	
	
	
	# Use the PFB to limit the range of SNPs we look at
	subsetPfb <- function(CHR, START, STOP, NPROBES, probePortion=0.5, DIR=getwd()){
		pfb <- read.table(paste("/lustre/scratch103/sanger/km5/ddd/snp-cnv/scotland/penncnv_in/pfb/scotland-PFB-",CHR,".txt",sep=""),h=F,colClasses=c("character","numeric","numeric","numeric"))
		pfb <- pfb[order(pfb[,3]), ]
		pfbSubset <- pfb[which(pfb[,3] >= START & pfb[,3] <= STOP),]
		nProbes <- length(pfbSubset[,1])
		# Check new start/stop positions against chromosome 
		pfbOut <- pfb[ (which(pfb[,1]==pfbSubset[1,1])-floor(NPROBES*probePortion)):(which(pfb[,1]==pfbSubset[length(pfbSubset[,1]),1])+ceiling(NPROBES*probePortion)), ]
		write.table(pfbOut, file=paste(DIR,"/pfb.txt",sep=""))
	}
	
	# Run PennCNV using this PFB file and chromosomal input files for trio as joint call
	pennCnvJoint <- function(DIR, CHR, START, STOP){
		PFB <- paste(" -pfb ",DIR,"/pfb.txt",sep="") 
		LIST <- paste(" --list ", DIR,"/","famList.txt ",sep="") 
		LOG <- paste(" -log ", DIR,"/",CHR,"-joint.log ",sep="") 
		OUTFILE <- paste(DIR,"/",CHR,"-",START,"-",STOP,".jointcnv",sep="")
		cmd <- paste("perl /nfs/team143/ddd/bin/penncnv/detect_cnv.pl -joint -hmm /nfs/team143/ddd/bin/penncnv/lib/hh550.hmm ",PFB,LIST,LOG," -out ",OUTFILE,sep="")
		system(cmd)
		if (file.exists(OUTFILE) ){
			cnCall <- try(read.table(OUTFILE,h=F,colClasses="character"))
			if (inherits(cnCall, 'try-error')){
				cnCall <- 0 
			} 
		}else{
			cnCall <- 0
		}
		return(cnCall)
	}
	

	# Use the PFB to limit the range of SNPs we look at
	subsetPfb <- function(CHR, START, STOP, NPROBES, probePortion=0.5, DIR=getwd()){
		pfb <- read.table(paste("/lustre/scratch103/sanger/km5/ddd/snp-cnv/scotland/penncnv_in/pfb/scotland-PFB-",CHR,".txt",sep=""),h=F,colClasses=c("character","numeric","numeric","numeric"))
		pfb <- pfb[order(pfb[,3]), ]
		pfbSubset <- pfb[which(pfb[,3] >= START & pfb[,3] <= STOP),]
		nProbes <- length(pfbSubset[,1])
		# Check new start/stop positions against chromosome 
		pfbOut <- pfb[ (which(pfb[,1]==pfbSubset[1,1])-floor(NPROBES*probePortion)):(which(pfb[,1]==pfbSubset[length(pfbSubset[,1]),1])+ceiling(NPROBES*probePortion)), ]
		write.table(pfbOut, file=paste(DIR,"/pfb.txt",sep=""))
	}
	
	# Run PennCNV using this PFB file and chromosomal input files for trio as joint call
	pennCnvJoint <- function(DIR, CHR, START, STOP){
		PFB <- paste(" -pfb ",DIR,"/pfb.txt",sep="") 
		LIST <- paste(" --list ", DIR,"/","famList.txt ",sep="") 
		LOG <- paste(" -log ", DIR,"/",CHR,"-joint.log ",sep="") 
		OUTFILE <- paste(DIR,"/",CHR,"-",START,"-",STOP,".jointcnv",sep="")
		cmd <- paste("perl /nfs/team143/ddd/bin/penncnv/detect_cnv.pl -joint -hmm /nfs/team143/ddd/bin/penncnv/lib/hh550.hmm ",PFB,LIST,LOG," -out ",OUTFILE,sep="")
		system(cmd)
		if (file.exists(OUTFILE) ){
			cnCall <- try(read.table(OUTFILE,h=F,colClasses="character"))
			if (inherits(cnCall, 'try-error')){
				cnCall <- 0 
			} 
		}else{
			cnCall <- 0
		}
		return(cnCall)
	}
	
	# Reformat output
	parseOutput <- function(cnCall){
		if (length(cnCall) == 1){
			return(0)
		}else{
			rowFilter <- function(cnCall){
				cnProbe <- as.numeric(strsplit(as.character(cnCall[2]),"=")[[1]][2])
				cnLength <- strsplit(strsplit(as.character(cnCall[3]),"=")[[1]][2],",")[[1]][1]
				if (length(strsplit(strsplit(as.character(cnCall[3]),"=")[[1]][2],",")[[1]]) > 1){
					for (i in 2:length(strsplit(strsplit(as.character(cnCall[3]),"=")[[1]][2],",")[[1]])){
						cnLength <- paste(cnLength,strsplit(strsplit(as.character(cnCall[3]),"=")[[1]][2],",")[[1]][i],sep="")
					}
				}
				cnLength <- as.numeric(cnLength)
				FILT <- FALSE
				if (filtType=="snp"){
					if (cnProbe >= filtThr){
						FILT <- TRUE
					}
				}else if (filtType=="length"){
					if (cnLength >= filtThr){
						FILT <- TRUE
					}
				}else if (filtType=="coverage"){
					if (cnProbe/cnLength >= filtThr){
						FILT <- TRUE
					}
				}
				if (FILT==TRUE){
					chr <- strsplit(strsplit(as.character(cnCall[1]),":")[[1]][1],"r")[[1]][2]
					start <- strsplit(strsplit(as.character(cnCall[1]),":")[[1]][2],"-")[[1]][1]
					stop <- strsplit(strsplit(as.character(cnCall[1]),":")[[1]][2],"-")[[1]][2]
					cn <- as.numeric(strsplit(as.character(cnCall[4]),"=")[[1]][2])
					startSnp <- strsplit(as.character(cnCall[6]),"=")[[1]][2] 
					stopSnp <- strsplit(as.character(cnCall[7]),"=")[[1]][2]
					temp <- cbind(chr,start,stop,cn,startSnp,stopSnp,cnProbe,cnLength,as.character(cnCall[8]))
					return(temp)
				}else{
					return(cbind(NA,NA,NA,NA,NA,NA,NA,NA,NA))
				}	
			}
			report <- apply(cnCall,1,rowFilter)
			filtList <- na.omit(data.frame(t(report)))
			names(filtList) <- c("CHR","START","STOP","COPYSTATE","STARTSNP","STOPSNP","NSNPS","LENGTH","IND")
			if (length(filtList[,1]) > 0){
				return(filtList)
			}else{
				return(0)
			}
		}
	}
	
	# Then run "validate" option to get logLik for each state
	# N.B. this has to be run on a list of INDIVIDUALS, not the family - so family information is only used to determine the CNV boundaries.
	pennCnvValidate <- function(){
		PFB <- paste(" -pfb ",DIR,"/pfb.txt",sep="") 
		LIST <- paste(" --list ", DIR,"/","indList.txt ",sep="") 
		LOG <- paste(" -log ", DIR,"/",CHR,"-joint.log ",sep="") 
		OUTFILE <- paste(DIR,"/",CHR,"-",START,"-",STOP,".jointcnv",sep="")
		cmd <- paste("perl /nfs/team143/ddd/bin/penncnv/detect_cnv.pl -joint -hmm /nfs/team143/ddd/bin/penncnv/lib/hh550.hmm ",PFB,LIST,LOG," -out ",OUTFILE,sep="")
		system(cmd)
		if (file.exists(OUTFILE) ){
			cnCall <- try(read.table(OUTFILE,h=F,colClasses="character"))
			if (inherits(cnCall, 'try-error')){
				cnCall <- 0 
			} 
		}else{
			cnCall <- 0
		}
		return(cnCall)
		
	}

	
# Write out PFB file for region
	write.table(currentRegion, file="region-pfb.txt", row.names=F, col.names=F, quote=F, sep = "\t")
	

	# Then compute likelihood ratio
	# Output summary info
	#	- CN state
	#	- number of probes
	#	- length
	#	- de novo status
		
	
	
	
}


