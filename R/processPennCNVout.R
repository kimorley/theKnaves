filterCalls <- function(cnList,minProbe=1,minLength=1,minCoverage=1){
	# Check
	if (minProbe==1 & minLength==1 & minCoverage==1){
		stop("Please provide filtering criteria, either minimum number of probes and minimum length, or minimum coverage.")
	}
	# Setup
	filtType <- c()
	filtThr <- c()
	if (minProbe > 1){
		filtType <- "snp"
		filtThr <- minProbe
		out <- paste("Filtering using criteria: minimum SNPs = ",minProbe,sep="")
	}else if (minLength > 1){
		filtType <- "length"
		filtThr <- minLength
		out <- paste("Filtering using criteria: minimum length = ",minLength,sep="")
	}else{
		filtType <- "coverage"
		filtThr <- minCoverage
		out <- paste("Filtering using criteria: minimum coverage = ",minCoverage,sep="")
	}
	print(out)
	# Execute
	filter <- function(cnCall){
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
	filtCalls <- lapply(cnList,filter)
}

writeFilteredList <- function(filtList){
	allChr <- data.frame()
	for (i in 1:length(names(filtList))){
		if (length(filtList[[i]]) != 1)
		allChr <- rbind(allChr,filtList[[i]])
	}
	return(allChr)
}

#splitInheritance <- function(data){
#	if (names(data) %in% c("CHR","START","STOP","COPYSTATE","STARTSNP","STOPSNP","NSNPS","LENGTH","IND")){
#		offspring <- subset(data, ID == "offspring")
#		parents <- subset(data, ID != "offspring")
#		deNovoTest <- offspring$V1 %in% parents$V1
#		deNovoCNV <- offspring[!deNovoTest,]
#		inheritCNV <- offspring[deNovoTest,]
#	}else{
#		stop("Incorrect headers.  Run filtering functions first.")
#	}
#}

