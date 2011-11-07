# FUNCTIONS FOR READING GAPI OUTPUT FILES
# N.B. This is written for Illuminus output, not GenoSNP

# [1] Read normalised intensity files
#	This generates a list with two objects:
#	(a) Map file with SNP name, position, and two alleles
#	(b) Intensities (A and B) for each individual
readRawIntu <- function(FILENAME){
	# Read and check
	data <- read.table(FILENAME,h=T,colClasses="character",row.names=1)
	if (sum(names(data)[1:2] == c("Coor", "Alleles")) != 2){
		stop("Coor and Alleles expected as headers in first two columns after SNP.")
	}
	map <- subset(data, select=c(Coor,Alleles))
	intu <- data.frame(t(data[3:length(data)]))
	if (length(row.names(intu))%%2 == 1){
		stop("Number of intensity values is not divisible by two; expecting two values (A and B) per sample.")
	}
	output <- list(map=map,intu=intu)
	return(output)
}

# [2] Read genotype call/confidence files
#	This generates a list with an object for each SNP.
#	Each object is a dataframe containing the genotype and call confidence for each individual
#	Default is to split IDs on underscore and take last two elements; select splitIDs=FALSE to avoid
readRawGtu <- function(FILENAME,splitIDs=TRUE){
	data <- read.table(FILENAME,h=T,colClasses="character",row.names=1)
	temp <- data.frame(t(data))
	splitter <- function(x){
		call <- substring(x,1,2)
		conf <- substring(x,4,9)
		callConf <- data.frame(call,conf)
		if (splitIDs){
			ids <- sapply(row.names(temp),idSingle)
			row.names(callConf) <- ids
			return(callConf)
		}else{
			row.names(callConf) <- row.names(temp)
			return(callConf)
		}
	}
	output <- lapply(temp,splitter)
	return(output)
}

# [3] Read normalised intensity file for specific SNP (for canonical cluster generation)
readSnpIntu <- function(SNP,clusterINTU){
	cmd <- paste( "awk '{ if (NR==1 || $1 == \"" ,SNP, "\") print $0 }' " ,clusterINTU, sep="")
	data <- read.table(pipe(cmd),h=T,row.names=1, colClasses="character")
	if (sum(names(data)[1:2] == c("Coor", "Alleles")) != 2){
		stop("Coor and Alleles expected as headers in first two columns after SNP.")
	}
	map <- subset(data, select=c(Coor,Alleles))
	intu <- data.frame(t(data[3:length(data)]))
	if (length(row.names(intu))%%2 == 1){
		stop("Number of intensity values is not divisible by two; expecting two values (A and B) per sample.")
	}
	output <- list(map=map,intu=intu)
	return(output)
}

# [4] Read genotype call/confidence for specific SNP (for canonical cluster generation)
readSnpGtu <- function(SNP,clusterGTU,splitIDs=T){
	cmd <- paste( "awk '{ if (NR==1 || $1 == \"" ,SNP, "\") print $0 }' " ,clusterGTU, sep="")
	data <- read.table(pipe(cmd),h=T,row.names=1,colClasses="character")
	temp <- data.frame(t(data))
	splitter <- function(x){
		call <- substring(x,1,2)
		conf <- substring(x,4,9)
		callConf <- data.frame(call,conf)
		if (splitIDs){
			ids <- sapply(row.names(temp),idSingle)
			row.names(callConf) <- ids
			return(callConf)
		}else{
			row.names(callConf) <- row.names(temp)
			return(callConf)
		}
	}
	output <- lapply(temp,splitter)
	return(output)
}
		