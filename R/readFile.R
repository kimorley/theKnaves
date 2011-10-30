# FUNCTIONS FOR READING GAPI OUTPUT FILES
# N.B. This is written for Illuminus output, not GenoSNP

# [1] Read normalised intensity files
#		This generates a list with two objects:
#		(a) Map file with SNP name, position, and two alleles
#		(b) Intensities (A and B) for each individual
readIntensity <- function(FILENAME){
	data <- read.table(FILENAME,h=T,colClasses="character",row.names=1)
	map <- data[1:2]
	intu <- data.frame(t(data[3:length(data)]))
	output <- list(map=map,intu=intu)
	return(output)
}

# [2] Read genotype call/confidence files
#		This generates a list with an object for each SNP.
#		Each object is a dataframe containing the genotype and call confidence for each individual
readCallConf <- function(FILENAME){
	data <- read.table(FILENAME,h=T,colClasses="character",row.names=1)
	temp <- data.frame(t(data))
	splitter <- function(x){
		call <- substring(x,1,2)
		conf <- substring(x,4,9)
		callConf <- data.frame(call,conf)
		ids <- sapply(row.names(temp),convertID)
		row.names(callConf) <- ids
		return(callConf)
	}
	output <- lapply(temp,splitter)
	return(output)
}
