extractID <- function(x){
	ID <- length(strsplit(x,"_")[[1]])
	return(as.character(substring(strsplit(x,"_")[[1]][ID],1,nchar(strsplit(x,"_")[[1]][ID])-1)))
}

convertID <- function(x){
	ID <- length(strsplit(x,"_")[[1]])
	return(as.character(substring(strsplit(x,"_")[[1]][ID],1,nchar(strsplit(x,"_")[[1]][ID]))))
}

readIntensity <- function(FILENAME){
	data <- read.table(FILENAME,h=T,colClasses="character",row.names=1)
	map <- data[1:2]
	intu <- data.frame(t(data[3:length(data)]))
	output <- list(map=map,intu=intu)
	return(output)
}

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
