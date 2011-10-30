# Split ID on underscore (_) and then take last element of vector as ID
extractID <- function(x){
	ID <- length(strsplit(x,"_")[[1]])
	return(as.character(substring(strsplit(x,"_")[[1]][ID],1,nchar(strsplit(x,"_")[[1]][ID])-1)))
}

# Split on underscore and return vector of elements
convertID <- function(x){
	ID <- length(strsplit(x,"_")[[1]])
	return(as.character(substring(strsplit(x,"_")[[1]][ID],1,nchar(strsplit(x,"_")[[1]][ID]))))
}