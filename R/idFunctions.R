# Get ID for intensity files where there are two measurements (A & B) for each individual
# Split ID on underscore (_) and then take last two elements of vector as ID
idDouble <- function(x){
	ID <- length(strsplit(x,"_")[[1]])
	return( paste(unlist(strsplit(x,"_"))[ID-1],substring(unlist(strsplit(x,"_"))[ID],1,nchar(unlist(strsplit(x,"_"))[ID])-1),sep="_") )
}

# Get ID for files where there is a single measurement for an invidual (no suffix on ID)
# Split ID on underscore (_) and then take last two elements of vector as ID
idSingle <- function(x){
	ID <- length(strsplit(x,"_")[[1]])
	return(paste(unlist(strsplit(x,"_"))[ID-1],unlist(strsplit(x,"_"))[ID],sep="_"))
}