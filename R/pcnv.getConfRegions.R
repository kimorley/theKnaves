# Function for extracting and writing list of regions for generating confidence estimates
pcnv.getConfRegions <- function(cnCall){
	# This assumes the input is the standard PennCNV output file
	startSnp <- unlist(strsplit(as.character(cnCall[6]),"="))[2]
	endSnp <- unlist(strsplit(as.character(cnCall[7]),"="))[2]
	region <- cbind(cnCall[1],startSnp,endSnp,0.01,0.01)
}
