convertIntensity <- function(LIST){
	SNP <- LIST$intu
	convert <- function(INT){
		index <- rep(c(0,1),(length(as.numeric(INT))/2))
		use <- index==0
		A <- INT[use]
		namesA <- sapply(names(A),extractID)
		use <- index==1
		B <- INT[use]
		namesB <- sapply(names(B),extractID)
		if (sum(namesA==namesB)==length(INT)/2){
			theta <- atan2(as.numeric(as.character(B)),as.numeric(as.character(A))) / (pi/2)
			R <- as.numeric(as.character(A)) + as.numeric(as.character(B))
			names(R) <- namesA
			names(theta) <- namesA
			results <- data.frame(R,theta)
			return(results)
		}else{
			stop("IDs for A allele data and B allele data do not match or have different order.")
		}
	}
	output <- lapply(SNP,convert)
	return(output)
}