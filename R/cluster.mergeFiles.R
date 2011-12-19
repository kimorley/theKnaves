cluster.mergeFiles <- function(SNP, gapiGTU, gapiINTU, confidence = 0.99, splitIDs = TRUE){
	print(SNP)
	GTU <- readSnpGtu(SNP, gapiGTU)
	GTU = GTU[[1]]
	int <- readSnpIntu(SNP, gapiINTU)
	INT <- int$intu
	if (!is.null(INT) && !is.null(GTU)){
		alleles <- unlist(int$map$Alleles)
		index <- rep(c(0, 1), (length(row.names(INT))/2))
		use <- index == 0
		A <- INT[use, ]
		use <- index == 1
		B <- INT[use, ]
		namesA <- names(A)
		namesB <- names(B)
		if (splitIDs) {
			namesA <- sapply(names(A), idDouble)
			namesB <- sapply(names(B), idDouble)
		}
		if (sum(namesA == namesB) == length(row.names(INT))/2){
			t <- atan2(as.numeric(as.character(B)), as.numeric(as.character(A)))/(pi/2)
			r <- as.numeric(as.character(A)) + as.numeric(as.character(B))
			names(r) <- namesA
			names(t) <- namesA
			if (sum(row.names(GTU) == namesA) == length(namesA)){
				return(data = cbind(data.frame(r, t), GTU))
			}else {
				stop(paste("Different samples in intensity and genotype files for SNP", SNP, sep = " "))
			}
		}else {
			stop(paste("IDs for A allele data and B allele data do not match or have different order for SNP", SNP, sep = " "))
		}
	}
}
