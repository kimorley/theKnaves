
convertIntGeno <- function(intFile,gtuFile){
	# Read in intensity and genotype files from GAPI
	gapiInt <- readIntensity(intFile)
	gapiGtu <- readCallConf(gtuFile)		
	# Function for converting intensities to R and theta, and merging with genotype call/confidence
	convert <- function(INT,GTU){
		# Split intensity data into A and B allele sets
		index <- rep(c(0,1),(length(as.numeric(INT))/2))
		use <- index==0
		A <- INT[use]
		namesA <- sapply(names(A),extractID)
		use <- index==1
		B <- INT[use]
		namesB <- sapply(names(B),extractID)
		if (sum(namesA==namesB)==length(INT)/2){
			# Calculate angle (theta) from signal intensities
			t <- atan2(as.numeric(as.character(B)),as.numeric(as.character(A))) / (pi/2) # Make sure values range from 0 to 1
			# Calculate ray (R) from signal intensities
			# r <- sqrt(as.numeric(as.character(A))^2 + as.numeric(as.character(B))^2)	# This is the Euclididean distance i.e. the hypotenuse of the triangle for conversion of Cartesian to polar BUT Illumina DOESN'T use this
			r <- as.numeric(as.character(A)) + as.numeric(as.character(B))	# This is the Manhattan distance, which is what Illumina use
			names(r) <- namesA
			names(t) <- namesA
			if (sum(row.names(GTU) == namesA) == length(namesA)){
				results <- cbind(data.frame(r,t), GTU)
				return(results)
			}else{
				stop("Different samples in intensity and genotype files.")
			}
		}else{
			stop("IDs for A allele data and B allele data do not match or have different order.")
		}
	}
	# Now apply the function to the two lists of intensity and genotype call data
	if ( sum(names(gapiInt$intu) == names(gapiGtu))/length(names(gapiInt$intu) == names(gapiGtu)) == 1){
		intu <- gapiInt$intu
		output <- mapply(convert,intu,gapiGtu,SIMPLIFY=F)
		allData <- list(map=gapiInt$map,data=output)
		return(output)
	}else{
		stop("Different SNPs in intensity and genotype call files.")
	}
}

