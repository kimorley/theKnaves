cluster.genModel <- function(gapiGTU, gapiINTU, outDIR, confidence = 0.99, splitIDs = TRUE){
	snp <- names(eval(as.list(sys.call(-1))[[2]]))[as.numeric(gsub("[^0-9]","", deparse(substitute(gapiGTU))))]
	SNP <- unlist(strsplit(snp, "X"))[length(unlist(strsplit(snp, "X")))]
	GTU <- gapiGTU
	INTall <- readSnpIntu(SNP, gapiINTU)
	INT <- INTall$intu
	if (is.na(sum(as.numeric(unlist(INT))))){
		print(paste("WARNING: No intensity data available for SNP", SNP, "in", gapiINTU, sep = " "))
		cData <- list(m = NA, s = NA)
		save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""), compression_level = 9)
		return(cbind(SNP = SNP, FLAG = "N", N_SAMPLES = 0))
	}else{
		alleles <- unlist(INTall$map$Alleles)
		# Reformat and merge genotype calls and intensity files
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
			t <- atan2(as.numeric(as.character(B)), as.numeric(as.character(A)))/(pi/2) # Calculate angle (theta) from signal intensities (standardised)
			r <- as.numeric(as.character(A)) + as.numeric(as.character(B)) # This is the Manhattan distance, which is what Illumina use (*should* be Euclidean)
			names(r) <- namesA
			names(t) <- namesA
			if (sum(row.names(GTU) == namesA) == length(namesA)){
				data = cbind(data.frame(r, t), GTU)
				rm(index, use, A, B, namesA, namesB, GTU, INTall, INT)
				# Process and return R object for use with generating LRR/BAF
				report = NULL
				if (is.null(data)) {
					cData <- list(m = NA, s = NA, a = alleles)
					save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""))
					report = rbind(report, (cbind(SNP = SNP, FLAG = "N", N_SAMPLES = 0)))
				}else{
					cc <- data[which(as.numeric(as.character(data$conf)) >= confidence & as.character(data$call) != "NN"), ]	# Remove data for which genotype could not be called with high confidence
					subsetGenotypes <- function(genotype, data) {
						xx <- subset(data, call == genotype)
						if (length(xx[,1]) >= 5){	# Remove genotypes with less than five samples in the cluster
							return(xx)
						}
					}
					dd <- lapply(unique(as.character(cc$call)), subsetGenotypes, data = cc)
					names(dd) <- unique(as.character(cc$call))
					cc <- do.call(rbind, dd[!sapply(dd, is.null)])	# This now contains data for clusters with > 5 points
					genEllipse <- function(xx){
						if (length(xx[,1]) >= 5){
							ee <- try(dataEllipse(xx$t,xx$r, levels=0.8, draw=F),silent=T)	# Draw ellipse for each genotype group that contains 80% of points
							if (inherits(data, 'try-error')){
								return(xx)
							}else{
								ff = apply(xx,1,function(x) point.in.polygon(x[2], x[1], ee[,1], ee[,2]))	# Check which observations fall within the ellipse...	
								gg = cbind(xx,ff)
								return( subset(xx, ff==1 ))		# ...and keep only those within the ellipse
							}
						}
					}
					ff <- lapply(dd,genEllipse)
					gg <- do.call(rbind, ff[!sapply(ff, is.null)])
					if (length(gg[, 1]) == 0){	# If removed all the data
						print(paste("WARNING: No data for", SNP))
						cData <- list(m = NA, s = NA, a = alleles)
						save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""), compression_level = 9)
						return(cbind(SNP = SNP, FLAG = "N", N_SAMPLES = 0))
					}else if ( (var(gg$t, na.rm=T) %in% c(0, NA)) || (var(gg$r,na.rm=T) %in% c(0, NA)) ){	# If we have no samples with intensity
						print(paste("WARNING: No response for", SNP))
						cData <- list(m = NA, s = NA, a = alleles)
						save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""), compression_level = 9)
						report = rbind(report, cbind(SNP = SNP, FLAG = "F", N_SAMPLES =  length(gg[,1])))
					}else{	
						summaryTable <- function(xx){
							s <- data.frame(rbind(0, tapply(xx$t, as.character(xx$call), median), tapply(xx$t, as.character(xx$call), min), tapply(xx$t, as.character(xx$call), max), tapply(xx$t, as.character(xx$call), var), table(as.character(xx$call))))
							row.names(s) <- c("baf","median","min","max", "variance","n")
							if (alleles %in% names(s)) {
								s[1, which(names(s) == alleles)] <- 0.5
							}
							if (paste(substr(alleles, 2, 2), substr(alleles, 2, 2), sep = "") %in% names(s)) {
								s[1, which(names(s) == (paste(substr(alleles, 2, 2), substr(alleles, 2, 2), sep = "")))] <- 1
							}
							return(s)
						}
						s = summaryTable(gg)
						cData = list(a = alleles, s = s, ss = summaryTable(cc))
						if (length(s) == 3) {
							if ((s[2, which(names(s) == alleles)] < s[2, which(s[2, ] == min(s[2, ]))]) || (s[2, which(names(s) == alleles)] > s[2, which(s[2, ] == max(s[2, ]))])) {	# Heterozygote is not between the two homozygotes - ignore this SNP
								cData <- list(m = NA, s = NA, a = alleles)
								save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""), compression_level = 9)
								return(cbind(SNP = SNP, FLAG = "F", N_SAMPLES = sum(s[6, ])))
							}else {
								temp <- lm(r ~ ns(t, knots = median(s[2,which(names(s) == alleles)]), Boundary.knots = c(s[2,which(s[2, ] == min(s[2, ]))], s[2, which(s[2, ] == max(s[2,]))])), data = gg, model = FALSE)
								clusterModel <- temp[match(c("coefficients", "rank", "qr", "terms"), names(temp))]
								class(clusterModel) <- class(temp)
								cData$m <- clusterModel
								rm(temp, clusterModel)
							}
						}else if (length(s) == 2) {
							temp <- lm(r ~ ns(t, Boundary.knots = c(s[2, which(s[2, ] == min(s[2, ]))], s[2, which(s[2, ] == max(s[2, ]))])), data = gg, model = FALSE)
							clusterModel <- temp[match(c("coefficients", "rank", "qr", "terms"), names(temp))]
							class(clusterModel) <- class(temp)
							cData$m <- clusterModel
							rm(temp, clusterModel)
						}else if (length(s) == 1) {
							temp <- lm(r ~ ns(t, Boundary.knots = c(min(gg$t), max(gg$t))), data = gg, model = FALSE)
							clusterModel <- temp[match(c("coefficients", "rank", "qr", "terms"), names(temp))]
							class(clusterModel) <- class(temp)
							cData$m <- clusterModel
							rm(temp, clusterModel)
						}
						rm(cc, dd, ff, gg)
						save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""), compression_level = 9)	# Cluster file object
						report = rbind(report, cbind(SNP = SNP, FLAG = "P", N_SAMPLES = sum(s[6, ])))	# Report on SNP
					}
				}
				return(report)
			}else {
				stop(paste("Different samples in intensity and genotype files for SNP", SNP, sep = " "))
			}
		}else {
			stop(paste("IDs for A allele data and B allele data do not match or have different order for SNP", SNP, sep = " "))
		}
	}
}
