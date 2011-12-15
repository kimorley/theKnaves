cluster.genModel <- function(gapiGTU, gapiINTU, outDIR, confidence = 0.99, splitIDs = TRUE){
	snp <- names(eval(as.list(sys.call(-1))[[2]]))[as.numeric(gsub("[^0-9]","", deparse(substitute(gapiGTU))))]
	SNP <- unlist(strsplit(snp, "X"))[length(unlist(strsplit(snp, "X")))]
	GTU <- gapiGTU
	INTall <- readSnpIntu(SNP, gapiINTU)
	INT <- INTall$intu
	if (is.na(sum(as.numeric(unlist(INT))))){
		print(paste("WARNING: No intensity data available for SNP", SNP, "in", gapiINTU, sep = " "))
		gsCluster <- clusterSnp(x = NA, SNP, DIR = outDIR)
	}else{
		alleles <- unlist(INTall$map$Alleles)
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
				data = cbind(data.frame(r, t), GTU)
				report = NULL
				if (is.null(data)) {
					cData <- list(m = NA, s = NA, a = alleles)
					save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""))
					report = rbind(report, (cbind(SNP = SNP, FLAG = "N", N_SAMPLES = 0)))
				}else{
					cc <- data[which(as.numeric(as.character(data$conf)) >= confidence & as.character(data$call) != "NN"), ]
					subsetGenotypes <- function(genotype, data) {
						xx <- subset(data, call == genotype)
						if (length(xx[,1]) >= 5){	# Remove genotypes with less than five samples in the cluster
							return(xx)
						}
					}
					cData <- lapply(unique(as.character(cc$call)), subsetGenotypes, data = cc)
					names(cData) <- unique(as.character(cc$call))
					cData <- cData[!sapply(cData, is.null)]
					cc <- do.call(rbind, cData)
					if (length(cc[, 1]) == 0) {
						print(paste("WARNING: No data for", SNP))
						cData <- list(m = NA, s = NA, a = alleles)
						save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""), compression_level = 9)
						return(cbind(SNP = SNP, FLAG = "N", N_SAMPLES = 0))
					}else if ( var(cc$t) %in% c(0, NA) || var(cc$r) %in% c(0, NA) ){
						print(paste("WARNING: No response for", SNP))
						cData <- list(m = NA, s = NA, a = alleles)
						save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""), compression_level = 9)
						report = rbind(report, cbind(SNP = SNP, FLAG = "F", N_SAMPLES =  length(cc[,1])))
					}else{
						s <- data.frame(rbind(0, tapply(cc$t, as.character(cc$call), median), tapply(cc$t, as.character(cc$call), min), tapply(cc$t, as.character(cc$call), max), tapply(cc$t, as.character(cc$call), var), table(as.character(cc$call))))
						row.names(s) <- c("baf","median","min","max", "variance","n")
						if (alleles %in% names(s)) {
							s[1, which(names(s) == alleles)] <- 0.5
						}
						if (paste(substr(alleles, 2, 2), substr(alleles, 2, 2), sep = "") %in% names(s)) {
							s[1, which(names(s) == (paste(substr(alleles, 2, 2), substr(alleles, 2, 2), sep = "")))] <- 1
						}
						cData$a <- alleles
						cData$s <- s
						if (length(s) == 3) {
							if ((s[2, which(names(s) == alleles)] < s[2, which(s[2, ] == min(s[2, ]))]) || (s[2, which(names(s) == alleles)] > s[2, which(s[2, ] == max(s[2, ]))])) {
								cData <- list(m = NA, s = NA, a = alleles)
								save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""), compression_level = 9)
								return(cbind(SNP = SNP, FLAG = "F", N_SAMPLES = sum(s[6, ])))
							}else {
								temp <- lm(r ~ ns(t, knots = median(s[2,which(names(s) == alleles)]), Boundary.knots = c(s[2,which(s[2, ] == min(s[2, ]))], s[2, which(s[2, ] == max(s[2,]))])), data = cc, model = FALSE)
								clusterModel <- temp[match(c("coefficients", "rank", "qr", "terms"), names(temp))]
								class(clusterModel) <- class(temp)
								cData$m <- clusterModel
							}
						}else if (length(s) == 2) {
							temp <- lm(r ~ ns(t, Boundary.knots = c(s[2, which(s[2, ] == min(s[2, ]))], s[2, which(s[2, ] == max(s[2, ]))])), data = cc, model = FALSE)
							clusterModel <- temp[match(c("coefficients", "rank", "qr", "terms"), names(temp))]
							class(clusterModel) <- class(temp)
							cData$m <- clusterModel
						}else if (length(s) == 1) {
							temp <- lm(r ~ ns(t, Boundary.knots = c(min(cc$t), max(cc$t))), data = cc, model = FALSE)
							clusterModel <- temp[match(c("coefficients", "rank", "qr", "terms"), names(temp))]
							class(clusterModel) <- class(temp)
							cData$m <- clusterModel
						}
						save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""), compression_level = 9)
						report = rbind(report, cbind(SNP = SNP, FLAG = "P", N_SAMPLES = sum(s[6, ])))
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
