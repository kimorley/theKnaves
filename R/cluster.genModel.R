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
					dd <- lapply(unique(as.character(cc$call)), subsetGenotypes, data = cc)
					names(dd) <- unique(as.character(cc$call))
					ee <- dd[!sapply(dd, is.null)]
					cc <- do.call(rbind, ee)
					genEllipse <- function(xx){
						if (length(xx[,1]) >= 5){
							ee = dataEllipse(xx$t,xx$r, levels=0.8, draw=F)
							ff = apply(xx,1,function(x) point.in.polygon(x[2], x[1], ee[,1], ee[,2]))
							gg = cbind(xx,ff)
							return( subset(xx, ff==1 ))
						}
					}
					ff <- lapply(dd,genEllipse)
					ff <- ff[!sapply(ff, is.null)]
					gg <- do.call(rbind, ff)
					if (length(gg[, 1]) == 0){
						print(paste("WARNING: No data for", SNP))
						cData <- list(m = NA, s = NA, a = alleles)
						save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""), compression_level = 9)
						return(cbind(SNP = SNP, FLAG = "N", N_SAMPLES = 0))
					}else if ( (var(gg$t, na.rm=T) %in% c(0, NA)) || (var(gg$r,na.rm=T) %in% c(0, NA)) ){
						print(paste("WARNING: No response for", SNP))
						cData <- list(m = NA, s = NA, a = alleles)
						save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""), compression_level = 9)
						report = rbind(report, cbind(SNP = SNP, FLAG = "F", N_SAMPLES =  length(gg[,1])))
					}else{	
						s <- data.frame(rbind(0, tapply(gg$t, as.character(gg$call), median), tapply(gg$t, as.character(gg$call), min), tapply(gg$t, as.character(gg$call), max), tapply(gg$t, as.character(gg$call), var), table(as.character(gg$call)), tapply(cc$t,as.character(cc$call), median), tapply(cc$t, as.character(cc$call), min), tapply(cc$t, as.character(cc$call), max), tapply(cc$t, as.character(cc$call), var), table(as.character(cc$call))))
						row.names(s) <- c("baf","median","min","max", "variance","n","median_a","min_a","max_a", "variance_a","n_a")
						if (alleles %in% names(s)) {
							s[1, which(names(s) == alleles)] <- 0.5
						}
						if (paste(substr(alleles, 2, 2), substr(alleles, 2, 2), sep = "") %in% names(s)) {
							s[1, which(names(s) == (paste(substr(alleles, 2, 2), substr(alleles, 2, 2), sep = "")))] <- 1
						}
						cData = list(a = alleles, s = s)
						if (length(s) == 3) {
							if ((s[2, which(names(s) == alleles)] < s[2, which(s[2, ] == min(s[2, ]))]) || (s[2, which(names(s) == alleles)] > s[2, which(s[2, ] == max(s[2, ]))])) {
								cData <- list(m = NA, s = NA, a = alleles)
								save(cData, file = paste(outDIR, "/clusterFile-", SNP, ".RData", sep = ""), compression_level = 9)
								return(cbind(SNP = SNP, FLAG = "F", N_SAMPLES = sum(s[6, ])))
							}else {
								temp <- lm(r ~ ns(t, knots = median(s[2,which(names(s) == alleles)]), Boundary.knots = c(s[2,which(s[2, ] == min(s[2, ]))], s[2, which(s[2, ] == max(s[2,]))])), data = gg, model = FALSE)
								clusterModel <- temp[match(c("coefficients", "rank", "qr", "terms"), names(temp))]
								class(clusterModel) <- class(temp)
								cData$m <- clusterModel
							}
						}else if (length(s) == 2) {
							temp <- lm(r ~ ns(t, Boundary.knots = c(s[2, which(s[2, ] == min(s[2, ]))], s[2, which(s[2, ] == max(s[2, ]))])), data = gg, model = FALSE)
							clusterModel <- temp[match(c("coefficients", "rank", "qr", "terms"), names(temp))]
							class(clusterModel) <- class(temp)
							cData$m <- clusterModel
						}else if (length(s) == 1) {
							temp <- lm(r ~ ns(t, Boundary.knots = c(min(gg$t), max(gg$t))), data = gg, model = FALSE)
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
