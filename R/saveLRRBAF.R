
saveLRRBAF <- function(data, CHR, outDIR=getwd()){
	if (missing(data)){
		stop("Must provide data argument.")
	}
	if (missing(CHR)){
		stop("Must provide chromosome argument.")
	}
	# Combine data across all list objects
	combo <- do.call(cbind,data)
	# Transpose data
	tCombo <- data.frame(t(combo))
	# Select LRR
	lrrUse <- row.names(tCombo)%in%grep("lrr",row.names(tCombo),value=T)
	lrrData <- tCombo[lrrUse,]
	# Replace row names with rsID only
	row.names(lrrData) <- unlist(strsplit(row.names(lrrData),"\\."))[seq(1,length(row.names(lrrData))*2,2)]
	# Select BAF
	bafUse <- row.names(tCombo)%in%grep("baf",row.names(tCombo),value=T)
	bafData <- tCombo[bafUse,]
	row.names(bafData) <- unlist(strsplit(row.names(bafData),"\\."))[seq(1,length(row.names(bafData))*2,2)]
	# Loop through all individuals in data frame and print out their LRR/BAF
	for (i in 1:length((names(tCombo)))){
		output <- data.frame(row.names(lrrData),subset(bafData,select=names(tCombo[i])),subset(lrrData,select=names(tCombo[i])))
		ID <- names(output)[2]
		names(output) <- c("Name",paste(ID,"-",CHR,"A.B Allele Freq",sep=""),paste(ID,"-",CHR,"A.Log R Ratio",sep=""))
		dir.create(file.path(outDIR, ID), showWarnings = FALSE)
		write.table(output,file=paste(outDIR,"/",ID,"/",ID,"-",CHR,".txt",sep=""),quote=F,row.names=F,sep = "\t")
	}
}

	