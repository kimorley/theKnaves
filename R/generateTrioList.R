
generateTrioList <- function(FAM,outDir=getwd()){
	# Setup
	setwd(outDir)
	famFile <- read.table(FAM,h=F,colClasses="character")
	# Check
	if (length(famFile[1,]) < 4){
		stop("Need minimum of four columns: Family ID, Child ID, Father ID, Mother ID")
	}
	# Execute
	kids <- famFile[which(famFile$V3 != 0 | famFile$V4 != 0),]
	writeTrio <- function(ids){
		trioList <- c(strsplit(ids[3],"_")[[1]][length(strsplit(ids[3],"_")[[1]])], strsplit(ids[4],"_")[[1]][length(strsplit(ids[4],"_")[[1]])], strsplit(ids[2],"_")[[1]][length(strsplit(ids[2],"_")[[1]])])
		cmd <- paste("mkdir",ids[1],sep=" ")
		system(cmd)
		write.table(trioList,file=paste(outDir,ids[1],"trio.list",sep="/"),quote=F,row.names=F,col.names=F)
	}
	apply(kids,1,writeTrio)
}
