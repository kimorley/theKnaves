processTrio <- function(trioList,inDir=getwd(),outDir=getwd(),chrStart=1,chrStop=22,writeOut=FALSE,returnData=TRUE,writeQC=FALSE){
	## Check
	if (length(trioList)!=3){
		stop("Trio list should consist of three individuals.")
	}
	if (!is.numeric(chrStart) || !is.numeric(chrStop)){
		stop("Start and stop chromosome should be numeric; recode X, Y, MT data as numbers.")
	}
	## Execute
	# Process final report for each individual (stored in list of lists)
	trioData <- mapply(processFinalReport,trioList,inDir=inDir,outDir=outDir,chrStart=chrStart,chrStop=chrStop,writeOut=writeQC)
	# If true, write out input files for PennCNV (data files and list files)
	if (writeOut){
		print("Writing input files for PennCNV")
		writePennCNVTrio(trioData,outDir=outDir)
	}
	if (returnData){
		return(trioData)
	}
}
