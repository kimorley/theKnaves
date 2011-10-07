mergeGenoThetaR <- function(GTU,RT){
	if (sum(names(RT) == names(GTU))==length(GTU)){
		mergeList <- c()
		for (i in 1:length(GTU)){
			if (names(RT)[i]==names(GTU)[i]){
				if (sum(row.names(RT[[i]])==row.names(GTU[[i]]))==length(GTU[[i]][[1]])){
					mergeList[[names(RT)[i]]] <- cbind(RT[[i]],GTU[[i]])
				}else{
					stop("Different sample IDs in intensity and call/confidence files.")
				}
			}else{
				stop("Different order or name format of markers in intensity and call/confidence files.")
			}
		}
		return(mergeList)
	}else{
		stop("Intensity and call/confidence files do not contain the same markers.")
	}
}