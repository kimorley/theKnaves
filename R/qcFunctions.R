summaryMarkers <- function(x){
	X1 <- x[1:(length(x)-1)]
	X2 <- x[2:length(x)]
	sumMark <- c(length(x),median(X2-X1))
	return(sumMark)
}

summaryBAF <- function(x){
	sumBaf <- c(round(median(x[which(x > 0.25 & x <0.75)]),3), round(mean(x[which(x > 0.25 & x <0.75)]),3), round(sd(x[which(x > 0.25 & x <0.75)]),3), round(mad(x,na.rm=T)/1.4826,3)) # BAF is NOT normally distributed so remove constant
	return(sumBaf)
}

summaryLRR <- function(x){
	sumLrr <- c(round(median(x[which(x > -2 & x <2)]),3), round(mean(x[which(x > -2 & x <2)]),3), round(sd(x[which(x > -2 & x <2)]),3), round(mad(x,na.rm=T)/1.4826,3)) # LRR is normal-ish but removed anyway  
	return(sumLrr)
} 
