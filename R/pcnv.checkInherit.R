pcnv.checkInherit <- function(x){
	x <- cbind(x,cnvid=paste(x$chr,x$start_bp,x$stop_bp,sep=":"))
	o <- subset(x, id=="offspring")	# Offspring calls
	if (length(o[,1])==0){
		return(NULL)
	}else{
		t <- subset(x, id!="offspring")	# Parental calls
		use <- o$cnvid %in% t$cnvid
		if (sum(!use)==0){	# No deNovo variants
			p <- subset(x, id=="father")
			m <- subset(x, id=="mother")
			i <- subset(x, id=="offspring")
			use_p <- i$cnvid %in% p$cnvid
			use_m <- i$cnvid %in% m$cnvid
			if ( (sum(use_p)!=0) & (sum(use_m)!=0) ){	# Both parents
				ip <- data.frame(i[use_p, ], dnv="P")
				im <- data.frame(i[use_m,], dnv="M")
				it <- rbind(ip,im)
				it$dnv <- as.character(it$dnv)
				a <- count(data.frame(it$cnvid))
				a <- as.vector(unlist(subset(a, freq==2, select=c(it.cnvid))))
				for (cnv_i in 1:length(a)){
					it$dnv[which(as.character(it$cnvid)==a[cnv_i])] <- "B"
				}
				it <- count(it)
				it$freq <- NULL
				return(it)
			}else if ( (sum(use_p)!=0) & (sum(use_m)==0) ){	# Paternal only
				ip <- data.frame(i[use_p,], dnv="P")
				ip$dnv <- as.character(ip$dnv)
				return(ip)
			}else if ( (sum(use_p)==0) & (sum(use_m)!=0) ){	# Maternal only
				im <- data.frame(i[use_m,], dnv="M")
				im$dnv <- as.character(im$dnv)
				return(im)
			}
		}else if (sum(use)==0){	# No inherited variants
			d <- data.frame(o[!use,], dnv="Y")
			d$dnv <- as.character(d$dnv)						
			return(d)
		}else{	# Both types
			d <- data.frame(o[!use,], dnv="Y")
			d$dnv <- as.character(d$dnv)
			p <- subset(x, id=="father")
			m <- subset(x, id=="mother")
			i <- subset(x, id=="offspring")
			use_p <- i$cnvid %in% p$cnvid
			use_m <- i$cnvid %in% m$cnvid
			if ( (sum(use_p)!=0) & (sum(use_m)!=0) ){	# Both
				ip <- data.frame(i[use_p,], dnv="P")
				im <- data.frame(i[use_m,], dnv="M")
				it <- rbind(ip,im)
				it$dnv <- as.character(it$dnv)
				a <- count(data.frame(it$cnvid))
				a <- as.vector(unlist(subset(a, freq==2, select=c(it.cnvid))))
				for (cnv_i in 1:length(a)){
					it$dnv[which(as.character(it$cnvid)==a[cnv_i])] <- "B"
				}
				it <- count(it)
				it$freq <- NULL
				all <- rbind(d,it)
				return(all)
			}else if ( (sum(use_p)!=0) & (sum(use_m)==0) ){	# Paternal only
				ip <- data.frame(i[use_p,], dnv="P")
				ip$dnv <- as.character(ip$dnv)
				all <- rbind(d,ip)
				return(all)
			}else if ( (sum(use_p)==0) & (sum(use_m)!=0) ){	# Maternal only
				im <- data.frame(i[use_m,], dnv="M")
				im$dnv <- as.character(im$dnv)
				all <- rbind(d,im)
				return(all)
			}
		}
	}
}
