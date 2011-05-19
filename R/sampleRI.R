# ChangeLog:
# 07.10.2008: the function was change to support the new Library format.

sampleRI <-
function(samples, Lib, r_thres=0.95, columns = c("SPECTRUM", "RETENTION_TIME_INDEX", "RETENTION_TIME"),
 method = "dayNorm", minPairObs = 5, showProgressBar = FALSE,
 makeReport = FALSE, pdfFile = "medianLibRep.pdf"){

	my.files  <- RIfiles(samples)
	Names     <- sampleNames(samples)
	refLib    <- refLib(Lib, w = 2, sel = TRUE)
	libId     <- libId(Lib, sel = TRUE)
	
	RES       <- FindPeaks(my.files, refLib, columns, showProgressBar)

	if(makeReport == TRUE)
	 	plotAllRIdev(Lib, RES, pdfFile)

	# you can't set this paremeter lower than 5.
	minPairObs <- max(minPairObs, 5)

	resInt    <- Intensity(RES)
	resRI     <- retIndex(RES)

	# normalise intensities using "method"
	res <- switch(method,
						dayNorm = dayNorm(samples, resInt),
						medianNorm = medianNorm(samples, resInt),
						none = resInt)

	res_log <- lapply(res,log2)

	met_cor <- list()
	options(warn = -1)
	
	if(showProgressBar) 
		pb <- ProgressBar(title="Correlating Masses...", label="File in processing...")
  
	for(i in 1:length(Lib)) {

  	if(showProgressBar)  
		  setProgressBar(pb, value=i/length(Lib),
						title=sprintf("Correlating Masses (%d %%)",round(100*i/length(Lib))),
						label=sprintf("Metabolite %d",i))

		  x <- which(libId == i)

      # don't perform calculation with 1 selective mass
      if (length(x) == 1) {
       	met_cor[[i]] <- 1
       	next
      }

      tmp <- cor(t(res_log[[i]]), use="pair")
      # this counts the number of pair values that were used to calculate the correlation
      # coefficient of every member of "tmp" and set to 0 the pairs with less than minPairObs
      tmp[is.finite(res_log[[i]]) %*% is.finite(t(res_log[[i]])) < minPairObs] <- 0

      # assume that the correlation of a metabolite with itself is always 1
      diag(tmp) <- 1

      tmp.max <- which.max(apply(tmp,1,function(x){ sum(x > r_thres, na.rm=T)}))
      tmp.sel <- tmp[tmp.max,]
      met_cor[[i]] <- which(tmp.sel > r_thres)
  }
  if(showProgressBar)  
		close(pb)  
	options(warn = 0)
	cor_RI <- matrix(ncol=length(my.files),nrow=length(Lib))
  colnames(cor_RI) <- Names
  rownames(cor_RI) <- rownames(libData(Lib))
  
	apply2 <- function(X, MARGIN, FUN, ...) {
		if(is.null(dim(X)))
			return(X)
		apply(X, MARGIN, FUN, ...)
	}
	
	if(showProgressBar)  
		pb <- ProgressBar(title="Getting RIs...", label="File in processing...")
  for(i in 1:length(Lib)){
  	if(showProgressBar)  
	   setProgressBar(pb, value=i/length(Lib),
				title=sprintf("Getting RIs (%d %%)",round(100*i/length(Lib))),
						label=sprintf("Metabolite %d",i))
       cor_RI[i,] <- apply2(resRI[[i]][met_cor[[i]],], 2, median, na.rm=T)
  }
  if(showProgressBar)  
		close(pb)
  return(cor_RI)
}

