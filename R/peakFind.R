
peakFind <-
function(samples, Lib, cor_RI, columns = c("SPECTRUM", "RETENTION_TIME_INDEX"),
	showProgressBar = FALSE) {

	my.files <- RIfiles(samples)
	my.names <- sampleNames(samples)
	
	resInt <- matrix(nrow=length(unlist(topMass(Lib))), ncol=length(my.files))
	colnames(resInt) <- my.names
	rownames(resInt)  <- libId(Lib, sel = FALSE)
	
	resRI <- resInt
  
	if(showProgressBar)  
		pb <- ProgressBar(title="Finding Peaks...", label="File in processing...")  
	for(i in 1:length(my.files)){
		if(showProgressBar)
			setProgressBar(pb, value=i/length(my.files),
				title=paste("Findind Peaks (", round(100*i/length(my.files)), "%)"),
				label=sprintf("Reading File %s", basename(my.files[i])))

		# use top masses
		tmpLib <- refLib(Lib, cor_RI[,i], w = 3, sel = FALSE)
		RES <- FindPeaks(my.files[i], tmpLib, columns)
		resInt[,i] <- Intensity(RES)
		resRI[,i]  <- retIndex(RES)
	}
  
	if(showProgressBar)  
	  close(pb)
	
	return(new("tsMSdata", RI = resRI, Intensity = resInt))
}
