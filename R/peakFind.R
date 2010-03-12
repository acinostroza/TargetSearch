
peakFind <-
function(samples, Lib, cor_RI, columns = c("SPECTRUM", "RETENTION_TIME_INDEX"),
	showProgressBar = FALSE) {

	my.files <- RIfiles(samples)
	my.names <- sampleNames(samples)
	
	refLib <- lapply(seq(ncol(cor_RI)), function(x) refLib(Lib, cor_RI[,x], w=3, sel=FALSE))
	RES    <- FindPeaks(my.files, refLib, columns, showProgressBar)
	resInt <- Intensity(RES)
	resRI  <- retIndex(RES)

    for(i in 1:length(resInt)) colnames(resRI[[i]]) <- colnames(resInt[[i]]) <- my.names
    names(resRI) <- names(resInt) <- rownames(cor_RI)

	return(new("tsMSdata", RI = resRI, Intensity = resInt))
}
