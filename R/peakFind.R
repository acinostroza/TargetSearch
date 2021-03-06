
peakFind <-
function(samples, Lib, cor_RI, columns = NULL, showProgressBar = FALSE)
{
	my.files <- RIfiles(samples)
	my.names <- sampleNames(samples)

	refLib <- lapply(seq(ncol(cor_RI)), function(x) refLib(Lib, cor_RI[,x], w=3, sel=FALSE))
	RES    <- FindPeaks(my.files, refLib, columns, showProgressBar)
	resInt <- Intensity(RES)
	resRI  <- retIndex(RES)
    resRT  <- retTime(RES)

    for(i in 1:length(resInt))
        colnames(resRT[[i]]) <- colnames(resRI[[i]]) <- colnames(resInt[[i]]) <- my.names

    names(resRI) <- names(resInt) <- names(resRT) <- rownames(cor_RI)

    return(new("tsMSdata", RI = resRI, Intensity = resInt, RT = resRT))
}
