`writeRIFile` <-
function(outFile, Peaks, riInde, massRange) {
	Header <- c("RETENTION_TIME","SPECTRUM","RETENTION_TIME_INDEX")
	# Get rid of scans width no peaks
	time_idx <- which(apply(Peaks$Peaks, 1, max) > 0)

	format_spectra <- function (x, minmass = massRange[1]) {
		tmp <- which( x > 0 )
		minmass <- minmass - 1
		return(paste(paste(tmp+minmass, x[tmp], sep = ":"), collapse = " "))
	}
	
	spectra <- apply(Peaks$Peaks[time_idx,], 1, format_spectra)
	write.table(data.frame(Peaks$Time[time_idx],spectra, riInde[time_idx]),
	 file = outFile, col.names = Header,row.names = F, sep="\t", quote=FALSE)
}

