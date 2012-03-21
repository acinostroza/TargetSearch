`writeRIFile` <-
function(outFile, Peaks, riInde, massRange, ftype=c("binary", "text")) {
	ftype <- pmatch(ftype[1], c("binary", "text"))

	if(ftype == 2) {
		Header <- "RETENTION_TIME\tSPECTRUM\tRETENTION_TIME_INDEX"
		res <- .C("writePeaksTXT", as.character(outFile), as.double(Peaks$Time),
			as.double(riInde), as.integer(Peaks$Peaks), as.integer(massRange),
			as.integer(length(Peaks$Time)), as.character(Header), PACKAGE="TargetSearch")
	} else if(ftype == 1) {
		swap <- pmatch(.Platform$endian, c("little", "big")) - 1
		res <- .C("writePeaksDAT", as.character(outFile), as.double(Peaks$Time),
			as.double(riInde), as.integer(Peaks$Peaks), as.integer(massRange),
			as.integer(length(Peaks$Time)), as.integer(swap), PACKAGE="TargetSearch")
	} else {
		stop("Error: invalid parameter ftype")
	}
	invisible(1)
}
