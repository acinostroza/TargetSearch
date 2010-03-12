RIcorrect <- function(samples, rimLimits = NULL, massRange, Window, IntThreshold,
	pp.method = "smoothing", showProgressBar = FALSE,
	baseline = FALSE, baseline.opts = NULL ) {
	
	manyFiles <- CDFfiles(samples)
	outFile   <- RIfiles(samples)
	Names     <- sampleNames(samples)

	if(is.null(rimLimits) == FALSE) {
		standard  <- rimStandard(rimLimits)
		mass      <- rimMass(rimLimits)
		rLimits   <- rimLimits(rimLimits)
		if(any(mass < massRange[1] | mass > massRange[2]))
			stop("'mass' marker out of Range")
		RIcheck <- matrix(nrow=dim(rLimits)[1], ncol=length(Names))
	}

	# check Files
	if(!all(file.exists(manyFiles))) {
	 	stop("These files don't exist: ", paste(manyFiles[!file.exists(manyFiles)], collapse = " "))
	}

	if(showProgressBar)
		pb <- ProgressBar(title="Extracting peaks...", label="File in processing...")
	for(i in 1:length(manyFiles)) {
		if(showProgressBar)
			setProgressBar(pb, value=i/length(manyFiles),
				title=paste("Extracting peaks (", round(100*i/length(manyFiles)), "%)"),
				label=basename(manyFiles[i]))
			
		Peaks  <- NetCDFPeakFinding(manyFiles[i], massRange, Window, IntThreshold, pp.method = pp.method,
				baseline = baseline, baseline.opts = baseline.opts)
		if(is.null(rimLimits) == FALSE) {
			fameTimes <- findRetentionTime(Peaks$Time, Peaks$Peaks[, mass - massRange[1] + 1], rLimits)
			RIcheck[,i] <- fameTimes
			riInde <- rt2ri(Peaks$Time, fameTimes, standard)
			writeRIFile(outFile[i], Peaks, riInde, massRange)
		} else {
		 	writeRIFile(outFile[i], Peaks, Peaks$Time, massRange)
		}
	}
	if(showProgressBar)	
		close(pb)
	
	if(is.null(rimLimits) == FALSE) {
		colnames(RIcheck) <- Names
		rownames(RIcheck) <- rownames(rLimits)
		return(RIcheck)
	} else {
		return(NULL)
	}
}
