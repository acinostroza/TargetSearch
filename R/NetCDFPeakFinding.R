# function to find peaks in chromatograms
# arguments:
#    - cdfFile : The file
#    - massRange: m/z range.
#    - Window: Window to be used by method
#    - IntThreshold: Intensity threshold
#    - method: peak picking method. PPC / smoothing

NetCDFPeakFinding <- function(cdfFile, massRange = NULL, Window = 15, IntThreshold = 10,
	pp.method = "ppc", baseline = FALSE, baseline.opts = NULL) {

	method <- pmatch(pp.method,c("smoothing", "ppc"))

	if(is.na(method))
		stop("Invalid peak picking method.")

	ncData <- .open.ncdf(cdfFile)

    if(any(ncData$scanindex < 0)) {
        message('Warning:')
        message('  The following file seems to be corrupted. TargetSearch will attempt to process it anyway...')
        message(paste('  ->', cdfFile))
        # removing negative values
        tmp <- ncData$scanindex >= 0
        ncData$scanindex <- ncData$scanindex[tmp]
        ncData$point_count <- ncData$point_count[tmp]
        ncData$rt        <- ncData$rt[tmp]
    }

	ncData <- .check.mz.precision(ncData)
	if(is.null(ncData)) {
		stop(paste("Error processing file '", cdfFile, "'. It seems to contains",
			"high mass accuracy data.", sep=" "))
	}

	if(baseline)
		ncData <- baseline(ncData, baseline.opts)

	if(is.null(massRange))
		massRange <- range(ncData$mz)

	if(method == 1)
		peaks <- .Call("peak_finding", ncData$mz, ncData$intensity, ncData$point_count,
			ncData$scanindex, Window, massRange, IntThreshold, PACKAGE="TargetSearch")

	if(method == 2)
		peaks <- .Call("ppc",  ncData$mz, ncData$intensity, ncData$point_count,
			ncData$scanindex, Window, massRange, IntThreshold, PACKAGE="TargetSearch")
	colnames(peaks) <- as.character(massRange[1]:massRange[2])
	return( list(Time = ncData$rt, Peaks = peaks, massRange=massRange) )
}

# imported functions
# open.ncdf,  get.var.ncdf, close.ncdf

.open.ncdf <- function(cdfFile) {
	nc <- open.ncdf(cdfFile)
	ncData <- list()
	ncData$point_count <- get.var.ncdf(nc, "point_count")
	ncData$scanindex   <- get.var.ncdf(nc, "scan_index")
	ncData$intensity   <- get.var.ncdf(nc, "intensity_values")
	ncData$mz          <- get.var.ncdf(nc, "mass_values")
	ncData$rt          <- get.var.ncdf(nc, "scan_acquisition_time")
	close.ncdf(nc)
	return(ncData)
}

# vim: set ts=4 sw=4:
