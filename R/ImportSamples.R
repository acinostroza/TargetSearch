`ImportSamples` <-
function(sampfile, CDFpath = ".", RIpath = ".", ...) {
	
	Samples <- read.delim(sampfile, ...)
	if(is.null(Samples$CDF_FILE)) stop("Column 'CDF_FILE' is missing!!")
	Samples$CDF_FILE <- as.character(Samples$CDF_FILE)	
	if(is.null(Samples$MEASUREMENT_DAY)) {
		Samples$MEASUREMENT_DAY <- "1"
		warning("Column 'MEASUREMENT_DAY' not found. Using default setting.")
	}
	cdf.idx <- grep("\\.cdf$", Samples$CDF_FILE, ignore.case = T)
	if(length(cdf.idx) != nrow(Samples)) {
		warning("Some CDF file names don't have a CDF extension. Those rows will be removed.")
		Samples <- Samples[cdf.idx,]
	}
	
	RIfiles <- sub("cdf$", "txt", paste("RI_" , Samples$CDF_FILE, sep = ""), ignore.case = T)
	Names   <- gsub(".cdf", "", Samples$CDF_FILE)
	new("tsSample", Names = Names, CDFfiles = Samples$CDF_FILE, days = as.character(Samples$MEASUREMENT_DAY),
		RIfiles = RIfiles, CDFpath = CDFpath, RIpath = RIpath, data = Samples)
}

# function to extract the /measurement day/ of a set of names.

getDays <- function(x) {
    y <- regexpr('[0-9]+', x)

     # check if one element doesn't have digits. If so, just return ones.
    if(any(y== -1)) {
        return(rep('1', length(x)))
    }

    d <- substring(x, y, y+attr(y, "match.length")-1)
    # count number of elements per day
    n <- sapply(unique(d), function(z) sum(z == d))
    if(any(n == 1)) {
        # warning TODO
        return(rep('1', length(x)))
    }
    return(d)
}

ImportSamplesFromDir <- function(CDFpath=".", RIfiles = FALSE, ignore.case = TRUE) {
    if(RIfiles == TRUE) {
        cdffiles <- dir(path=CDFpath, pattern='^RI_', ignore.case=ignore.case)
        cdffiles <- gsub("RI_", "", gsub("txt", "cdf", cdffiles, ignore.case=ignore.case))
    } else {
        cdffiles <- dir(path=CDFpath, pattern='\\.cdf$', ignore.case=ignore.case)
    }
    if(length(cdffiles) == 0)
        stop('Error: No CDF files were find in the directory')

    days <- TargetSearch:::getDays(cdffiles)
    new('tsSample', CDFfiles=cdffiles,RIpath=CDFpath,CDFpath=CDFpath)
}
