`ImportSamples` <-
function(sampfile, CDFpath = ".", RIpath = ".", ftype=c("binary", "text"), ...) {
	ftype <- pmatch(ftype[1], c("binary", "text"))
	stopifnot(!is.na(ftype))

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

	ext     <- c("dat","txt")[ftype]
	RIfiles <- sub("cdf$", ext, paste("RI_" , Samples$CDF_FILE, sep = ""), ignore.case = T)
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

`ImportSamplesFromDir` <-
function(CDFpath=".", RIfiles=FALSE, ignore.case=TRUE, ftype=c("binary", "text")) {

    if(RIfiles == TRUE) {
        cdffiles <- dir(path=CDFpath, pattern='^RI_', ignore.case=ignore.case)
        ext <- unique(substring(tolower(cdffiles), nchar(cdffiles)-2))
        if(!all(ext %in% c("dat","txt")))
            stop("RI file extension error. It must be either 'dat' or 'txt'")

        if(length(ext) == 1) {
            cdffiles <- gsub("RI_", "", gsub(ext, "cdf", cdffiles, ignore.case=ignore.case))
            ftype <- if(ext=="dat") "binary" else "text"
        } else {
            ft <- pmatch(ftype, c("binary", "text"))[1]
            if(is.null(ft) | is.na(ft) | length(ft)==0)
                stop(paste("Binary and text RI files found. Please select one file type",
                           "with the parameter 'ftype'."))
            pattern <- c("\\.dat$","\\.txt$")[ft]
            cdffiles <- gsub("RI_", "", gsub(pattern, ".cdf", cdffiles, ignore.case=ignore.case))
            cdffiles <- grep("cdf", cdffiles, value=TRUE, ignore.case=ignore.case)
        }
    } else {
        cdffiles <- dir(path=CDFpath, pattern='\\.cdf$', ignore.case=ignore.case)
    }
    if(length(cdffiles) == 0)
        stop('Error: No CDF files were find in the directory')

    new('tsSample', CDFfiles=cdffiles,RIpath=CDFpath,CDFpath=CDFpath, ftype=ftype)
}
