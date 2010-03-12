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

