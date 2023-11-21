`ImportFameSettings` <-
function(tmp.file=NA, mass = NA, ...) {
  ## next line specifies the correct RI-value's of your marker metabolites
  rim.perfect <- c(262320, 323120, 381020, 487220, 582620, 668720, 747420,
                   819620, 886620, 948820, 1006900, 1061700, 1113100)

  ## here you set up the search frames for these markers (in seconds)
  ## if you specify a path for 'tmp.file' you can load settings from a tab-delimited file
  if (is.na(tmp.file)) {
    rim.limits <- rbind(c(230,280), # lower/upper/RI_correct_value limit of first Marker
                        c(290,340),       # lower/upper limit of second Marker
                        c(350,400),       # ...
                        c(450,500),
                        c(540,590),
                        c(630,665),
                        c(705,745),
                        c(775,820),
                        c(845,885),
                        c(905,940),
                        c(965,1015),
                        c(1020,1060),
                        c(1070,1110))
  } else {
    rim.limits <- as.matrix(read.delim(tmp.file, ...))
    if(ncol(rim.limits) != 3 & ncol(rim.limits) != 4)
	  	stop("Error reading FAME file. The file doesn't have 3 or 4 columns",
			" (LowerLimits, UpperLimit, RIperfect, [mz marker] )")

		if(any(rim.limits[,1] > rim.limits[,2]))
			stop("Error: LowerLimits are greater than UpperLimits. Please check your file (rows ", 
				paste(which(rim.limits[,1] > rim.limits[,2]), collapse = " "), ")")

		if(ncol(rim.limits) == 4)
			mass <- as.numeric(rim.limits[,4])

		rim.perfect <- rim.limits[,3]
		rim.limits  <- rim.limits[,1:2]
  }
  if(any(is.na(mass)))
  	mass <- 87

  colnames(rim.limits) <- c("LowerLimit", "UpperLimit")
	if(is.null(rownames(rim.limits)) | all(rownames(rim.limits) == 1:nrow(rim.limits) ))
		rownames(rim.limits) <- paste("RI.Marker", 1:nrow(rim.limits))
  
  new("tsRim", limits = rim.limits, standard = rim.perfect, mass = mass)
}
