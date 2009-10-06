# this function saves all results in a text format
Write.Results <- function(Lib = NA, peakData = NA, finalProfile = NA, prefix = NA) {
	if(missing(prefix))
   	prefix <- paste("TargetSearch-", Sys.Date(), sep = "")
	
	msg <- NA
	
	if(!missing(Lib) & !missing(peakData)) {
		# save peakData
		file.peak.int <- paste(prefix, ".peak.intensity.txt", sep ="")
		file.peak.ri  <- paste(prefix, ".peak.RI.txt", sep = "")

    med_RI    <- format(apply(retIndex(peakData), 1, median, na.rm = T), digits = 3)
    libId     <- libId(Lib, sel = FALSE)
    Name      <- libName(Lib)[libId]
    libRI     <- libRI(Lib)[libId]
    mass      <- unlist(topMass(Lib))
    is_sel    <- unlist(lapply(seq(Lib), function(x) topMass(Lib)[[x]] %in% selMass(Lib)[[x]]))
    
    Lib2 <- data.frame(Lib_RI = libRI, Mass = mass, IS_SEL = is_sel)
		
		write.table( data.frame(libId, Name, Lib2, med_RI, Intensity(peakData), row.names = NULL,
            check.names = FALSE), file = file.peak.int,	sep = "\t", quote = F, row.names = F)
		write.table( data.frame(libId, Name, Lib2, med_RI, retIndex(peakData), row.names = NULL,
            check.names = FALSE), file = file.peak.ri, sep = "\t", quote = F, row.names = F)
		
		msg <- 1
		message(
			"The following files have been saved:\n",
		  " - ",file.peak.int, " RAW intensities of the found masses.\n",
			" - ",file.peak.ri,  " RIs of the RAW intensities found.")

		if(length(retTime(peakData)) > 0) {
			file.peak.rt  <- paste(prefix, ".peak.RT.txt", sep = "")
			write.table( data.frame(cbind(Name, Lib2, med_RI, retTime(peakData)), check.names = FALSE),
                file = file.peak.rt, sep = "\t", quote = F, row.names = F)
			message(" - ", file.peak.rt, " RTs of th RAW intensties.\n")
  	}
	}
	
	if(!missing(finalProfile)) {
		# save finalProfile		
		file.pro.info <- paste(prefix, ".profile.info.txt", sep ="")	
		file.pro.int  <- paste(prefix, ".profile.intensities.txt", sep ="")
		file.pro.ri   <- paste(prefix, ".profile.ri.txt", sep ="")
	
		write.table(profileInfo(finalProfile), file = file.pro.info, sep = "\t", quote = F)
		write.table(Intensity(finalProfile), file = file.pro.int, sep = "\t", quote = F)
		write.table(retIndex(finalProfile), file = file.pro.ri,	sep = "\t", quote = F)
		
	  if(is.na(msg)) message("The following files have been saved:\n")
	  message(" - ",file.pro.info, " Information about metabolites.\n",
		" - ",file.pro.int,  " Normalised intensities of the metabolites.\n",
		" - ",file.pro.ri,   " RIs of the metabolites.")
		
		if(length(retTime(finalProfile)) > 0) {
			# saving Retention Times
			file.pro.rt   <- paste(prefix, ".profile.rt.txt", sep ="")
			write.table(retTime(finalProfile), file = file.pro.rt,	sep = "\t", quote = F)
			message(" - ",file.pro.rt,   " RTs of the metabolites.\n")
		}
	}
}

# function to create a MSP file that can be viewed with NIST

writeMSP <- function(lib, peaks, file, append=FALSE) {

    if(is.logical(append) == FALSE) {
        stop("append must be logical")
    }
    file <- file(file, ifelse(append[1], "a", "w"))
    tmp <- t(Intensity(peaks))
    
    for(i in 1:length(lib)) {

        x <- tmp[, colnames(tmp) == i]
        if(all(is.na(x))) {
            next
        }
   		# remove samples with no data.
		x <- x[apply(x, 1, function(x) all(is.na(x))) == FALSE,]

   		# remove masses with no data
		bar      <- apply(x, 2, function(x) all(is.na(x))) == FALSE
		x        <- x[,bar]

   		x.median <- apply(x, 2, median, na.rm = T)
		x.median <- round(999 * x.median / max(x.median))
		mz <- topMass(lib)[[i]][bar]

		cat(sprintf("Name: %s", libName(lib)[i]), file = file, sep = "\n")
		cat(sprintf("Synon: RI: %.1f", medRI(lib)[i]), file = file, sep = "\n")
		cat(sprintf("Synon: MST SEL MASS: %s", substring(paste(selMass(lib)[[i]], collapse = ";"), 1, 230)), file = file, sep = "\n")
		cat(sprintf("Num Peaks: %d", length(x.median)), file = file, sep = "\n")
#		cat(paste(mz, " ", x.median, ";" , sep = ""), file = file, sep = "\n")
        foo <- paste(mz, " ", x.median, ";", sep = "")
        foo <- split(foo, rep(1:ceiling(length(foo)/6), each = 6)[1:length(foo)])
        cat (sapply(foo, function(x) paste(x, collapse = " ")), file = file, sep = "\n")
		cat("", file = file, sep = "\n")
	}
	close(file)
}
