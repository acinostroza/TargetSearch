# this function saves all results in a text format
Write.Results <- function(Lib = NA, peakData = NA, finalProfile = NA, prefix = NA) {
	if(missing(prefix))
   	prefix <- paste("TargetSearch-", Sys.Date(), sep = "")
	
	msg <- NA
	
	if(!missing(Lib) & !missing(peakData)) {
		# save peakData
		file.peak.int <- paste(prefix, ".peak.intensity.txt", sep ="")
		file.peak.ri  <- paste(prefix, ".peak.RI.txt", sep = "")

    med_RI    <- apply(retIndex(peakData), 1, median, na.rm = T)
    libId     <- libId(Lib, sel = FALSE)
    Name      <- libName(Lib)[libId]
    libRI     <- libRI(Lib)[libId]
    mass      <- unlist(topMass(Lib))
    is_sel    <- F
    
    Lib2 <- data.frame(Lib_RI = libRI, Mass = mass, IS_SEL = is_sel)
		
		write.table( data.frame(cbind(Name, Lib2, med_RI, Intensity(peakData))), file = file.peak.int,
			sep = "\t", quote = F, row.names = F)
		write.table( data.frame(cbind(Name, Lib2, med_RI, retIndex(peakData))), file = file.peak.ri,
			sep = "\t", quote = F, row.names = F)
		
		msg <- 1
		message(
			"The following files have been saved:\n",
		  " - ",file.peak.int, " RAW intensities of the found masses.\n",
			" - ",file.peak.ri,  " RIs of the RAW intensities found.")

		if(length(retTime(peakData)) > 0) {
			file.peak.rt  <- paste(prefix, ".peak.RT.txt", sep = "")
			write.table( data.frame(cbind(Name, Lib2, med_RI, retTime(peakData))), file = file.peak.rt,
				sep = "\t", quote = F, row.names = F)
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
