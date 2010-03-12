# this function saves all results in a text format
Write.Results <- function(Lib, metabProfile, prefix = NA) {
    if(missing(prefix))
        prefix <- paste("TargetSearch-", Sys.Date(), sep = "")

    file.peak.int <- paste(prefix, ".peak.intensity.txt", sep ="")
    file.peak.ri  <- paste(prefix, ".peak.RI.txt", sep = "")
    file.pro.info <- paste(prefix, ".profile.info.txt", sep ="")
    file.pro.int  <- paste(prefix, ".profile.intensities.txt", sep ="")
    file.pro.ri   <- paste(prefix, ".profile.ri.txt", sep ="")

    write.table(profileInfo(metabProfile), file = file.pro.info, sep = "\t", quote = F)
    write.table(profileInt(metabProfile), file = file.pro.int, sep = "\t", quote = F)
    write.table(profileRI(metabProfile), file = file.pro.ri,	sep = "\t", quote = F)

    metInfo   <- profileInfo(metabProfile)
    metPrf    <- as.list(metabProfile)
    libIndex  <- attr(metPrf$RI, "index")
    libId     <- rownames(metInfo)

    Name      <- as.character(metInfo[libIndex, "Name"])
    libRI     <- metInfo[libIndex, "Lib_RI"]
    mass      <- unlist(topMass(Lib)[libId])
    is_sel    <- unlist(lapply(libId, function(x) topMass(Lib)[[x]] %in% selMass(Lib)[[x]]))
    corMass   <- lapply(strsplit(metInfo$Masses, ";"), function(x) as.numeric(x))
    is_cor    <- unlist(lapply(seq(libId), function(x) topMass(Lib)[[libId[x]]] %in% corMass[[x]]))
    
    Out <- data.frame(libIndex=libIndex, Name= Name, Mass = mass, IS_SEL = is_sel, IS_COR=is_cor)
		
		write.table( data.frame(Out, metPrf$Intensity, row.names = NULL,
            check.names = FALSE), file = file.peak.int,	sep = "\t", quote = F, row.names = F)
		write.table( data.frame(Out, metPrf$RI, row.names = NULL,
            check.names = FALSE), file = file.peak.ri, sep = "\t", quote = F, row.names = F)
		
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
