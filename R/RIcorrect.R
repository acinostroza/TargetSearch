RIcorrect <- function(samples, rimLimits = NULL, massRange=NULL, Window, IntThreshold,
	pp.method="ppc", showProgressBar=FALSE, baseline=FALSE, writeCDF4path=TRUE, ...)
{
	manyFiles <- CDFfiles(samples)
	outFile   <- RIfiles(samples)
	Names     <- sampleNames(samples)
	ftype     <- fileFormat(samples)
	if(is.na(ftype)) ftype <- "binary"

	if(is.null(rimLimits) == FALSE) {
		standard  <- rimStandard(rimLimits)
		mass      <- rimMass(rimLimits)
		rLimits   <- rimLimits(rimLimits)
		RIcheck <- matrix(nrow=dim(rLimits)[1], ncol=length(Names))
	}

	# check Files
	if(!all(file.exists(manyFiles))) {
		stop("These files don't exist: ", paste(manyFiles[!file.exists(manyFiles)], collapse = " "))
	}

	# get cdf4Files depending on the value of writeCDF4path.
	cdf4Files <- .mk_cdf_files(manyFiles, writeCDF4path)
	writeCDF4 <- !is.null(cdf4Files)

	if(showProgressBar)
		pb <- ProgressBar(title="Extracting peaks...", label="File in processing...")

	for(i in 1:length(manyFiles)) {
		ncdfInfo <- .ncdf_info(manyFiles[i])
		ncdf <- peakCDFextraction(manyFiles[i])

		# baseline correction
		if(baseline) {
			if(ncdfInfo$baseline_corrected == 1)
				warning("File '", manyFiles[i], "' already baseline corrected. Skipping")
			else
				ncdf <- baseline(ncdf, ...)
		}

		# writing NCDF-4
		if(writeCDF4) {
			if(ncdfInfo$creator == 'TargetSearch') {
				if(ncdfInfo$baseline_corrected == 0 & baseline)
					ncdf4_write(cdf4Files[i], ncdf)
			} else {
				ncdf4_write(cdf4Files[i], ncdf)
			}
		}

		Peaks  <- NetCDFPeakFinding(ncdf, massRange, Window, IntThreshold, pp.method = pp.method)

		# mass range has no effect, so we ignore
		massRange <- Peaks$massRange

		if(is.null(rimLimits) == FALSE) {
			# check that the mass of rimLimits is within the mass range
			if(any(mass < massRange[1] | mass > massRange[2]))
				stop(sprintf(
					paste("m/z of markers out of range:",
						" => file: '%s' | m/z range: %d, %d | m/z out of range: %s",
						sep="\n"),
					manyFiles[i], massRange[1], massRange[2],
					paste(mass[mass < massRange[1] | mass > massRange[2]], collapse=", ")))

			fameTimes <- findRetentionTime(Peaks$Time, Peaks$Peaks[, mass - massRange[1] + 1], rLimits)
			RIcheck[,i] <- fameTimes
			riInde <- rt2ri(Peaks$Time, fameTimes, standard)
			writeRIFile(outFile[i], Peaks, riInde, massRange, ftype)

			# update RI on CDF file
			if(writeCDF4)
				ncdf4_update_ri(cdf4Files[i], fameTimes, standard)

		} else {
			writeRIFile(outFile[i], Peaks, Peaks$Time, massRange, ftype)
		}

		if(showProgressBar)
			setProgressBar(pb, value=i/length(manyFiles),
				title=paste("Extracting peaks (", round(100*i/length(manyFiles)), "%)"),
				label=basename(manyFiles[i]))
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

#' make NCDF4 from CDF files.
#'
#' Generates a list of NCDF4 files from a list of CDF files (formats 3 or 4).
#' The output is controlled by the parameter path, for which there are several
#' options:
#'  - If path is FALSE or NULL or NA, then return NULL
#'  - If path is TRUE, return a list of NCDF4 files. If the files are already
#'    NCDF4 files, then output them.
#-  - If path is a string, then update the files' path to this value. Recycle
#'    as necessary.
#' @param files A list of NCDF3 or NCDF4 files
#' @param path A new path or FALSE.
#' @exts A list of valid extension which will be used to remove from the input files.
#' get a list nc4 files. If the files are already nc4 files,
.mk_cdf_files <- function(files, path, exts=c('nc4','cdf'))
{
	if(is_nullOrNA(path))
		return(NULL)
	if(is.logical(path)) {
		if(all(path))
			return(sprintf("%s.nc4", .trim_file_ext(files, exts)))
		return(NULL)
	}
	if(is.character(path)) {
		files <- .setpath(files, path)
		return(sprintf("%s.nc4", .trim_file_ext(files, exts)))
	}
	stop("Invalid parameter `path`. logical or character argument expected")
}

# vim: set ts=4 sw=4:
