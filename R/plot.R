
# function to make a boxplot of the RI deviations of a given metabolite.

plotRIdev <- function(Lib, peaks, libID = 1, ...) {
	equal <- function(x, y) isTRUE(all.equal(x, y))
    # check if `libId` was given
    if('libId' %in% names(list(...))) {
        warning('The parameter "libId" is deprecated, use "libID" instead.')
        libID <- list(...)$libId
    }
	n <- length(libID)
	j <- ceiling(sqrt(n))
	i <- ceiling(n/j)

	op <- par(no.readonly=TRUE)
	on.exit(par(op))

	par(mfrow = c(i,j), mar = c(3,3,3,2), mgp = c(1.5,0.75,0))
	for(id in libID) {
		if(is.null(retIndex(peaks)[[id]]))
		    stop(sprintf("Error: libID=%d not found", id))

		x <- t(retIndex(peaks)[[id]])
		mz <- as.numeric(colnames(x))
		lib.name <- libName(Lib)[id]
		sel_mz <- selMass(Lib)[[id]]
		top_mz <- topMass(Lib)[[id]]

		if(!equal(mz, sel_mz) && !equal(mz, top_mz))
			stop("LibraryID Error: Library object and tsMSdata object don't match.")

		if(all(is.na(x))) {
		    warning(sprintf("All values are NAs for libID=%d", id))
	    	plot.NA(main = lib.name, xlab = "mz", ylab = "RI")
		} else {
			boxplot(data.frame(x), names = mz, main = lib.name, xlab = "mz", ylab = "RI", ...)
			abline( h = medRI(Lib)[id], col = "red")
		}
	}
    invisible()
}

# a wrapper function to make a PDF of all the RIdev of all metabolites
plotAllRIdev <- function(Lib, peaks, pdfFile, width = 8, height = 8,...) {
	n <- length(Lib)
	pdf(pdfFile, width, height, ...)
	on.exit(dev.off())
	for(i in 1:ceiling(n/9)) {
		x <- 9*i - 8
		y <- min(9*i, n)
		plotRIdev(Lib, peaks, x:y)
	}
    invisible()
}

# function to plot a empty box with "NA" text.
plot.NA <- function(...) {
	plot(1, type = "n", axes = FALSE, ...)
	text(1,1, "NA", cex = 2)
	box()
}

# function to plot a chromatographic peak.

plotPeakSimple <- function(rawpeaks, time.range, masses, cdfFile = NULL, useRI = FALSE, rimTime = NULL,
	standard = NULL, massRange = NULL, ...) {

	if(is.null(cdfFile) == FALSE)
		rawpeaks <- peakCDFextraction(cdfFile, massRange)

	if(length(time.range) != 2)
		stop("time.range length must be equal to 2.")

	if(is.null(massRange)) {
		if(!is.null(rawpeaks$massRange))
			massRange <- rawpeaks$massRange
		else
			stop("mass range definition is missing.")
	}

	if(useRI == TRUE) {
	    tm <- rt2ri(rawpeaks$Time,rimTime, standard)
	    xlab <- "RI"
	} else {
	    tm <- rawpeaks$Time
		xlab <- "RT"
	}

	op <- par(no.readonly=TRUE)
	on.exit(par(op))

	ms  <- masses - massRange[1] + 1
	idx <- which(tm > time.range[1] & tm < time.range[2] )
	par(mar = c(5,4,5,2)+0.1, mgp = c(2,0.8,0))
	matplot(tm[idx], rawpeaks$Peaks[idx, ms ], type = 'l', xlab = xlab, ylab = "Intensity", ...)
	if(useRI == TRUE) {
        axis(3, at = tm[idx[1:(length(idx/10)/10)*10]], labels = rawpeaks$Time[idx[1:(length(idx/10)/10)*10]])
	}
	legend("topright", legend = masses, col = 1:6, lty = 1:5)
    invisible()
}

# a new version
plotPeak <- function(samples, Lib, metProf, rawpeaks, which.smp=1, which.met=1, massRange=NULL, corMass=FALSE)
{
	grep2 <- function(pattern, x) {
		out <- grep(tolower(pattern), tolower(x), fixed=TRUE)
		if(length(out) == 0) # trying regular expression
			out <- grep(tolower(pattern), tolower(x), ignore.case=TRUE)
		out
	}

	if(is.character(which.smp)) {
		which.smp <- grep2(which.smp[1], sampleNames(samples))
		if(length(which.smp) > 1) {
			message("Multiple samples found. Using the first match.")
			message(sprintf("[%d] %s\n", which.smp, sampleNames(samples)[which.smp]))
		} else if(length(which.smp) == 0) {
			stop("No samples matching pattern found")
		}
	}
	which.smp <- which.smp[1]

	if(is.character(which.met)) {
		which.met <- grep2(which.met[1], libName(Lib))
		if(length(which.met) > 1) {
			message("Multiple metabolites found. Using the first match.")
			message(sprintf("[%d] %s\n", which.met, libName(Lib)[which.met]))
		} else if(length(which.met) == 0) {
			stop("No metabolites found.")
		}
	}
	which.met <- which.met[1]

	cdfFile   <- CDFfiles(samples)[which.smp]
	riFile    <- RIfiles(samples)[which.smp]

	if(missing(rawpeaks) || is.null(rawpeaks))
		rawpeaks <- peakCDFextraction(cdfFile, massRange)
	if(is.null(massRange)) {
		if(!is.null(rawpeaks$massRange))
			massRange <- rawpeaks$massRange
		else
			stop("mass range definition is missing.")
	}

	# code to transform from RI to RT
	if(is.null(rawpeaks$Index)) {
		ftype <- file_type(riFile)
		if(ftype == 1) {
			if(is.integer(cols <- get.columns.name()))
				cols <- cols + 1
			ri_col <- cols[2]
			rt_col <- cols[3]
			tmp  <- read.delim(riFile, as.is = TRUE)
			ri   <- rt2ri(rawpeaks$Time, tmp[, rt_col], tmp[, ri_col])
		} else if(ftype == 0) {
			tmp <- readRIBin(riFile)
			ri  <- rt2ri(rawpeaks$Time, tmp$retTime, tmp$retIndex)
		} else {
			stop('Unable to determine file type')
		}
		rawpeaks$Index <- ri
	}

	id <- as.character(which.met)
	topMz <- topMass(Lib)[[which.met]]
	topMz <- topMz[topMz >= massRange[1] & topMz <= massRange[2]]
	corMz <- profileInfo(metProf)[id, "Masses"]
	corMz <- as.numeric(unlist(strsplit(corMz,";")))

	if(corMass) topMz <- corMz

	font <- lwd <- rep(1, length(topMz))
	lwd[topMz %in% corMz] <- 3
	font[topMz %in% corMz] <- 2

	libOriRI <- libRI(Lib)[which.met]
	libMedRI <- medRI(Lib)[which.met]
	smpMedRI <- median(retIndex(metProf)[[which.met]][, which.smp], na.rm=TRUE)
	rdev <- RIdev(Lib)[which.met,]

	riRange <- range(c(c(libOriRI,libMedRI,smpMedRI)-rdev,
		c(libOriRI,libMedRI,smpMedRI)+rdev), na.rm=TRUE)

	idx <- which( ri >= riRange[1] & ri <= riRange[2])
	mz  <- topMz - massRange[1] + 1

	intRange <- range(rawpeaks$Peaks[idx, mz ])

	main <- sprintf("Sample: %s | Metab: %s",  sampleNames(samples)[which.smp], libName(Lib)[which.met])

	plot(NA, type="n", xlab = 'Retention Index (RI)', ylab = "Intensity", xlim=range(ri[idx]), ylim=intRange)
	mtext(main, 3, font=2, line=2)
	axis(3, at = ri[idx[1:(length(idx/10)/10)*10]], labels = rawpeaks$Time[idx[1:(length(idx/10)/10)*10]])

	rect(libMedRI-rdev[2], -intRange[2]*2, libMedRI+rdev[2], intRange[2]*2, col=gray(0.95), lty=0)
	rect(smpMedRI-rdev[3], -intRange[2]*2, smpMedRI+rdev[3], intRange[2]*2, col=gray(0.9), lty=0)

	matlines(ri[idx], rawpeaks$Peaks[idx, mz ], lwd=lwd)
	legend("topright", legend = topMz, col = 1:6, lty = 1:5, lwd=lwd, text.col=font)
	box()

	invisible(rawpeaks)
}

# function to plot the median intensities across all the samples with the reference spectrum
# for a given metabolite.

plotSpectra <- function(Lib, peaks, libID = 1, type = c("ht", "ss", "diff")) {

	type <- match.arg(type, c("ht", "ss", "diff"))

    id <- libID
	x <- t(Intensity(peaks)[[id]])

	if(all(is.na(x))) {
		plot.NA(main = libName(Lib)[id], xlab = "mz", ylab = "Intensity")
	} else {
		# remove samples with no data.
		x <- x[apply(x, 1, function(x) all(is.na(x))) == FALSE,,drop=FALSE]

		# remove masses with no data
		bar      <- apply(x, 2, function(x) all(is.na(x))) == FALSE
		x        <- x[,bar,drop=FALSE]

		x.median <- apply(x, 2, median, na.rm = TRUE)
		x.median <- 999 * x.median / max(x.median)

		mz <- topMass(Lib)[[id]][bar]
		if(length(spectra(Lib)) > 0 ) {
    		sp.mz  <- spectra(Lib)[[id]][,1]
	       	sp.int <- spectra(Lib)[[id]][,2]
     	} else {
            sp.mz <- mz
            sp.int <- rep(0, length(sp.mz))
     	}

		if(type == "ht") {
			plot  (mz, x.median, type = 'h', col = 'blue', ylim = c(-1000,1000), main = libName(Lib)[id],
				xlab = "mz", ylab = "Intensity", yaxt = "n")
			points(sp.mz, - sp.int, type = 'h', col = 'red')
			axis(2, at = c(-1000,-500,0,500,1000), labels =  c(1000,500,0,500,1000))

			o1 <- order(x.median, decreasing = TRUE)[1:min(6,length(x.median))]
			o2 <- order(sp.int, decreasing = TRUE)[1:min(6,length(sp.int))]
			text(mz[o1], x.median[o1], as.character(mz[o1]), cex = 0.7)
			text(sp.mz[o2], -sp.int[o2], as.character(sp.mz[o2]), cex = 0.7)
			legend("bottomright", "reference spectrum", box.lty = 0, cex = 0.8)
			legend("topright", "median spectrum", box.lty = 0, cex = 0.8)
		} else if (type == "ss") {
			plot  (mz, x.median, type = 'h', col = 'blue', ylim = c(0,1000), main = libName(Lib)[id],
				xlab = "mz", ylab = "Intensity")
			points(sp.mz+0.5, sp.int, type = 'h', col = 'red')
			o1 <- order(x.median, decreasing = TRUE)[1:min(6,length(x.median))]
			text(mz[o1], x.median[o1], as.character(mz[o1]), cex = 0.7)
			legend("topright", c("reference spectrum", "median spectrum"), text.col = c("red","blue"), box.lty = 0, cex = 0.8)

		} else if (type == "diff") {

			plot(mz, x.median, type = 'n', ylim = c(-1000,1000), main = libName(Lib)[id],
				xlab = "mz", ylab = "Intensity", yaxt = "n")
			axis(2, at = c(-1000,-500,0,500,1000), labels =  c(1000,500,0,500,1000))

			foo <- mz %in% sp.mz
			points(mz[!foo], x.median[!foo], col = "blue", type = 'h')

			bar <- sp.mz %in% mz
			points(sp.mz[!bar], -sp.int[!bar], col = "red", type = 'h')

			mz <- mz[foo]
			x.median <- x.median[foo]
			mz.diff.int <- apply(cbind(mz, x.median), 1, function(x) x[2] - sp.int[sp.mz == x[1]])
			points(mz, mz.diff.int, col = "darkgreen", type = 'h')
			legend("bottomright", "reference spectrum", box.lty = 0, cex = 0.8)
			legend("topright", "median spectrum", box.lty = 0, cex = 0.8)
		}
	}
    invisible()
}

# a wrapper function to plotSpectrum to plot all spectra
plotAllSpectra <- function(Lib, peaks, type = "ht", pdfFile, width = 8, height = 8, ...) {
	pdf(pdfFile, width, height, ...)
	on.exit(dev.off())
	for(i in 1:length(Lib)) {
		plotSpectra(Lib, peaks, i, type)
	}
    invisible()
}

# vim: set ts=4 sw=4 noet:
