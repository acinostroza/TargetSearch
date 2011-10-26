
# function to make a boxplot of the RI deviations of a given metabolite.

plotRIdev <- function(Lib, peaks, libId = 1) {

	n <- length(libId)
	j <- ceiling(sqrt(n))
	i <- ceiling(n/j)
	par(mfrow = c(i,j), mar = c(3,3,3,2), mgp = c(1.5,0.75,0))
	for(id in libId) {
		if(is.null(retIndex(peaks)[[id]]))
		    stop(sprintf("Error: libId=%d not found", id))

		x <- t(retIndex(peaks)[[id]])
		lib.name <- libName(Lib)[id]
	    mz <- selMass(Lib)[[id]]
    	if(any(mz != colnames(x)))
        	mz <- topMass(Lib)[[id]]
	    if(any(mz != colnames(x)))
    	    stop("LibraryID Error: Library object and tsMSdata object don't match.")

		if(all(is.na(x))) {
		    warning(sprintf("All values are NAs for libId=%d", id))
	    	plot.NA(main = lib.name, xlab = "mz", ylab = "RI")
		} else {
			boxplot(data.frame(x), names = mz, main = lib.name, xlab = "mz", ylab = "RI")
			abline( h = medRI(Lib)[id], col = "red")
		}
	}
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
}

# function to plot a empty box with "NA" text.
plot.NA <- function(...) {
	plot(1, type = "n", axes = FALSE, ...)
	text(1,1, "NA", cex = 2)
	box()
}

# function to plot a chromatographic peak.

plotPeak <- function(rawpeaks, time.range, masses, cdfFile = NULL, useRI = FALSE, rimTime = NULL,
	standard = NULL, massRange = c(85, 500), ...) {

	if(is.null(cdfFile) == FALSE)
		rawpeaks <- peakCDFextraction(cdfFile, massRange)

	if(length(time.range) != 2)
		stop("time.range length must be equal to 2.")

	if(useRI == TRUE) {
	    tm <- rt2ri(rawpeaks$Time,rimTime, standard)
	    xlab <- "RI"
	} else {
	    tm <- rawpeaks$Time
		xlab <- "RT"
	}

	ms  <- masses - massRange[1] + 1
	idx <- which(tm > time.range[1] & tm < time.range[2] )
	par(mar = c(5,4,5,2)+0.1, mgp = c(2,0.8,0))
	matplot(tm[idx], rawpeaks$Peaks[idx, ms ], type = 'l', xlab = xlab, ylab = "Intensity", ...)
	if(useRI == TRUE) {
        axis(3, at = tm[idx[1:(length(idx/10)/10)*10]], label = rawpeaks$Time[idx[1:(length(idx/10)/10)*10]])
	}
	legend("topright", legend = masses, col = 1:6, lty = 1:5)
}

# function to plot the median intensities across all the samples with the reference spectrum
# for a given metabolite.

plotSpectra <- function(Lib, peaks, libId = 1, type = "ht") {

	ptype <- pmatch(type, c("ht", "ss", "diff"))
	if(is.na(ptype))
		stop("Unknown type parameter ", type)

    id <- libId
	x <- t(Intensity(peaks)[[id]])

	if(all(is.na(x))) {
		plot.NA(main = libName(Lib)[id], xlab = "mz", ylab = "Intensity")
	} else {
		# remove samples with no data.
		x <- x[apply(x, 1, function(x) all(is.na(x))) == FALSE,,drop=FALSE]

		# remove masses with no data
		bar      <- apply(x, 2, function(x) all(is.na(x))) == FALSE
		x        <- x[,bar,drop=FALSE]
		# x        <- apply(x, 1,median,na.rm = T), FUN = "/")

		x.median <- apply(x, 2, median, na.rm = T)
		x.median <- 999 * x.median / max(x.median)

		mz <- topMass(Lib)[[id]][bar]
		if(length(spectra(Lib)) > 0 ) {
    		sp.mz  <- spectra(Lib)[[id]][,1]
	       	sp.int <- spectra(Lib)[[id]][,2]
     	} else {
            sp.mz <- mz
            sp.int <- rep(0, length(sp.mz))
     	}

		if(ptype == 1) {
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
		} else if (ptype == 2) {
			plot  (mz, x.median, type = 'h', col = 'blue', ylim = c(0,1000), main = libName(Lib)[id],
				xlab = "mz", ylab = "Intensity")
			points(sp.mz+0.5, sp.int, type = 'h', col = 'red')
			o1 <- order(x.median, decreasing = TRUE)[1:min(6,length(x.median))]
			text(mz[o1], x.median[o1], as.character(mz[o1]), cex = 0.7)
			legend("topright", c("reference spectrum", "median spectrum"), text.col = c("red","blue"), box.lty = 0, cex = 0.8)

		} else if (ptype == 3) {

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
}

# a wrapper function to plotSpectrum to plot all spectra
plotAllSpectra <- function(Lib, peaks, type = "ht", pdfFile, width = 8, height = 8, ...) {
	pdf(pdfFile, width, height, ...)
	on.exit(dev.off())
	for(i in 1:length(Lib)) {
		plotSpectra(Lib, peaks, i, type)
	}
}
