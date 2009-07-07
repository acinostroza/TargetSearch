# a wrapper for baseline Correction Function

baseline <- function(ncData, baseline.opts = NULL) {

	foo <- function(x, y, si, pc) y[si[x] + 1:pc[x]]
	if( all(ncData$point_count == ncData$point_count[1]) ) {
    	int <- sapply(1:length(ncData$scanindex), foo, ncData$intensity, ncData$scanindex, ncData$point_count)
    	stopifnot(is.matrix(int))
	} else {
		warning("Baseline Correction: It seems that the data is already baseline corrected.")
		int <- t(.Call("peakExtraction", ncData$mz, ncData$intensity,
            ncData$point_count, ncData$scanindex, range(ncData$mz),
            PACKAGE = "TargetSearch"))
	}

	int <- do.call(baselineCorrection, append(list(int = int), baseline.opts))
	ncData$intensity <- as.vector(int)
	ncData
}

# baseline correction
# args:
#   int : An intensity matrix where columns represent scan times and rows
#         differents mass values.
#   threshold: a threshold for the baseline. A value of one returns a baseline
#              just above the noise, 0.5 in the middle of the noise and 0 below 
#              the noise.
#   alpha: alpha parameter for the high pass filter
#   segments: the number of segments in which the time range is cut
#   bfraction: a 0-1 value representing a percentage of the segments
#         considered to be noise.
#   signalwindow:

baselineCorrection <- function(int, threshold = 0.5, alpha = 0.95, bfraction = 0.2,
    segments = 100, signalWindow = 10) {

        n <- ncol(int)
        bfrac <- round(bfraction * segments)

        baseline <- apply(int, 1, function(x) {
                # apply high pass filter
                xf <- hpf(x, alpha)
                # divide xf in segments
                np <- ceiling( n / segments )
                segment.idx  <- rep(1:segments, each = np)[1:n]
                # compute standard deviation of every segment
                segment.sd <- sapply( split(xf, segment.idx), sd)
                # get standard deviation of the 'bfraction' segments with
                # the lowest standard deviation. (An estimation of
                # the sd of the noise)
                o  <- order(segment.sd)[1:bfrac]
                stdn <- sd(xf[segment.idx %in% o])
                # points with absolute values higher than 2*stdn are considered as signal
                # and the center of a signal window of width signalwindow
                tmp  <- which(abs( xf ) > 2*stdn)
                # windowing step: apply signalwindow to the signal points obtained before
                sm   <- windowingStep(tmp, n, signalWindow)
                sm[1] <- sm[n] <- FALSE
                # Fit a cubic smoothing spline of the baseline using the points 
                # considerated noise.
                sp <- smooth.spline((1:n)[!sm], x[!sm])
                spx <- predict(sp, 1:n)
                # returns the smoothed points plus offset
                spx$y + 2 * (threshold - 0.5) * 2 * stdn
        })
        int <- int - t(baseline)
	int[int < 0] <- 0
	int
}

# high pass filter
# y[i] <- alpha * (y[i-1] + x[i] - x[i-1])

hpf <- function(x, alpha) {
        n <- length(x)
        .C("hpf", x = as.numeric(x), y = double(n), n = as.integer(n),
		 a = as.numeric(alpha), PACKAGE="TargetSearch")$y
}

windowingStep <- function(x, n, wm) {
        x.len <- length(x)
        out <- .C("windowing", s = integer(n), x = as.integer(x), wm = as.integer(wm),
                n = as.integer(n), nidx = as.integer(x.len),  PACKAGE="TargetSearch" )
        as.logical(out$s)
}
