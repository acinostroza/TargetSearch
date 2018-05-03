# a wrapper for baseline Correction Function

baseline <- function(ncData, baseline.opts = NULL) {

    foo <- function(x, y, si, pc) y[si[x] + 1:pc[x]]

    if(all(ncData$point_count == ncData$point_count[1]))
    {
        int <- sapply(seq_along(ncData$scanindex), foo, ncData$intensity, ncData$scanindex, ncData$point_count)
        if(!is.matrix(int))
            stop("Unexpected Error: Object 'int' is not a matrix")
    }
    else
    {
        massRange <- range(ncData$mz)
        int <- t(.Call("ncdfToMatrix", ncData, massRange, PACKAGE="TargetSearch"))
    }

    int <- do.call(baselineCorrection, append(list(int = int), baseline.opts))
    ## update ncData components...
    int.id  <- lapply(1:length(ncData$scanindex), function(x) which( int[,x] > 0 ))
    ncData$intensity <- as.vector( int[int > 0] )
    ncData$mz <- range(ncData$mz)
    ncData$mz <- ncData$mz[1]:ncData$mz[2]
    ncData$mz <- unlist(sapply(int.id, function(x) ncData$mz[x]))
    ncData$point_count  <- sapply(int.id, length)
    ncData$scanindex  <- c(0, cumsum(ncData$point_count))[1:length(ncData$point_count)]
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
    segments = 100, signalWindow = 10, method = "linear") {

        n <- ncol(int)
        bfrac <- round(bfraction * segments)
        met <- pmatch(method, c("linear", "spline"))[1]
        if(is.na(met))
            stop("Invalid method '", method, "' for baseline correction")

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
                if(met == 2) {
                # Fit a cubic smoothing spline of the baseline using the points
                # considerated noise.
                    sp <- smooth.spline((1:n)[!sm], x[!sm])
                    spx <- predict(sp, 1:n)
                } else if (met == 1)
                    spx <- approx( (1:n)[!sm], x[!sm], 1:n )
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

# vim: set ts=4 sw=4 et:
