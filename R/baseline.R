# a wrapper for baseline Correction Function
#
baseline <- function(ncdf, bsline_method = c('classic', 'quantiles', 'none'), ...)
{
    # assume `ncdf` is a file name
    ncdf <- .peakExtractWrap(ncdf)

    method <- match.arg(bsline_method)

    if(method == 'none')
        return(ncdf)
    if(method == 'classic')
        return(baselineCorrection(ncdf, ...))
    if(method == 'quantiles')
        return(baselineCorrectionQuant(ncdf, ...))
    return(NULL)
}

# baseline correction
# args:
#   int : a named two component list: `Time` and `Peaks` returned by a call to
#         peakCDFextraction.
#   threshold: a threshold for the baseline. A value of one returns a baseline
#              just above the noise, 0.5 in the middle of the noise and 0 below
#              the noise.
#   alpha: alpha parameter for the high pass filter
#   segments: the number of segments in which the time range is cut
#   bfraction: a 0-1 value representing a percentage of the segments
#         considered to be noise.
#   signalwindow:

baselineCorrection <- function(peaks, threshold = 0.5, alpha = 0.95, bfraction = 0.2,
    segments = 100, signalWindow = 10, method = "linear") {

    # core function for baseline correction. It expects a single vector
    # of intensity values
    base_fun <- function(x, threshold, alpha, bfrac, segments, signalWindow, method)
    {
        # apply high pass filter
        xf <- hpf(x, alpha)
        n  <- length(x)
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
        if(method == 'spline') {
            # Fit a cubic smoothing spline of the baseline using the points
            # considerated noise.
            sp <- smooth.spline((1:n)[!sm], x[!sm])
            spx <- predict(sp, 1:n)
        } else {
            spx <- approx( (1:n)[!sm], x[!sm], 1:n )
        }
        # returns the smoothed points plus offset
        spx$y + 2 * (threshold - 0.5) * 2 * stdn
    }

    int <- NULL
    method <- match.arg(method)
    bfrac <- round(bfraction * segments)

    if(is.list(peaks))
        int <- peaks$Peaks
    if(is.matrix(peaks))
        int <- peaks
    if(!is.matrix(int))
        stop("Invalid input object. Expecting a netCDF object or matrix")

    baseline <- apply(int, 2, base_fun, threshold, alpha, bfrac, segments, signalWindow, method)
    int <- int - baseline
    int[ int < 0 ] <- 0
    if(is.matrix(peaks))
        return(int)
    peaks$Peaks <- int
    peaks$baselineCorrected <- TRUE
    peaks
}

# high pass filter
# y[i] <- alpha * (y[i-1] + x[i] - x[i-1])

hpf <- function(x, alpha) {
        n <- length(x)
        .C(c_hpf, x = as.numeric(x), y = double(n), n = as.integer(n),
            a = as.numeric(alpha), PACKAGE="TargetSearch")$y
}

windowingStep <- function(x, n, wm) {
        x.len <- length(x)
        out <- .C(c_windowing, s = integer(n), x = as.integer(x), wm = as.integer(wm),
                n = as.integer(n), nidx = as.integer(x.len),  PACKAGE="TargetSearch" )
        as.logical(out$s)
}

#' baseline correction of MS data. (quantiles)
#'
#' Use the quantiles method to estimate the baseline
#'
#' @param peaks The netCDF object or a matrix
#' @param time the retention time in seconds. used if `peaks` is a matrix.
#' @param smooth integer. Smooth each signal by this number of points.
#'   Smoothing is disabled if this value is less or equal than 1.
#' @param qntl numeric scalar. The quantile for baseline estimation.
#'   The value must be in [0, 1].
#' @param width numeric scalar. The size of the window centered around
#'   a scan for baseline estimation.
#'   The size depends on the parameter `unit` below.
#' @param unit The width unit, which can be in seconds or points.
#' @param steps integer scalar greater than zero. To speed up computation, the baseline
#'   algorithm does not compute the estimate in each single scan, but in intervals of `steps`
#'   steps. The intermediate points are estimated by simple linear regression.
#' @return a netCDF object list or a baseline corrected matrix
baselineCorrectionQuant <- function(peaks, time, smooth=0, qntl=0.50, width=30, unit=c("seconds", "points"), steps=10)
{
    bslinefun <- function(x, ...) .Call(c_baseline, x, ...)

    assert_that(is.scalar(smooth), noNA(smooth))
    assert_that(is.scalar(qntl), noNA(qntl))
    assert_that(is.scalar(width), noNA(width))
    assert_that(is.count(steps))

    int <- NULL
    unit <- match.arg(unit)

    if(is.list(peaks)) {
        int <- peaks$Peaks
        t <- switch(unit, seconds=peaks$Time, points=seq(peaks$Time))
    }
    if(is.matrix(peaks)) {
        int <- peaks
        t <- if(unit == 'points') seq(nrow(int)) else time
    }
    if(!is.matrix(int))
        stop("Invalid input object. Expecting a netCDF object or matrix")

    # Compute the baseline
    bsline <- apply(int, 2, bslinefun, t, qntl, width, steps)

    int <- int - bsline
    int[ int < 0 ] <- 0
    smooth <- round(smooth)
    if(smooth > 1)
        int <- apply(int, 2, filter, rep(1, smooth) / smooth)

    if(is.matrix(peaks))
        return(int)

    keep <- rowSums(is.na(int)) == 0
    peaks$Peaks <- int[keep, ]
    peaks$Time  <- peaks$Time[ keep ]
    if(!is.null(peaks$Index))
        peaks$Index <- peaks$Index[ keep ]
    peaks$baselineCorrected <- TRUE
    peaks
}

# vim: set ts=4 sw=4 et:
