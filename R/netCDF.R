# auxiliary functions to NetCDFPeakFinding and peakCDFextraction

# test whether the mass values are integers of double.
# if double, then checks whether they correspond with high mass
# accuracy data or not. If not, the mz values are rounded
# to the closes integer value.

.check.mz.precision <- function(x) {
	mz <- sample(x$mz, min(c(5000, length(x$mz))))
	z <- sum(abs(mz - round(mz)))
	w <- all.equal(z, 0L)
	if(is.logical(w)) {
		if(w) {
			return(x)
		}
	}

	k <- which.max(x$point_count)
	mz <- x$mz[ x$scanindex[k] + c(1:x$point_count[k])]
	z  <- sum(diff(mz) < 0.1)
	if( z / x$point_count[k] < 0.05 ) {
		res <- .cdffix(x)
		res$tic <- x$tic
		res$rt  <- x$rt
		return(res)
	}
	else {
		return(NULL)
	}
}

.cdffix <- function(x) {
	res <- .Call("cdffix", x$mz, x$intensity, x$scanindex, x$point_count,
		PACKAGE="TargetSearch")
	if(!is.null(res)) {
		res <- append(res, list(length(res[[3]]), length(res[[4]])))
		names(res) <- c("mz","intensity","scanindex", "point_count","ns","np")
	}
	res
}

# vim: set ts=4 sw=4:
