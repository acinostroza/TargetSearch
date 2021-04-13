#' Extract peak data from a RI file
#'
#' A convenient function to extract peak data from a RI file. This function is mostly
#' intended for debugging purposes.
#'
#' @param RIfile A path to a RI file format (binary or text).
#' @param massValues A numeric vector representing m/z values.
#' @param timeRange A numeric vector of even length describing the lower and upper time limits
#' @param useRT Logical. If `TRUE`, the time range is in seconds.
#' @param ... extra parameters passed to internal functions
#' @return A four column matrix with columm names RI, RT, Intensity, mz.
#'
`ri_data_extract` <- function(RIfile, massValues, timeRange, useRT=FALSE, ...)
{
    assert_that(is.string(RIfile))
    assert_that(is.numeric(massValues))
    assert_that(is.numeric(timeRange), length(timeRange) %% 2 == 0)
    assert_that((length(massValues) %% (length(timeRange) / 2)) == 0)
    assert_that(is.flag(useRT))
    if(!is.matrix(timeRange)) timeRange <- matrix(timeRange, ncol=2)
    assert_that(ncol(timeRange) == 2)

    ref <- cbind(mz=massValues, minRI=timeRange[,1], maxRI=timeRange[,2])
    opt <- as.integer(get.file.format.opt(RIfile, ...))
    result <- .Call(c_find_peaks, RIfile, as.integer(ref[,'mz']), NULL, ref[,'minRI'],
                                  ref[,'maxRI'], opt, useRT, 1L)
    result <- cbind(RI=result[[2]], RT=result[[3]], Intensity=result[[1]],
                    mz=massValues[result[[4]] + 1])
    result
}

# vim: set ts=4 sw=4 et:
