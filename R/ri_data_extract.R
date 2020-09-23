#' Extract peak data from a RI file
#'
#' A convenient function to extract peak data from a RI file. This function is mostly
#' intended for debugging purposes.
#'
#' @param RIfile A path to a RI file format (binary or text).
#' @param massValues A numeric vector representing m/z values.
#' @param timeRange A numeric vector of length 2 (lower and upper time limits)
#' @param useRT Logical. If `TRUE`, the time range is in seconds.
#' @param ... extra parameters passed to internal functions
#' @return A four column matrix with columm names RI, RT, Intensity, mz.
#'
`ri_data_extract` <- function(RIfile, massValues, timeRange, useRT=FALSE, ...)
{
    assert_that(is.numeric(timeRange), length(timeRange) == 2)
    assert_that(is.numeric(massValues))
    ref <- cbind(minRI=timeRange[1], mz=massValues, maxRI=timeRange[2])
    result <- getAllPeaks(RIfile, ref, useRT, ...)
    result <- result[, c('RI', 'RT', 'Int', 'mz'), drop=FALSE]
    colnames(result)[3] <- 'Intensity'
    result
}

# vim: set ts=4 sw=4 et:
