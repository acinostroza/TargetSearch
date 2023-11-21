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
    result <- .c_find_peaks(RIfile, as.integer(ref[,'mz']), NULL, ref[,'minRI'],
                                  ref[,'maxRI'], useRT, 1L, ...)
    assert_that(!is.null(result), msg="An error occurred during processing")
    result <- cbind(RI=result[[2]], RT=result[[3]], Intensity=result[[1]],
                    mz=massValues[result[[4]] + 1])
    result
}

#
# R implementation for the above function, used for unit-testing
#

#' Extract peak data from a RI file (R version)
#'
#' A function to extract peak data from a RI file. This function is
#' basically identical to `ri_data_extract` except it uses pure R instead
#' of C. The idea is to use this function for unit testting
#'
#' @note Only text files are supported.
#'
#' @param RIfile A path to a RI file format (text only).
#' @param massValues A numeric vector representing m/z values.
#' @param timeRange A numeric vector of even length describing the lower and upper time limits
#' @param useRT Logical. If `TRUE`, the time range is in seconds.
#' @param ... extra parameters passed to internal functions
#' @return A four column matrix with columm names RI, RT, Intensity, mz.
#'
`ri_data_extract_text` <- function(RIfile, massValues, timeRange, useRT=FALSE, ...)
{
    assert_that(is.string(RIfile))
    assert_that(is.numeric(massValues))
    assert_that(is.numeric(timeRange), length(timeRange) %% 2 == 0)
    assert_that((length(massValues) %% (length(timeRange) / 2)) == 0)
    assert_that(is.flag(useRT))
    if(!is.matrix(timeRange)) timeRange <- matrix(timeRange, ncol=2)
    assert_that(ncol(timeRange) == 2)

    ref <- cbind(mz=massValues, minRI=timeRange[,1], maxRI=timeRange[,2])
    .r_find_peaks(RIfile, ref, useRT, ...)
}

#' R implementation of c_find_peaks, but only works on TXT files
#'
#' @param RIfile A path to a RI file format (text only)
#' @param ref a three-column matrix with m/z, min RI and max RI to search
#' @param useRT logical. if TRUE, search by RT, otherwise search by RI
#' @param ... potentially pass the column options
#' @return A four column matrix with columm names RI, RT, Intensity, mz.
#'
.r_find_peaks <- function(RIfile, ref, useRT, ...)
{
    search <- function(m, t_min, t_max) {
        t <- if(useRT) rt else ri
        j <- which(t > t_min & t < t_max)
        x <- vapply(sp[j], function(x) x[ match(m, x[,1]), 2], 0)
        n <- !is.na(x)
        mm <- if(any(n)) m else numeric(0)
        cbind(RI=ri[j][n], RT=rt[j][n], Intensity=x[n], mz=mm)
    }

    if(is.integer(cols <- get.columns.name(...)))
        cols <- cols + 1

    z <- read.delim(RIfile)
    sp <- Spectra( z[, cols[1] ] )
    ri <- z[, cols[2] ]
    rt <- z[, cols[3] ]
    rs <- mapply(search, ref[,'mz'], ref[,'minRI'], ref[,'maxRI'], SIMPLIFY=FALSE, USE.NAMES=FALSE)
    do.call('rbind', rs)
}

# vim: set ts=4 sw=4 et:
