#' interface to C function c_find_peaks
#'
#' @param file the RI file to read (text or bin)
#' @param mz the mass values to search (integer). It is coerced internally.
#' @param exp_ri the expected RI or NULL. It is only used for choosing the peak
#'      closer to the expected RI (instead of the maximum intensity)
#' @param min_ri the lower RI of the search window.
#' @param max_ri the upper RI of the search window.
#' @param useRT logical. Whether to use RT instead of RI. defaults to FALSE.
#' @param search integer. the search type. Values are 1 (all peaks) 2 (min RI
#'      error) or 3 (max intensity).
#' @param columns the columns for spectrum, RI and RT. Used only for text files.
#' @return a list with intensity, RI, RT and index of found values. The index
#'      is an integer that works as an index for the m/z values (zero-based).
`.c_find_peaks` <-
function(file, mz, exp_ri, min_ri, max_ri, useRT, search, columns)
{
    ad <- function(x) if(!is.null(x)) as.double(x) else x
    columns <- get.columns.name(columns)
    search <- c(all=1L, minRI=2L, maxInt=3L)[ search ]
    .Call(c_find_peaks, file, as.integer(mz), ad(exp_ri), ad(min_ri), ad(max_ri),
                        useRT, search, columns)
}

# low level search of peaks in a RI file.
`getAllPeaks` <-
function(file, ref, useRT=FALSE, searchType=c("all", "minRI", "maxInt"), columns=NULL)
{
    if(!all(c('mz', 'minRI', 'maxRI') %in% colnames(ref)))
        stop("Error: missing columns in 'ref'")

    searchType <- match.arg(searchType)
    assert_that(is.flag(useRT))

    z <- .c_find_peaks(file, ref[,'mz'], NULL, ref[, 'minRI'], ref[, 'maxRI'],
                       useRT, searchType, columns)
    assert_that(!is.null(z), msg=paste0('Error processing file ', file))
    z <- do.call('cbind', z)
    z[,4] <- z[,4] + 1
    colnames(z) <- c('Int', 'RI', 'RT', 'rowid')
    cbind(z, mz=ref[z[, 'rowid'],2])
}
