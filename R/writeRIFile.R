#' Function to write RI files
#'
#' Internal function to write the so-called RI files which contains the
#' extracted peaks from a chromatogram. The file format is binary (custom for
#' Targetsearch) or text (an old tab-delimited format). The actual writing
#' process is performed by a C function.
#'
#' @param outFile a file path in which the data will be stored (string)
#' @param peaks a list with two named components:
#'    - 'Peaks': a matrix where a non-zero value indicates a peak on a scan
#'      (row) on a m/z value (column).
#'    - 'Time': the RT of each scan. Its length matches the rows of 'Peaks'
#' @param riInde the retention time index RI. Its length should be equal to the
#'    'Time' component.
#' @param massRange a two-compenent vector containing the m/z range. The number
#'    of columns of 'Peaks' should match the m/z range (one nominal mass per
#'    column.
#' @param ftype a string to indicate the file format: 'binary' or 'text'
#' @return TRUE on success, FALSE otherwise. A message is printed if an error
#'    happens on the C side, but no error is raised. In this case, is up to
#'    the caller to check for its return value.
#' @note If the input paremeters are invalid, an assertion error is thrown.
`writeRIFile` <-
function(outFile, Peaks, riInde, massRange, ftype=c("binary", "text"))
{
    op <- switch(match.arg(ftype),
                    text=list(type=1L, head=get.file.header()),
                    binary=list(type=0L, head=NULL))
    assert_that(is.string(outFile))
    assert_that(is.numeric(Peaks$Peaks))
    assert_that(is.numeric(riInde), is.numeric(massRange))
    assert_that(is.matrix(Peaks$Peaks))
    assert_that(length(massRange) == 2)
    assert_that(ncol(Peaks$Peaks) == diff(massRange) + 1)
    assert_that(length(riInde) == nrow(Peaks$Peaks))
    assert_that(length(riInde) == length(Peaks$Time))

    ret <- .Call(c_write_peaks, outFile, Peaks$Time, riInde,
               as.integer(Peaks$Peaks), as.integer(massRange), op$type, op$head)
    invisible(ret == 0)
}
