\name{baseline}
\alias{baseline}
\title{Baseline correction - wrapper function}
\description{
    This function perform baseline correction by wrapping around the methods
    implemented on \code{\link{baselineCorrection}} and \code{\link{baselineCorrectionQuant}}.
}
\usage{
    baseline(ncdf, bsline_method = c('classic', 'quantiles', 'none'), ...)
}
\arguments{
  \item{ncdf}{A list containing the raw chromatogram data. The list can be generated by
  \code{\link{peakCDFextraction}}. Expected elements are \code{"Peaks"} which is matrix
  of intensities where the rows are retention times and columns are mass traces, \code{"Time"}
  which is a vector of retention time in seconds.}
  \item{bsline_method}{A string to select the baseline retention method. Options
  are \code{"classic"} which implements Chang's method (the old or classic TargetSearch method),
  \code{"quantiles"} the quantiles based method, and \code{"none"} which does nothing
  (returns the same input).}
  \item{\dots}{Extra parameters to be passed to \code{\link{baselineCorrection}} or
  \code{\link{baselineCorrectionQuant}}}.
}
\details{
This is a wrapper function around the different baseline correction algorithms. It
is not intended to be executed by the average user. Please refer to the respective
man pages for details.}
\value{
  Returns a list with the same elements as the the input, but the element \code{"Peaks"}
  containing baseline corrected values.
}
\author{Alvaro Cuadros-Inostroza}
\seealso{ \code{\link{RIcorrect}}, \code{\link{baselineCorrection}}, \code{\link{baselineCorrectionQuant}} }
