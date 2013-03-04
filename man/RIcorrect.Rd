\name{RIcorrect}
\alias{RIcorrect}
\title{ Peak picking from CDF files and RI correction }
\description{
	This function reads from CDF files, finds the apex intensities, converts
	the retention time to retention time index (RI), and writes RI corrected
	text files (a.k.a. RI files).
}
\usage{
RIcorrect(samples, rimLimits = NULL, massRange = NULL, Window, IntThreshold,
	pp.method = "ppc", showProgressBar = FALSE, baseline = FALSE,
	baseline.opts = NULL )
}
\arguments{
  \item{samples}{ A \code{tsSample} object created by \code{ImportSamples} function. }
  \item{rimLimits}{ A \code{tsRim} object. If set to \code{NULL}, no retention
	time will be performed. See \code{ImportFameSettings}.}
  \item{massRange}{ A two component vector of m/z range used by the GC-MS machine
	or \code{NULL}, in which case the range will be detected automatically. }
  \item{Window}{ The window used for smoothing. The number of points actually used
	is \code{2*Window + 1}. }
  \item{IntThreshold}{ Apex intensities lower than this value will be removed from the RI files. }
  \item{pp.method}{ Peak picking method. Options are either "smoothing" or "ppc". See details. }
  \item{showProgressBar}{Logical. Should the progress bar be displayed?}
  \item{baseline}{Logical. Should baseline correction be performed?}
  \item{baseline.opts}{A list of options passed to \code{\link{baselineCorrection}}}.
}
\details{
	There are two pick picking methods available: "ppc" and "smoothing".

	The "ppc" method (default) implements the peak detection method described in the \code{ppc}
	package. It looks for the local maxima within a \code{2*Window + 1} scans for
	every mass trace.

	The "smoothing" method calculates a moving average of \code{2*Window + 1} points
	for every mass trace. Then it looks for a change of sign (from positive to negative) of
	the difference between two consecutive points. Those points will be returned as
	detected peaks.

	The RI file type is determined by the output of \code{\link{fileFormat}}
	method applied to the \code{\linkS4class{tsSample}} input object. To
	choose between the available formats ("binary" and "text"), select it
	with \code{\link{fileFormat}} method before calling \code{RIcorrect}.
}
\value{
	A retention time matrix of the found retention time markers. Every column represents a
	sample and rows RT markers.
}
\examples{
require(TargetSearchData)
# import refLibrary, rimLimits and sampleDescription.
data(TargetSearchData)
# get the CDF files
cdfpath <- file.path(find.package("TargetSearchData"), "gc-ms-data")
cdfpath
list.files(cdfpath)
# update the CDF path
CDFpath(sampleDescription) <- cdfpath
# run RIcorrect (massScanRange = 85-320; Intensity Threshold = 50;
# peak detection method = "ppc", window = 15)
RImatrix <- RIcorrect(sampleDescription, rimLimits, massRange = c(85,320), 
            Window = 15, pp.method = "ppc", IntThreshold = 50)

# you can try other parameters and other peak picking algorithm.
RImatrix <- RIcorrect(sampleDescription, rimLimits, massRange = c(85,320), 
            Window = 15, pp.method = "smoothing", IntThreshold = 10)

RImatrix <- RIcorrect(sampleDescription, rimLimits, massRange = c(85,320), 
            Window = 15, pp.method = "ppc", IntThreshold = 100)

}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{ImportSamples}}, \code{\link{ImportFameSettings}}, 
\code{\link{NetCDFPeakFinding}}, \code{\link{FAMEoutliers}},
\code{\linkS4class{tsSample}}, \code{\linkS4class{tsRim}}.}
