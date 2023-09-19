\name{RIcorrect}
\alias{RIcorrect}
\title{ Peak picking from CDF files and RI correction }
\description{
    This function reads from CDF files, finds the apex intensities, converts
    the retention time to retention time index (RI), and writes RI corrected
    text files (a.k.a. RI files). In addition, it can perform baseline
    correction and also convert files to the new NetCDF-4 TargetSearch format.
}
\usage{
RIcorrect(samples, rimLimits = NULL, massRange=NULL, Window, IntThreshold,
         pp.method="ppc", showProgressBar=FALSE, baseline=FALSE,
         writeCDF4path=TRUE, ...)
}
\arguments{
    \item{samples}{A \code{\linkS4class{tsSample}} object usually created by
        \code{\link{ImportSamples}} function.}
    \item{rimLimits}{ A \code{\linkS4class{tsRim}} object. If set to \code{NULL}, no retention
        time will be performed. See \code{ImportFameSettings}.}
    \item{massRange}{Deprecated. It is completely ignored but it is kept for
        compatibility with old scripts.}
    \item{Window}{The window used for smoothing. The number of points actually used
        is \code{2*Window + 1}. It must be an integer. See details.}
  \item{IntThreshold}{ Apex intensities lower than this value will be removed from the RI files. }
  \item{pp.method}{ Peak picking method. Options are "smoothing", "gaussian" and "ppc". See details. }
  \item{showProgressBar}{Logical. Should the progress bar be displayed?}
  \item{baseline}{Logical. Should baseline correction be performed?}
    \item{writeCDF4path}{Whether or not convert a CDF-3 into a CDF-4. It can
        take a logical value or a character vector representing file paths.
        See details below.}
  \item{\dots}{A list of options passed to \code{\link{baseline}}}.
}
\details{
	There are three pick picking methods available: "ppc", "smoothing", "gaussian".

	The "ppc" method (default) implements the peak detection method described in the \code{ppc}
	package. It looks for the local maxima within a \code{2*Window + 1} scans for
	every mass trace.

	The "smoothing" method calculates a moving average of \code{2*Window + 1} points
	for every mass trace. Then it looks for a change of sign (from positive to negative) of
	the difference between two consecutive points. Those points will be returned as
	detected peaks.

	The "gaussian" method behaves similar to the "smoothing" method, but instead
	a gaussian smoother is used instead of the moving average.

	To work out a suitable \code{Window} value, the following might be useful:
	\code{Window = (SR * PW - 1) / 2}, where SR is the scan rate of the MS instrument and
	PW is the peak width. Because \code{Window} is an integer, the resulting value must be
	rounded. For example, for SR = 20 scans per second, a PW = 1.5 seconds, then
	\code{Window} = 14.5, which can be rounded to 15.

	The RI file type is determined by the output of \code{\link{fileFormat}}
	method applied to the \code{\linkS4class{tsSample}} input object. To
	choose between the available formats ("binary" and "text"), select it
	with \code{\link{fileFormat}} method before calling \code{RIcorrect}.

    The parameter \code{writeCDF4path} is used to convert CDF-3 files into a
    custom CDF-4 format. It can be logical or a character vector. If \code{TRUE},
    the default, the CDF-3 files will be converted automatically to CDF-4
    format. The files will be saved in the same directory as the original CDF-3.
    If it is character vector representing file paths, then the CDF-4 will be
    saved in those paths instead (re-cycled to the length of samples).
    Finally, if \code{FALSE}, then no CDF conversion will be performed. This is not
    recommended, but if possible if you want to match the old \code{\link{TargetSearch}}
    behaviour. Note that also the baseline correction or the retention
    time indices will not be updated as well.

    If \code{baseline} is \code{TRUE}, the CDF files will be baseline corrected
    by passing the extra parameters to \code{\link{baseline}}. See there for
    details.
}
\value{
	A retention time matrix of the found retention time markers. Every column represents a
	sample and rows RT markers.
}
\examples{
require(TargetSearchData)
# import refLibrary, rimLimits and sampleDescription.
data(TSExample)
# get the CDF files
cdfpath <- tsd_data_path()
cdfpath
list.files(cdfpath)
# update the CDF path
CDFpath(sampleDescription) <- cdfpath

# change file format of RI files as bin
fileFormat(sampleDescription) <- 'binary'

# Parameters: Intensity Threshold = 50 peak detection method = "ppc", window = 15
# To match the old behavior, the do not create CDF-4 Files (not recommended)
RImatrix <- RIcorrect(sampleDescription, rimLimits, writeCDF4path=FALSE,
            Window = 15, pp.method = "ppc", IntThreshold = 50)

# Convert to CDF-4 (recommended) with same parameters
# Note: save files in same directory (as
RImatrix <- RIcorrect(sampleDescription, rimLimits, writeCDF4path=".",
            Window = 15, pp.method = "ppc", IntThreshold = 50)

# we need to update the sampleDescription to use the new files
# this is not done automatically
sampleDescription <- ncdf4Convert(sampleDescription, ".")

# you can try other parameters and other peak picking algorithm
RImatrix <- RIcorrect(sampleDescription, rimLimits,
            Window = 15, pp.method = "smoothing", IntThreshold = 10)

RImatrix <- RIcorrect(sampleDescription, rimLimits,
            Window = 15, pp.method = "ppc", IntThreshold = 100)
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{ImportSamples}}, \code{\link{ImportFameSettings}},
\code{\link{NetCDFPeakFinding}}, \code{\link{FAMEoutliers}},
\code{\linkS4class{tsSample}}, \code{\linkS4class{tsRim}}.}

% vim: set ts=4 sw=4 et:
