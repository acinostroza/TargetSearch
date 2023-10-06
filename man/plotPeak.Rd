\name{plotPeak}
\alias{plotPeak}
\title{ Plot peaks }
\description{
	Plot a peak of a given metabolite in a given sample showing the search windows.
}
\usage{
plotPeak(samples, Lib, metProf, rawpeaks, which.smp=1, which.met=1,
	massRange=NULL, corMass=FALSE)
}
\arguments{
  \item{samples}{A \code{tsSample} object created by \code{ImportSamples} function.}
  \item{Lib}{A \code{tsLib} object created by \code{\link{ImportLibrary}} function.}
  \item{metProf}{A metabolite profile object. See \code{\link{Profile}}}
  \item{rawpeaks}{ A three component list containing the retention time, the intensity matrix,
	and the mass range. See \code{\link{peakCDFextraction}}. }
  \item{which.smp}{ A numeric value indicating the sample. }
  \item{which.met}{ A numeric value indicating the metabolite. }
  \item{massRange}{ A two component numeric vector with the scan mass range to extract.
	Use \code{NULL} for automatic detection.}
  \item{corMass}{ Logical. If TRUE, show only correlating masses for the selected
	metabolite. Show all masses otherwise.}
}
\note{
    This function was completely rewritten. For the old function, see \code{\link{plotPeakSimple}}.

    In case the RI files are in text format and their column names are not standard (for
    example, when the files were generated with another software), use the global option
    \code{'TS_RI_columns'} or transform the RI files to \code{\link{TargetSearch}} binary
    format. See the documentation in function \code{\link{text2bin}}.
}
\value{
  The output is the same as the function \code{\link{peakCDFextraction}}, that is, a list
  containing the retention time (\code{Time}), the intensity matrix (\code{Peaks}),
  the retention index (\code{Index}), the m/z range (\code{massRange}), and a flag indicating
  whether the chromatogram was baseline corrected (\code{baselineCorrected}).

  The list can be recycled as the \code{rawpeaks} parameter for further plots (for example
  in a loop), so the CDF file doesn't need to be read again.
}
\examples{
require(TargetSearchData)
data(TSExample)

# update CDF and RI paths
CDFpath(sampleDescription) <- tsd_data_path()
RIpath(sampleDescription) <- tsd_data_path()

# Plot the peak "Valine" for sample number 1
grep("Valine", libName(refLibrary)) # answer: 3

# plot peak from the cdf file. The rawpeaks object can be recycled in order to plot
# other metabolites.
rawpeaks <- plotPeak(sampleDescription, refLibrary, metabProfile, which.smp=1,
        which.met=3, massRange=c(85,500), corMass=FALSE)
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{\code{\link{plotPeakSimple}}, \code{\link{RIcorrect}}, \code{\linkS4class{tsMSdata}}, \code{\linkS4class{tsRim}},
\code{\link{peakCDFextraction}}, \code{\link{matplot}} }
\keyword{ hplot }

