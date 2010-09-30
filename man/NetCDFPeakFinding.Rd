\name{NetCDFPeakFinding}
\alias{NetCDFPeakFinding}
\title{ Peak picking algorithm from CDF files }

\description{
	This function reads a netcdf chromatogram file, finds the apex intensities
	and returns a list containing the retention time and the intensity matrices.
}

\usage{
NetCDFPeakFinding(cdfFile, massRange = c(85, 500), Window = 5, IntThreshold = 10,
                  pp.method = "ppc", baseline = FALSE, baseline.opts = NULL)
}

\arguments{
  \item{cdfFile}{A character string naming a netcdf file.}
  \item{massRange}{A two component numeric vector with the scan mass range to extract.}
  \item{Window}{The window used by peak picking method. The number of points actually used
		is \code{2*Window + 1}.  }
  \item{IntThreshold}{Apex intensities lower than this value will be removed from the RI files. }
  \item{pp.method}{The pick picking method to be used. Options are \code{"smoothing"} and \code{"ppc"}.}
  \item{baseline}{Logical. Should baseline correction be performed?}
  \item{baseline.opts}{A list of options passed to \code{\link{baselineCorrection}}}.
}

\details{
	The function expects the following NetCDF variables: \code{intensity_values},
	\code{mass_values}, \code{scan_index}, \code{point_count} and \code{scan_acquisition_time}.
	Otherwise, an error will be displayed.
	
	The \code{massRange} parameter is a numeric vector with two components: lower
	and higher masses. All masses in that range will be extracted. Note that it is
	not possible to extract a discontinuous mass range.
	
	There are two peak picking algorithms that can be used. The \code{"smoothing"} method
	smooths the m/z curves and then looks for a change of sign of the intensity difference
	between two consecutive points. The \code{"ppc"} uses a sliding window and looks
	for the local maxima. This method is based on R-package \code{ppc}.

}

\value{
	A two component list.
  \item{Time}{The retention time vector.}
  \item{Peaks}{The intensity matrix. Rows are the retention times and columns are masses. The first
	column is the lower mass value and the last one is the higher mass. }
}

\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\examples{
require(TargetSearchData)
data(TargetSearchData)
CDFpath <- file.path(.find.package("TargetSearchData"), "gc-ms-data")
CDFfiles <- dir(CDFpath, pattern = ".cdf$", full.names = TRUE)
CDFfiles

# extrac peaks of first chromatogram
peaks.1 <- NetCDFPeakFinding(CDFfiles[1], massRange = c(85, 320), Window = 15,
			 IntThreshold = 10, pp.method = "smoothing")
# scan acquisition times
head(peaks.1$Time)
# peaks in matrix form. first column is mass 85, last one is mass 320.
head(peaks.1$Peaks)

}
\seealso{ \code{\link{peakCDFextraction}} }
