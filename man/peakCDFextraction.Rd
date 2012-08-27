\name{peakCDFextraction}
\alias{peakCDFextraction}
\title{ NetCDF to R }
\description{
	This function reads a netcdf chromatogram file and returns a list containing
	the retention time and the intensity matrices.
}
\usage{
peakCDFextraction(cdfFile, massRange = NULL)
}
\arguments{
  \item{cdfFile}{ A character string naming a netcdf file. }
  \item{massRange}{ A two component numeric vector with the scan mass range to extract. }
}
\details{
	The function expects the following NetCDF variables: \code{intensity_values},
	\code{mass_values}, \code{scan_index}, \code{point_count} and \code{scan_acquisition_time}.
	Otherwise, an error will be displayed.\cr
	
	The \code{massRange} parameter is a numeric vector with two components: lower
	and higher masses. All masses in that range will be extracted. Note that it is
	not possible to extract a discontinuous mass range. It also can be \code{NULL},
	in which case the mass range will be determined by scanning the cdf file.
}
\value{
	A two component list.
  \item{Time}{The retention time vector.}
  \item{Peaks}{The intensity matrix. Rows are the retention times and columns are masses. The first
	column is the lower mass value and the last one is the higher mass. }
  \item{massRange}{The mass range.}
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\note{
	This function does not look for peaks, just extracts all the raw intensity values
	of the chromatogram file. Use \code{\link{NetCDFPeakFinding}} instead.
}
\seealso{ \code{\link{NetCDFPeakFinding}} }
