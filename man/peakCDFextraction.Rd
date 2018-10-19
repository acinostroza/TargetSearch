\name{peakCDFextraction}
\alias{peakCDFextraction}
\title{ NetCDF to R }
\description{
	This function reads a netcdf chromatogram file and returns a list containing
	the retention time and the intensity matrices.
}
\usage{
peakCDFextraction(cdfFile, massRange)
}
\arguments{
  \item{cdfFile}{ A character string naming a netcdf file. }
  \item{massRange}{Deprecated. The mass range is extracted automatically. }
}
\details{
	The function expects the following NetCDF variables: \code{intensity_values},
	\code{mass_values}, \code{scan_index}, \code{point_count} and \code{scan_acquisition_time}.
	Otherwise, an error will be displayed.\cr

	The \code{massRange} parameter used to be a numeric vector with two components: lower
	and higher masses. Now, the mass rangeis determined by scanning the CDF file. All masses
	in the CDF file will be extracted.
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
