\name{Profile}
\alias{Profile}
\title{ Average the correlating masses for each metabolite }
\description{
  This function makes a profile from the masses that correlate for each metabolite.
}
\usage{
Profile(samples, Lib, peakData, r_thres = 0.95, method = "dayNorm", minPairObs = 5)
}
\arguments{
  \item{samples}{A \code{tsSample} object created by \code{ImportSamples} function. }
  \item{Lib}{ A \code{tsLib} object created by \code{\link{ImportLibrary}} function with corrected RI values.
		See \code{\link{medianRILib}}. }
  \item{peakData}{ A \code{tsMSdata} object. See \code{\link{peakFind}}.}
  \item{r_thres}{ A correlation threshold. }
  \item{method}{ Normalisation method. Options are \code{"dayNorm"}, a day based median
	normalisation, \code{"medianNorm"}, normalisation using the median of all the intensities
	of a given mass, and \code{"none"}, no normalisation at all. }
  \item{minPairObs}{Minimum number of pair observations. Correlations between two variables are
		computed using all complete pairs of observations in those variables. If the number
		of observations is too small, you may get high correlations values just by chance,
		so this parameters is used to avoid that. Cannot be set lower than 5.}
}
\value{
	A \code{tsProfile} object. The slots are:
	\item{Info}{A data frame with a profile of all masses that correlate.}
	\item{Intensity}{A list containing peak-intensity matrices, one matrix per metabolite.}
	\item{RI}{A list containing RI matrices, one matrix per metabolite.}
	\item{profInt}{A matrix with the averaged intensities of the correlating masses.}
	\item{profRI}{A matrix with the averaged RI of the correlating masses.}
}
\examples{
require(TargetSearchData)
data(TargetSearchData)

# get RI file path
RI.path <- file.path(.find.package("TargetSearchData"), "gc-ms-data")
# update RI file path
RIpath(sampleDescription) <- RI.path
# Import Library
refLibrary        <- ImportLibrary(file.path(RI.path,'library.txt'))
# update median RI
refLibrary        <- medianRILib(sampleDescription, refLibrary)
# get the sample RI
corRI             <- sampleRI(sampleDescription, refLibrary, r_thres = 0.95)
# obtain the peak Intensities of all the masses in the library
peakData          <- peakFind(sampleDescription, refLibrary, corRI)
# make a profile of the metabolite data
metabProfile      <- Profile(sampleDescription, refLibrary, peakData, r_thres = 0.95)

# same as above, but with different thresholds.
metabProfile      <- Profile(sampleDescription, refLibrary, peakData,
                     r_thres = 0.9, minPairObs = 5)

}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{\code{\link{ImportSamples}}, \code{\link{ImportLibrary}}, \code{\link{medianRILib}},
	 \code{\link{peakFind}}, \code{\linkS4class{tsProfile}}  }
