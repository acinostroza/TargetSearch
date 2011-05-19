\name{ProfileCleanUp}
\alias{ProfileCleanUp}
\title{ Reduce redundancy of the profile }
\description{
	This function reduces/removes redundancy in a profile.
}
\usage{
ProfileCleanUp(Profile, timeSplit=500, r_thres=0.95, minPairObs=5)
}
\arguments{
  \item{Profile}{ A \code{tsProfile} object. See \code{\link{Profile}}. }
  \item{timeSplit}{ A RI window. }
  \item{r_thres}{ A correlation threshold. }
  \item{minPairObs}{Minimum number of pair observations. Correlations between two variables are
                computed using all complete pairs of observations in those variables. If the number
                of observations is too small, you may get high correlations values just by chance,
                so this parameters is used to avoid that. Cannot be set lower than 5.}
}
\details{
	Metabolites that are inside a \code{timeSplit} window will be correlated
	to see whether the metabolites are the same or not, by using \code{r_thres}
	as a cutoff. If so, the intensities and RI will be averaged and the metabolite
	with more correlating masses will be suggested.
}
\value{
	A \code{tsProfile} object with a non-redundant profile of the masses that
	were searched and correlated, and intensity and RI matrices of the correlating masses.
	\item{slot "Info"}{A data frame with a profile of all masses that correlate and the metabolites
		that correlate in a \code{timeSplit} window.}
	\item{slot "profInt"}{A matrix with the averaged intensities of the correlating masses.}
	\item{slot "profRI"}{A matrix with the averaged RI of the correlating masses.}
	\item{slot "Intensity"}{A list containing peak-intensity matrices, one matrix per metabolite.}
	\item{slot "RI"}{A list containing RI matrices, one matrix per metabolite.}
}
\examples{
# load example data
require(TargetSearchData)
data(TargetSearchData)

RI.path <- file.path(.find.package("TargetSearchData"), "gc-ms-data")
refLibrary <- ImportLibrary(file.path(RI.path,"library.txt"))
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
metabProfile      <- Profile(sampleDescription, refLibrary, peakData, r_thres = 0.95) 

# here we use the metabProfile previously calculated and return a "cleaned" profile.
metabProfile.clean <- ProfileCleanUp(metabProfile, timeSplit = 500,
                      r_thres = 0.95) 

# Different cutoffs could be specified
metabProfile.clean <- ProfileCleanUp(metabProfile, timeSplit = 1000,
                      r_thres = 0.9) 

}
\author{ Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{Profile}}, \code{\linkS4class{tsProfile}} }
