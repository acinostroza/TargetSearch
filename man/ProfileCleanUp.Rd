\name{ProfileCleanUp}
\alias{ProfileCleanUp}
\title{ Reduce redundancy of the profile }
\description{
	This function reduces/removes redundancy in a profile.
}
\usage{
ProfileCleanUp(Profile, timeSplit = 500, r_thres = 0.95)
}
\arguments{
  \item{Profile}{ A \code{tsProfile} object. See \code{\link{Profile}}. }
  \item{timeSplit}{ A RI window. }
  \item{r_thres}{ A correlation threshold. }
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
	\item{slot "Intensity"}{A matrix with the averaged intensities of the correlating masses.}
	\item{slot "RI"}{A matrix with the averaged RI of the correlating masses.}
}
\examples{
# load example data
require(TargetSearchData)
data(TargetSearchData)

# here we use the metabProfile previously calculated and return a "cleaned" profile.
metabProfile.clean <- ProfileCleanUp(metabProfile, timeSplit = 500,
                      r_thres = 0.95) 

# Different cutoffs could be specified
metabProfile.clean <- ProfileCleanUp(metabProfile, timeSplit = 1000,
                      r_thres = 0.9) 

}
\author{ Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{Profile}}, \code{\linkS4class{tsProfile}} }
