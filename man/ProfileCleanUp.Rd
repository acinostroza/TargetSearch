\name{ProfileCleanUp}
\alias{ProfileCleanUp}
\title{ Reduce redundancy of the profile }
\description{
	This function reduces/removes redundancy in a profile.
}
\usage{
ProfileCleanUp(Profile, timeSplit=500, r_thres=0.95, minPairObs=5,
    prioritization=c('mass','score'), corMass=1, score=0,
    show=c('unidentified','knowns','full'))
}
\arguments{
  \item{Profile}{ A \code{tsProfile} object. See \code{\link{Profile}}. }
  \item{timeSplit}{ A RI window. }
  \item{r_thres}{ A correlation threshold. }
  \item{minPairObs}{Minimum number of pair observations. Correlations between two variables are
                computed using all complete pairs of observations in those variables. If the number
                of observations is too small, you may get high correlations values just by chance,
                so this parameters is used to avoid that. Cannot be set lower than 5.}
  \item{prioritization}{Selects whether the metabolite suggestion should be based on
                the number of correlation masses (\code{mass}) or the score (\code{score}).}
  \item{corMass}{Metabolites with a number of correlation masses lower than \code{score}
                will be marked as 'Unidentified RI'}
  \item{score}{Metabolites with a score lower than \code{score} will be marked as
                unidentified.}
  \item{show}{A character vector. If \code{unidentified}, all non-redundant metabolites will
                be returned; if \code{knowns}, only returns those metabolites with
                correlation masses and score greater than the given values; and if
                \code{full}, it shows all redundant metabolites, which may be useful
                to retrieve the data from misidentified metabolites. }
}
\details{
	Metabolites that are inside a \code{timeSplit} window will be correlated
	to see whether the metabolites are potentially the same or not, by using
	\code{r_thres} as a cutoff. If so, the best candidate will be chosen
	according to the value of \code{prioritization}: If 'mass', then
	metabolites will be suggested based on number of correlating masses, and
	if 'score', then the score will be used. Metabolites that don't have
	al least \code{corMass} correlating masses and \code{score} score will
	be marked as 'unidentified' and not will be suggested, unless all the
	metabolites in group are unidentified.

	For example, suppose that three metabolites A (CM=3, S=900),
	B (CM=6, S=700), C (CM=5, S=800) correlate within the same time group,
	where CM is the number of correlating masses and S is the score.

	\itemize{
	\item If prioritization='mass', corMass=3, score=650, then the
	suggested order is B, C, A.
	\item If prioritization='mass', corMass=3, score=750, then the
        suggested order is C, A, B.
	\item If prioritization='mass', corMass=3, score=850, then the
        suggested order is A, B, C.
	\item If prioritization='score', corMass=3, score=650, then the
        suggested order is A, C, B.
	\item If prioritization='score', corMass=4, score=650, then the
        suggested order is C, B, A.
	\item If prioritization='score', corMass=4, score=850, then the
        suggested order is C, A, B.
	}
	Note that by chosing \code{prioritization='mass'}, \code{score=0}, and
	\code{corMass=1} you will get the former behavior (\code{TargetSearch <= 1.6}).
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

RI.path <- file.path(find.package("TargetSearchData"), "gc-ms-data")
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
