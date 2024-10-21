\name{FAMEoutliers}
\alias{FAMEoutliers}
\title{ FAME outlier detection }
\description{
    A function to detect retention time marker (FAME) outliers.
}
\usage{
FAMEoutliers(samples, RImatrix, pdffile = NA, startDay = NULL, endDay = NULL,
             threshold = 3, group.threshold = 0.05)
}
\arguments{
  \item{samples}{ A \code{tsSample} object created by \code{ImportSamples} function. }
  \item{RImatrix}{ A retention time matrix of the found retention time markers. }
  \item{pdffile}{ A character string naming a PDF file where the FAMEs report
    will be saved. }
  \item{startDay}{ A character vector with the starting days of your day groups. }
  \item{endDay}{ A character vector with the ending days of your day groups. }
  \item{threshold}{ A standard deviations cutoff to detect outliers. }
  \item{group.threshold}{ A numeric cutoff to detect day groups based on hierarchical
    clustering. Must be between \code{0..1}.}
}
\details{
    If no \code{pdffile} argument is given, the report will be saved on a file called
    \code{"TargetSearch-YYYY-MM-DD.FAME-report.pdf"}, where \code{YYYY-MM-DD} is a
    date.

    If both \code{startDay} and \code{endDay} are not given (both set to either
    \code{NULL} or \code{NA}), the function will try
    to detect day groups using a hierarchical clustering approach by cutting the tree
    using \code{group.threshold} as cutoff height. Otherwise, both must have the
    same length, must not contain \code{NA}s, and must match the measurement days
    of the object \code{samples}. See example below.

    Retention time markers that deviate more than \code{threshold} standard
    deviations from the    mean of their day group will be identified as outliers.
}
\value{
    A logical matrix of the same size of \code{RImatrix}. A \code{TRUE} value indicates
    that the retention time marker in that particular sample is an outlier.
}
\examples{
# load pre-calculated example data and objects
data(TSExample)

# find the retention marker outliers of the example data and save it in "outlier.pdf"
outliers <- FAMEoutliers(sampleDescription, RImatrix, pdffile = "outlier.pdf")

# find the outliers (although they are reported in the output PDF file)
apply(outliers, 1, which)

# using `startDay` and `endDay` to set day groups manually
sampleDays(sampleDescription) <- rep(letters[1:3], each=5) # make fake days a, b, c

# set days 'a' and 'b' in a group and day 'c' by itself
startDay <- c('a', 'c')
endDay <- c('b', 'c')

outliers <- FAMEoutliers(sampleDescription, RImatrix, pdffile="outlier.pdf",
            startDay=startDay, endDay=endDay)
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{RIcorrect}}, \code{\link{ImportSamples}}, \code{\link{TSExample}} }
