\name{FindAllPeaks}
\alias{FindAllPeaks}
\title{Extract peaks from chromatogram files - low level function}
\description{
  This function extracts all peaks of a given metabolite in a given
  RI window. It is intented to fine-tune metabolite search parameters.
}
\usage{
FindAllPeaks(samples, Lib, libID, dev=NULL, mz=NULL, RI=NULL,
         mz_type = c('selMass', 'quantMass', 'topMass'),
         columns = c("SPECTRUM", "RETENTION_TIME_INDEX", "RETENTION_TIME"))
}
\arguments{
  \item{samples}{ A \code{tsSample} object created by \code{\link{ImportSamples}} function. }
  \item{Lib}{ A \code{tsLib} object created by \code{\link{ImportLibrary}} function. }
  \item{libID}{ An index (integer or character) value representing the respective metabolite in
    the reference library \code{Lib}. }
  \item{dev}{ The allowed retention index (RI) deviation or \code{NULL}. }
  \item{mz}{A list of m/z values to search or \code{NULL}. }
  \item{RI}{ The expected retention index  or \code{NULL}. }
  \item{mz_type}{ whether to search for the selective, quant or top masses of the respective
    metabolite.}
  \item{columns}{ A numeric vector with the positions of the columns \code{SPECTRUM},
    \code{RETENTION_TIME_INDEX}, and \code{RETENTION_TIME} or a character vector
    with the header names of those columns.}
}
\details{
   The function searches for all peaks of a metabolite in all samples within a RI
   window. The parameters \code{dev}, \code{mz}, and \code{RI} have preference over
   the settings of the metabolite indexed by \code{libID}. In fact, if all of
   these parameters are not \code{NULL}, then \code{refLib} and \code{libID} are not
   used.

   The \code{columns} parameter is only needed for custom text RI files. There is
   no need to change it.
}
\note{
    This is an internal function not intended to be invoked directly, but
    it is exposed for convenience and advanced used.

    In the future it will replace \code{\link{FindPeaks}}.
}
\value{
   It returns a matrix in which each row represent a hit. Note that there can be zero
   rows if no hits are found. The columns are named and these are:
   \item{Int}{Peak intensity}
   \item{RI}{Retention Index}
   \item{RI}{Retention Time}
   \item{mz}{the searched m/z value}
   \item{fid}{the respective file or sample index. An integer value.}
}

\examples{
# load pre-calculated example data files and objects
require(TargetSearchData)
data(TSExample)

# get and set the RI file path
RIpath(sampleDescription) <- file.path(find.package("TargetSearchData"), "gc-ms-data")

# search all peaks of Valine (GC.3) and selective masses
peaks <- FindAllPeaks(sampleDescription, refLibrary, 'GC.3')
head(peaks)

# a numeric index is also allowed
peaks <- FindAllPeaks(sampleDescription, refLibrary, 3)
head(peaks)

# search arbitrary masses at arbitrary RI. the reference library and ID
# must be set set to NULL.
peaks <- FindAllPeaks(sampleDescription, NULL, NULL, dev=3000, RI=270000, mz=c(144, 100))
head(peaks)

}
\author{Alvaro Cuadros-Inostroza }
\seealso{\code{\link{FindPeaks}}}

% vim: set ts=4 sw=4 et:
