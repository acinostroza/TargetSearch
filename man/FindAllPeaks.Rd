\name{FindAllPeaks}
\alias{FindAllPeaks}
\title{Extract peaks from chromatogram files - low level function}
\description{
  This function extracts all peaks of a given metabolite in a given
  RI window. This function is intended for fine-tuning metabolite search parameters.
}
\usage{
FindAllPeaks(samples, Lib, libID, dev=NULL, mz=NULL, RI=NULL,
         mz_type = c('selMass', 'quantMass', 'topMass'),
         columns = NULL)
}
\arguments{
  \item{samples}{ A \code{tsSample} object created by \code{\link{ImportSamples}} function. }
  \item{Lib}{ A \code{tsLib} object created by \code{\link{ImportLibrary}} function. }
  \item{libID}{ An index (integer or character) value representing the respective metabolite in
    the reference library \code{Lib}. It must be a scalar. }
  \item{dev}{ The allowed retention index (RI) deviation or \code{NULL}. See details below. }
  \item{mz}{A list of m/z values to search or \code{NULL}. }
  \item{RI}{ The expected retention index  or \code{NULL}. }
  \item{mz_type}{ whether to search for the selective, quant or top masses of the respective
    metabolite.}
  \item{columns}{Either \code{NULL}, a \code{character} vector, or an \code{integer} vector. In
    most cases, leave it as \code{NULL}. This parameter is used to configure the column names or
    positions of RI text files. See the documentation on the \code{\link{text2bin}} function for
    further details.}
}
\details{
   The function searches for all peaks of a metabolite in all samples within a RI
   window. The parameters \code{dev}, \code{mz}, and \code{RI} have preference over
   the settings of the metabolite indexed by \code{libID}. In fact, if all of
   these parameters are not \code{NULL}, then \code{refLib} and \code{libID} are not
   used.

   The \code{dev} parameter can either be \code{NULL} (as stated above), or a numeric vector.
   If the length is equal to 1, then the RI search window is the given \code{RI}
   plus or minus \code{dev}. If the length is 2, then the search window is from \code{RI + dev[1]}
   to \code{RI + dev[2]}. \emph{Note} than in this case \code{dev[1]} is usually a negative value.
   An error is raised if the length is greater than 2

   The \code{columns} parameter is only needed for custom text RI files. There is
   no need to change it.
}
\note{
    This is an internal function not intended to be invoked directly, but
    it is exposed for convenience and advanced users.

    In the future it may replace \code{\link{FindPeaks}}.

    See also the function \code{\link{ri_data_extract}} which offers a similar functionality
    but with different input parameters.
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

# an asymmetric deviation search
peaks <- FindAllPeaks(sampleDescription, refLibrary, 'GC.3', dev=c(-1000, 2000))
head(peaks)

# search arbitrary masses at arbitrary RI. the reference library and ID
# must be set set to NULL.
peaks <- FindAllPeaks(sampleDescription, NULL, NULL, dev=3000, RI=270000, mz=c(144, 100))
head(peaks)

}
\author{Alvaro Cuadros-Inostroza }
\seealso{
    \code{\link{FindPeaks}},
    \code{\link{ri_data_extract}},
}

% vim: set ts=4 sw=4 et:
