\name{plotPeakRI}
\alias{plotPeakRI}
\title{ Plot peak RI across samples }
\description{
    Plot peak RI of the quant mass of a given metabolite across samples.
    This function can be used to visualize the elution time of a metabolite
    in the RI and RT dimension, and in combination with the function
    \code{\link{sampleRI}} for fine tunning.
}
\usage{
plotPeakRI(samples, Lib, libID, dev=NULL, mz=NULL, RI=NULL,
           method=c('RI', 'Intensity'), useRI=TRUE, main=NULL,
           col=NULL, int_range=c(2,6), cex_range=c(.7,6), key_width=2)
}
\arguments{
  \item{samples}{A \code{tsSample} object created by \code{\link{ImportSamples}} function.}
  \item{Lib}{A \code{tsLib} object created by \code{\link{ImportLibrary}} function.}
  \item{libID}{ An index (integer or character) value representing the respective metabolite in
    the reference library \code{Lib}. }
  \item{dev}{ The allowed retention index (RI) deviation or \code{NULL}. }
  \item{mz}{A list of m/z values to search or \code{NULL}. }
  \item{RI}{ The expected retention index  or \code{NULL}. }
  \item{method}{ A character vector used to decided what peak should be chosen in case there are
      ambiguous peaks. If \code{'RI'}, then the closest peak to the expected RI is chosen. If
      \code{'Intensity'}, then the highest is taken. }
  \item{useRI}{Logical. Should the RI or RT be displayed in the y-axis? }
  \item{main}{ The title for the plot. If \code{NULL}, then the metabolite name is displayed. }
  \item{col}{ A color vector (length > 2) to show different levels of peak intensity. }
  \item{int_range}{ A length-two vector. The limits of the intensity for the color key. Note the
      the intensity is expressed in log10. Eg, a value of 2 represents a intensity of 100. }
  \item{cex_range}{ The 'cex' range value of the points. lower-intense peaks are represented
      as points of smaller size. }
  \item{key_width}{ The width in cm of the area allocated for the color key. }
}
\details{
    This function uses internally \code{\link{FindAllPeaks}}, so the same rules apply,
    i.e., the parameters \code{dev}, \code{mz}, and \code{RI} have preference over
    the settings of the metabolite indexed by \code{libID}.

    In the plot, the x-axis are samples as defined by the object \code{samples}. The
    y-axis is retention index (RI) as shown on the left-hand-size. On the right-hand-size
    y-axis the approximate retention time (RT) is shown. This is because the RT varies
    across samples, therefore it is averaged for displays purposes. If \code{useRI==FALSE},
    then the RT is displayed on the left hand size and the RI is averaged and shown on
    the left. Note that in either case, the RI is used for searching.

    The point size is proportional to the log10 of the peak intensity. Their size is
    controlled by the parameters \code{int_range} and \code{cex_range}. By default,
    intensities of 100 (log10 => 2) or lower are shown with \code{cex=0.7}, while intensities
    greater than 1000000 (log10 => 6) as displaued with \code{cex=6}. This also affects the
    scaling of the color key.

    The best peaks, selected acording to \code{method}, are shown with a black border,
    while the other are shown with no border and slightly transparent.

    The output is the RI of the best peaks or invisible. Note that if no peak is found,
    then no plot is drawn.
}
\value{
    Returns invisible or a numeric vector with the corresponding RI of the best peak chosen by
    \code{method}.
}
\examples{
def.par <- par(no.readonly = TRUE) # save parameters for resetting

# load pre-calculated example data files and objects
require(TargetSearchData)
data(TSExample)

# get and set the RI file path
RIpath(sampleDescription) <- tsd_data_path()

# search all peaks of Valine (GC.3) and selective masses. Retention index
ri <- plotPeakRI(sampleDescription, refLibrary, 'GC.3')

# increase deviation, change m/z to search, change colors and title
main <- 'Valine'
cols <- c('red', 'blue', 'green')
ri <- plotPeakRI(sampleDescription, refLibrary, 'GC.3', dev=4000, mz=144,
                 main=main, col=cols)

# plot by RT instead. Note the RI is still returned
ri <- plotPeakRI(sampleDescription, refLibrary, 'GC.3', useRI=FALSE)

par(def.par) # reset to default
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{\code{\link{FindAllPeaks}}, \code{\link{sampleRI}}, \code{\link{ImportSamples}}, \code{\link{ImportLibrary}}}

% vim: set ts=4 sw=4 et:
