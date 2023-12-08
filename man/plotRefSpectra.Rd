\name{plotRefSpectra}
\alias{plotRefSpectra}
\title{ Plot reference spectrum of a compound }
\description{
    A simple function to plot the reference spectrum of a compound in the style of
    common visualization tools. This is useful when comparing the spectrum of a metabolite
    with such software tools.
}
\usage{
plotRefSpectra(lib, libID, col='darkblue', main=NULL, xlab='m/z', ylab='Intensity',
               label=TRUE, cex=0.8, label_col=NULL, sel_font=4, type='h', ...)
}
\arguments{
    \item{lib}{ A \code{tsLib} object created by \code{\link{ImportLibrary}} function.}
    \item{libID}{ An index (integer or character) value representing the respective metabolite in
        the reference library \code{lib}. If missing, the first metabolite in the library
        is taken. }
    \item{col}{ The color(s) of the vertical lines. This option is passed to
        \code{\link[base]{plot}}. }
    \item{main}{ The title for the plot. If \code{NULL}, then the metabolite name is displayed. }
    \item{xlab}{ The \emph{x} axis label. Passed to \code{\link[base]{plot}}. }
    \item{ylab}{ The \emph{y} axis label. Passed to \code{\link[base]{plot}}. }
    \item{label}{ Logical. It controls whether the \emph{m/z} labels are displayed. }
    \item{cex}{ The \code{cex} values (numerical vector) used to scale the \emph{m/z} labels of
        the plot. }
    \item{label_col}{ The color vector for the \emph{m/z} labels. }
    \item{sel_font}{ The font for the \emph{m/z} labels of the \emph{selective} masses. This should
        be an integer like the parameter \code{font} documented in \code{\link[graphics]{par}}. }
    \item{type}{ A character that specifies the plot type. Defaults to \code{'h'} (for vertical
        bars) and \strong{cannot} be changed by the user (it has no effect). }
    \item{\dots}{ Extra graphical parameters passed to \code{\link[base]{plot}}. }
}
\details{
    This is a simple function to plot the spectrum of a reference metabolite. It is expected
    that the spectrum for that metabolite is not empty, in which case nothing is plotted
    and a warning is thrown.

    The \emph{m/z} label positions are determined automatically such that they do not intersect
    with each other and do not touch the vertical lines of the plot.
}
\value{
    Returns \code{invisible()}.
}
\examples{
def.par <- par(no.readonly = TRUE) # save parameters for resetting

# load pre-calculated example data files and objects
data(TSExample)

# plot the reference spectrum of the compound ID 'GC.3'
plotRefSpectra(refLibrary, 'GC.3')

# that is equivalent to
plotRefSpectra(refLibrary['GC.3'])

# change some graphical parameter and add a grid via 'panel.first'
plotRefSpectra(refLibrary, 'GC.3', col='darkred', main='Valine', xlab='mz',
               ylab='intensity', cex=0.7, label_col='blue', sel_font=2,
               panel.first=grid())

par(def.par) # reset to default
}
\author{Alvaro Cuadros-Inostroza}
\seealso{See the function \code{\link{ImportLibrary}} for importing a metabolite library.
    See \code{\link[base]{plot}} and \code{\link[graphics]{par}} for graphical parameters.}

% vim: set ts=4 sw=4 et:
