% Please see refer to file R/checkRimLim.R
\name{checkRimLim}
\alias{checkRimLim}
\title{Visually check retention index marker limits}
\usage{
checkRimLim(samples, rim, layout, show = TRUE, single = TRUE, extend = 0.5,
  rect.col = "#e7e7e7", mar = c(2, 2, 2, 2), oma = c(3, 3, 2, 0.5),
  cex.main = 1, type = "l", ...)
}
\arguments{
\item{samples}{A  \code{\linkS4class{tsSample}} object created by
\code{\link{ImportSamples}}.}

\item{rim}{A \code{\linkS4class{tsRim}} object describing the retention index
markers. See \code{\link{ImportFameSettings}}.}

\item{layout}{A vector of the form \code{c(nr, nc)} to arrange the panel by
\code{nr} rows and \code{nc} columns. If missing then the layout is created
automatically.}

\item{show}{Logical. If \code{FALSE} the plot is not shown, but the data points
can be used for further inspection or for custom plots.}

\item{single}{Logical. If \code{TRUE} a \emph{single} sample will be selected randomly
for plotting. This was the old default behavior. If \code{FALSE} all samples
will be used for plotting (but see note below).}

\item{extend}{a numeric coefficient to extend the time window search of the
respective time marker. Defaults to \code{0.5}.}

\item{rect.col}{the color for the background rectangle which indicates
the current retention time limits.}

\item{mar}{the subplots margins, passed to \code{\link[graphics]{par}()}.}

\item{oma}{the outer plot margins, passed to \code{\link[graphics]{par}()}.}

\item{cex.main}{The magnification to be used for main titles, passed to
\code{\link[graphics]{par}()}.}

\item{type}{A character vector indicating the type of plots. Default \code{"l"} for
lines. Passed to \code{\link[graphics]{matplot}()}.}

\item{\ldots}{extra plotting arguments passed to \code{\link[graphics]{matplot}()}
such as \code{col}, \code{lty}, \code{pch}, \code{lwd}.}
}
\value{
    The output value is either \code{invisible} or it depends on the option \code{single}.
    If this option is \code{TRUE} (the old behavior), then the output will be a list
    of \code{n} times \code{2} matrices, each element corresponding to a retention marker.
    Columns are retention time and intensities of the respective marker's \emph{m/z}. The rows
    can be as many data points are within the search window.

    If \code{single=FALSE}, the output is a list whose length is equal
    to \code{length(samples)} and its names are equal to \code{sampleNames(samples)}.
    Each element is in turn a list with exactly the same structure described in the
    paragraph above.
}
\description{
    A function to visually check if the retention time search limits of the
    retention index markers (aka FAMEs) are correct or need some adjustment.
}
\details{
    The function takes a \code{\linkS4class{tsSample}} object and creates a panel
    plot of the m/z traces around the area in which a marker is expected to be,
    one panel for each marker. By default, a single (random) sample is chosen (see
    option \code{single}) for plotting, but it is also possible to visualize many
    samples at the same time (but see note below).

    For multiple sample visualization, it is recommended to use sub-setting for the
    \code{samples} and \code{rim} arguments in order to avoid over crowded plots,
    specially when there are several samples and several retention index markers.
    See the examples below.

    Plotting options such as \code{col}, \code{lty}, \code{pch}, \code{lwd} can
    be passed; these are in turn passed to \code{\link[graphics]{matplot}()}.
    Note that these are passed to each panel; it not possible to pass different
    options to different panels (for example to have different line colors).
    Moreover, using options such as \code{xlim} and \code{ylim} may result in
    \sQuote{empty} plots (because each panel needs its own limits). If different
    styles are required, then either make your own function from the output or
    use subsets of the \code{rim} (\code{\linkS4class{tsRim}}) object.
}
\note{
    If \code{single=TRUE}, all \code{CDF} files will be scanned, which can
    take a significant amount of time to parse, in particular, if there are
    hundred of them. Also due to this, each panel of the plot could be extremely
    crowded. It is therefore recommended to use sample sub-setting to reduce
    this number, as shown in the examples below.
}
\examples{
require(TargetSearchData)

# get the cdf path TargetSearchData
cdfpath <- file.path(find.package("TargetSearchData"), "gc-ms-data")

# import samples (see ImportSamples() for details)
samples <- ImportSamples(file.path(cdfpath, "samples.txt"), CDFpath = cdfpath)

# Import RI markers (see ImportFameSettings())
rim <- ImportFameSettings(file.path(cdfpath, "rimLimits.txt"))

# choose a sample at random and plot the m/z traces around the retention time window
ret <- checkRimLim(samples, rim)

# to choose a specific sample and marker, use subsetting
ret <- checkRimLim(samples[3], rim[1:2])

# to display many samples at the same time, set the option `single` to `FALSE` and select a subset
# of samples (recommended). In this example the first three samples are chosen.
ret <- checkRimLim(samples[1:3], rim, single=FALSE)

# in general, visualizing all samples (as shown below) at the same time it is not recommended
# for large number of samples.
\dontrun{
  ret <- checkRimLim(samples, rim, single=FALSE)
}
}
\seealso{
  \code{\link[graphics]{matplot}} for plotting parameters,
  \code{\link[graphics]{par}} for inner and outer margins parameters,
  \code{\linkS4class{tsSample}},
  \code{\link{sampleNames}},
  \code{\linkS4class{tsRim}},
  \code{\link{ImportFameSettings}}
}
