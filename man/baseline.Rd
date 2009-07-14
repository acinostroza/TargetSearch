\name{baselineCorrection}
\alias{baseline}
\alias{baselineCorrection}
\title{Baseline correction algorithm }
\description{
    Functions for baseline correction of GC-MS chromatograms.
}
\usage{
    baselineCorrection(int, threshold = 0.5, alpha = 0.95, bfraction = 0.2,
           segments = 100, signalWindow = 10, method = "linear")
    
    baseline(ncData, baseline.opts = NULL)
}
\arguments{
  \item{int}{A matrix object of spectra peak intensities to be baseline corrected,
  where the columns are retention times and rows mass traces.}
  \item{threshold}{A numeric value between 0 and 1. A value of one sets the baseline
  above the noise, 0.5 in the middle of the noise and 0 below the noise.}
  \item{alpha}{The alpha parameter of the high pass filter.}
  \item{bfraction}{The percentage of the fragments with the lowest intensities of
  the filtered signal that are assumed to be baseline signal.}
  \item{segments}{The number of segments in which the filtered signal is divided.}
  \item{signalWindow}{The window size (number of points) used in the signal windowing step.}
  \item{method}{The method used to approximate the baseline. \code{"linear"} (default) uses
  linear interpolation. \code{"spline"} fits a cubic smoothing spline (warning: really slow). }
  \item{ncData}{A list returned by the function \code{xcms:::netCDFRawData}}.
  \item{baseline.opts}{A list with parameters to be passed to \code{baselineCorrection}
  function. For example \code{baseline.opts = list(threshold = 0.5, alpha = 0.95)}.}
}
\details{
The baseline correction algorithm is based on the work of Chang et al, and it works as
follows. For every mass trace, i.e., rows of matrix \code{int}, the signal intensity is filtered
by a first high pass filter: \emph{y[i] = alpha * (y[i-1] + x[i] - x[i-1])}. The
filtered signal is divided into evenly spaced segments (\code{segments})
and the standard deviation of each segment is calculated. A percentage (\code{bfraction})
of the segments with the lowest values are assumed to be baseline signal and the
standard deviation (\emph{stdn}) of the points within those segments is calculated.

Once \emph{stdn} has been determined, the points with absolute filtered values larger than
\emph{2 * stdn} are considered signal. After that, the signal windowing step takes
every one of the points found to be signal as the center of a signal window (\code{signalWindow)}
and marks the points within that window as signal. The remaining points are now considered
to be noise.

The baseline signal is obtained by either using linear interpolation (default) or fitting a cubic
smoothing spline taking only
the noise. The baseline can be shifted up or down by using the parameter (threshold),
which is done by the formula: \emph{B' = B + 2*(threshold - 0.5)*2*stdn}, where
\emph{B} is the fitted spline, \emph{stdn} the standard deviation of the noise,
and \code{threshold} a value between 0 and 1. Finally, the corrected signal
is calculated by subtracting \emph{B'} to the original signal.
                                                                                      
The \code{baseline} function is called by the function \code{\link{NetCDFPeakFinding}}
before the peak picking algorithm is performed. Since it is an internal function,
it is not intended to be executed directly.
}
\value{
A matrix of the same dimensions of \code{int} with the baseline corrected intensities.
}
\references{
David Chang, Cory D. Banack and Sirish L. Shah, Robust baseline correction algorithm
for signal dense NMR spectra. \emph{Journal of Magnetic Resonance 187 (2007) 288-292}
}
\author{Alvaro Cuadros-Inostroza}
\seealso{ \code{\link{RIcorrect}} }
