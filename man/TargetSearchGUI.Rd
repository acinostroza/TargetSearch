\name{TargetSearchGUI}
\alias{TargetSearchGUI}
\alias{TargetSearchGUI-deprecated}
\title{A GUI for TargetSearch}
\description{
    Opens a Graphical User Interface (GUI, written using Tcl/Tk) to allow easy
    setting and manipulation of most processing parameters which control GC-MS
    Data Evaluation with \emph{TargetSearch}.

    NOTE: This function has been deprecated. Use the command line instead.
}
\usage{
    TargetSearchGUI()
}
\details{
    The GUI is intended to facilitate the use of \emph{TargetSearch} for users unfamiliar
    with R otherwise. Many parameters that would be set calling the individual
    \emph{TargetSearch} Functions as described in the manual can be set here 'in one go'
    before running the complete analysis.
    
	Important Note: Please select the folder where you store your GC-MS Data
    (NetCDF or Apex) as the \code{Working Directory}. It is not yet possible to process
    data files from other/different locations.
	
	The parameters:
	
	\itemize{
	\item \code{Working Directory}: Use the \emph{Browse}-button to select the folder
    on your hard drive containing all your GC-MS data files. The output of
    \emph{TargetSearch} will be written to this folder too.

	\item \code{File Import}: Clicking \emph{NetCDF Data} or \emph{Apex Data} radio buttons will open
    a file select dialog. Choose the files you would like to be processed. You
    may check your selection pressing the \emph{Show}-button.

	\item \code{Baseline Correction}: Clicking \emph{on/off} button will perform
	baseline correction before peak detection. If selected, the \code{threshold}
	parameter is a numeric value between 0 and 1. A value of one returns a baseline
    above the noise, 0.5 in the middle of the noise and 0 below the noise. See
   \code{\link{baselineCorrection}} for further details.

	\item \code{Retention Index Correction}: Retention Index Correction is neccessary and
    applied only if you supply NetCDF Data (Apex Data contain already Retention
    Indices). You may \emph{Load} or \emph{Create} the search windows for your RI-Markers
    here. 

	\item \code{Peak Detection}: \emph{Search Windows} refers to the allowed RI deviation of your
    metabolites which are narrowed in 3 consecutive searches. \emph{Intensity Counts
    threshold} defines the minimum apex intensity incorporated in the analysis.
    A value of 1 would include all peaks. \emph{Mass Range} allows to limit the mass
    values (m/z) to be included in the analysis. \emph{Smoothing} averages raw data
    to eliminate some inherent noise leading to multiple peaks otherwise. 

	\item \code{Library}: A Library (to detect metabolites) usable by \emph{TargetSearch} contains
    at least information about the metabolite 'Name', its expected 'RI' and the
    selective masses in its spectrum 'SEL_MASS'. You may \emph{Load} or \emph{Create} one
    yourself using the respective buttons.
    The parameter \emph{no. of top masses} is the number of most intensive masses
    that will be taken from the spectrum, and \emph{excluded masses} is a list of
    masses that will be excluded.
    A more detailed description of the file formats can be found in \code{\link{ImportLibrary}}.

	\item \code{Normalization}: This selects how the data will be normalized during
    the metabolite search. Options are \code{"dayNorm"}, a day based median
	normalization, \code{"medianNorm"}, normalization using the median of all the intensities
	of a given mass, and \code{"none"}, no normalization at all.

	\item \code{Final Profiles}: Here you may set the parameters used by the functions
	\code{\link{Profile}} and \code{\link{ProfileCleanUp}}. \emph{timesplit} sets an
    RI window that will be used to look for metabolites that could have been redundanly
    identified. \emph{correl. thr.} is the correlation threshold and \emph{min. number
    of correlation samples} is a threshold used to make sure that correlations are
    computed with at least said number of observations.

	\item \code{Parameters}: You may \emph{Save} the current parameters as an \code{*.RData}
    file or \emph{Load} previously saved parameters to compare the outcome of different settings or
    just repeat the analysis.

	\item \code{Program}: \emph{Run} starts to process all currently selected files using the current
    parameters and saving output to \emph{Working Directory}. \emph{Quit} closes the GUI.
    }
}
\author{Jan Lisec}
