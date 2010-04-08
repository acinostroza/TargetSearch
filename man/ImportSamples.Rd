\name{ImportSamples}
\alias{ImportSamples}
\alias{ImportSamplesFromDir}
\title{ Sample definitions }
\description{
	This function imports a sample list that will be processed from a tab
	delimited file.
}
\usage{
ImportSamples(sampfile, CDFpath = ".", RIpath = ".", ...)

ImportSamplesFromDir(CDFpath = ".", RIfiles = FALSE, ignore.case = TRUE)
}
\arguments{
  \item{sampfile}{ A character string naming a sample file. See details. }
  \item{CDFpath}{ A character string naming a directory where the CDF files are
	located. }
  \item{RIpath}{ A character string naming a directory where the RI corrected text
	files are/will be located. }
  \item{RIfiles}{ Logical. If \code{TRUE}, the function will look for for RI files
  (RI_*) instead of CDF files (the default).}
  \item{ignore.case}{Logical. Should pattern-matching be case-insensitive?}
  \item{\dots}{ Other options passed to \code{\link[utils]{read.delim}} function. }
}
\details{
	The sample file is a tab-delimited text file with at least two columns:
	
	\itemize{
	\item \code{CDF_FILE} - The list of baseline corrected CDF files.
  \item \code{MEASUREMENT_DAY} - The day when the sample was measured.
	}
	  
  The column names must be exactly those indicated, but the column order doesn't
	matter. Other columns could be included in that file. They won't be used by the script,
	but will be included in the sample R object.
}
\value{
	A \code{tsSample} object.
}
\examples{
# get the sample definition definition file
cdfpath <- file.path(.find.package("TargetSearchData"), "gc-ms-data")
sample.file  <- file.path(cdfpath, "samples.txt")

# set a path where the RI files will be created
RIpath <- "."

# import samples
sampleDescription <- ImportSamples(sample.file, CDFpath = cdfpath, RIpath = RIpath)

# change the sample names
sampleNames(sampleDescription) <- paste("Sample", 1:length(sampleDescription), sep = "_")

# change the file paths (relative to the working path)
CDFpath(sampleDescription) <- "my_cdfs/"
RIpath(sampleDescription)  <- "my_RIs/"

}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{ImportLibrary}}, \code{\linkS4class{tsSample}} }
