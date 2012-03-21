\name{ImportSamples}
\alias{ImportSamples}
\alias{ImportSamplesFromDir}
\title{ Sample definitions }
\description{
	This function imports a sample list that will be processed from a tab
	delimited file.
}
\usage{
ImportSamples(sampfile, CDFpath=".", RIpath=".", ftype=c("binary", "text"), ...)

ImportSamplesFromDir(CDFpath=".", RIfiles=FALSE, ignore.case=TRUE,
	ftype=c("binary", "text"))
}
\arguments{
  \item{sampfile}{ A character string naming a sample file. See details. }
  \item{CDFpath}{ A character string naming a directory where the CDF files are
	located. }
  \item{RIpath}{ A character string naming a directory where the RI corrected text
	files are/will be located. }
  \item{RIfiles}{ Logical. If \code{TRUE}, the function will look for for RI files
  (RI_*) instead of CDF files (the default).}
  \item{ftype}{A character string giving the file type for RI files. Options are
	"binary" and "text".}
  \item{ignore.case}{Logical. Should pattern-matching be case-insensitive?}
  \item{\dots}{ Other options passed to \code{\link[utils]{read.delim}} function. }
}
\details{
	The sample file is a tab-delimited text file with at least two columns:
	
	\itemize{
	\item \code{CDF_FILE} - The list of baseline corrected CDF files.
  \item \code{MEASUREMENT_DAY} - The day when the sample was measured.
	}

	The column names must be exactly as those indicated, but the column
	order doesn't matter. Other columns could be included in that file.
	They won't be used by the script, but will be included in the "sample"
	R object.

	The \code{ftype} parameter sets the file format of the RI files, which
	are created by the function \code{\link{RIcorrect}}. "text" is the
	\emph{classic} text format and "binary" is a binary version, designed
	for speed, of the text format (TargetSearch >= 1.12.0). The file format
	can be identified by the file extension (".txt" for text, ".dat" for
	binary), but this is just an indicator: the file format is detected
	dinamically during file reading. Use the method \code{\link{fileFormat}}
	to set or get this parameter.
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
