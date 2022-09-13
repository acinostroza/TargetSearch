\name{ImportSamples}
\alias{ImportSamples}
\alias{ImportSamplesFromDir}
\title{ Sample definitions }
\description{
	These functions import sample information from either a tab delimited file
	or a file system path. The sample information contain sample names, CDF
	files and their measurement day, also known as, measurement run or measurement
	batch.
}
\usage{
ImportSamples(sampfile, CDFpath, RIpath, ftype=c("binary", "text"), ...)

ImportSamplesFromDir(CDFpath=".", RIfiles=FALSE, ftype=c("binary", "text"),
		verbose=FALSE, ...)
}
\arguments{
  \item{sampfile}{ A character string naming a sample file or a \code{data.frame}. See details. }
  \item{CDFpath}{ A character string naming a directory where the CDF files are
	located (by default is the current directory) }
  \item{RIpath}{ A character string naming a directory where the RI corrected text
	files are/will be located. }
  \item{RIfiles}{ Logical. If \code{TRUE}, the function will look for for RI files
  (RI_*) instead of CDF files (the default).}
  \item{ftype}{A character string giving the file type for RI files. Options are
	"binary" and "text".}
  \item{verbose}{Logical. Show more output if \code{TRUE}.}
  \item{\dots}{Extra options whose meaning depend on the function: for \code{ImportSamples},
  the options will be passed to \code{\link[utils:read.table]{read.delim}}, while
  for \code{ImportSamplesFromDir}, these are passed to \code{\link[base:list.files]{dir}}}.
}
\details{
	The sample file is a tab-delimited text file or, alternatively, a
	\code{data.frame} like the one would get by calling
	\code{\link[utils:read.table]{read.delim}}. The following columns
	are expected, though only the first one below is mandatory.

	\itemize{
	\item \code{CDF_FILE} - The list of baseline corrected CDF files (mandatory).
	\item \code{MEASUREMENT_DAY} - The day or batch when the sample was measured (optional).
	\item \code{SAMPLE_NAME} - A unique sample identifier (optional).
	}

	The function first looks for columns matching those names. If they are not
	found, then it looks for columns with the substrings 'cdf', 'day' and 'name' respectively.
	If the 'cdf' column is not found, the function raises an error, but if 'day'
	is not found, then it assumes the measurement day is the same for all.
	Moreover, if the column for sample names/idenfiers is missing, then the function uses
	the 'cdf' file names as identifiers.

	During the column matching, the first match is taken in case of ambiguity.
	The column order does not matter. Other columns could be included in that file.
	They won't be used by the script, but will be included in the "sample"
	R object.

	The files given in the 'cdf' column may optionally include a relative path,
	such as \code{mypath/myfile.cdf}, which will be relative to the working
	directory or relative to the argument \code{CDFpath}. Optionally, the
	file extension can be omitted because \sQuote{TargetSearch} adds the \code{.cdf}
	extension automatically if missing. Note: it may fail in case sensitive
	file systems as a lower case extension is appended.

	If the parameter \code{CDFpath} is missing, then the current directory
	is used. If \code{RIpath} is missing, the value of \code{CDFpath} will
	be used.

	The \code{ftype} parameter sets the file format of the RI files, which
	are created by the function \code{\link{RIcorrect}}. "text" is the
	\emph{classic} text format and "binary" is a binary version, designed
	for speed, of the text format (TargetSearch >= 1.12.0). The file format
	can be identified by the file extension (".txt" for text, ".dat" for
	binary), but this is just an indicator: the file format is detected
	dinamically during file reading. Use the method \code{\link{fileFormat}}
	to set or get this parameter.

	The function \code{ImportSamplesFromDir} scans for cdf files in the current
	or given directory. The search can be made recursive by passing the
	parameter \code{recursive=TRUE} (actually passed to \code{\link[base:list.files]{dir}}).

	The parameter \code{verbose} will print a list of detected files
	is set to \code{TRUE}. This can be used for debugging.
}
\value{
	A \code{tsSample} object.
}
\examples{
# get the sample definition definition file
require(TargetSearchData)
cdfpath <- tsd_data_path()
sample.file  <- tsd_file_path("samples.txt")

# set a path where the RI files will be created
RIpath <- "."

# import samples
sampleDescription <- ImportSamples(sample.file, CDFpath = cdfpath, RIpath = RIpath)

# change the sample names. It is required that the names must be unique.
sampleNames(sampleDescription) <- paste("Sample", 1:length(sampleDescription), sep = "_")

# change the file paths (relative to the working path)
CDFpath(sampleDescription) <- "my_cdfs"
RIpath(sampleDescription)  <- "my_RIs"

}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{ImportLibrary}}, \code{\linkS4class{tsSample}} }
