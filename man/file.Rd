\name{text2bin}
\alias{text2bin}
\alias{bin2text}
\title{Convert RI files from text to binary format and viceversa.}
\description{
	This function converts a list of RI files (peak list files)
	in text (binary) format to binary (text) format.
}
\usage{
text2bin(in.files, out.files=NULL, columns=NULL)

bin2text(in.files, out.files=NULL)
}
\arguments{
	\item{in.files}{A character string naming the input files. }

	\item{out.files}{A character string naming the output files. If
		\code{NULL}, the input file extensions will be changed
		accordingly ("txt" to "dat" and viceversa).}

	\item{columns}{A numeric vector with the positions of the columns
		\code{SPECTRUM}, \code{RETENTION_TIME_INDEX}, and
		\code{RETENTION_TIME} or a character vector with the header
		names of those columns. This parameter is required only
		for \code{text2bin}.}
}
\details{
	The format of the input files is detected dinamically and error
	will be issued if it is incorrect.

	Note that the respective \code{\linkS4class{tsSample}} object
	may need to be updated by using the method \code{\link{fileFormat}}.
}
\value{
	A character vector of the created files.
}
\author{Alvaro Cuadros-Inostroza}
\seealso{\code{\link{ImportSamples}}, \code{\linkS4class{tsSample}},
	\code{\link{RIcorrect}}}
