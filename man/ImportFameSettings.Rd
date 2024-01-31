\name{ImportFameSettings}
\alias{ImportFameSettings}
\title{ Retention time markers settings }
\description{
  This function imports a list of retention standard markers for RI correction.
}
\usage{
ImportFameSettings(file, mass=NULL, standard=NULL, colnames=FALSE, ...)
}
\arguments{
  \item{file}{ A character string naming a tab-delimited text file with standard
    markers OR a matrix-like object. }
  \item{mass}{ Optional. The m/z standard marker(s). Either a scalar, which
    applies to every marker, or a vector. }
  \item{standard}{ Optional. The RI values of the standards. A numeric vector. }
  \item{colnames}{ Logical flag. If \code{TRUE}, the column names will be used
    to look up for the settings, otherwise the column positions will be used.
    See details below. }
  \item{\dots}{ Options passed to \code{\link[utils:read.table]{read.delim}},
    which is used internally to read the settings file. }
}
\details{
  The standard marker file is a tab-delimited text file with 2, 3 or 4
  columns. Additional columns are ignored. If the parameter \code{colnames}
  is \code{FALSE}, then the column names are ignored and they must be in the
  following specified below. If \code{TRUE}, use the column names below (case
  is ignored). Note that some columns are optional as they can be specified by
  the options \code{mass} and \code{standard}.

  \itemize{
    \item \code{LowerLimit} - The Retention time lower limit in seconds
      (required).
    \item \code{UpperLimit} - The Retention time upper limit in seconds
      (required).
    \item \code{RIstandard} - The RI value of that standard. This is optional
      (option \code{standard}).
    \item \code{mass} - The m/z standard marker. This is optional
      (option \code{mass}).
  }

  Instead of passing a file via the \code{file} option, one can pass a
  \code{matrix} or a \code{data.frame}, The same aforementioned conditions
  apply in terms of column names and number.

  In general, the \code{mass} and the \code{standard} should be to \code{NULL}.
  If not, then they take precedence over the values defined by the \code{file}
  parameter.

  If no arguments are given, a default object will be returned.
}
\value{
    A \code{tsRim} object.
}
\note{
  If the parameter \code{colnames} is set to \code{TRUE}, then care needs to be
  taken regarding the column names: they must be named exactly as described in
  the \emph{details} section, otherwise they can be silently ignored and a
  default value might be used instead. This has been the default behavior.
}

\examples{
require(TargetSearchData)
# get the RI marker definition file
rim.file <- tsd_file_path("rimLimits.txt")

# set the mass marker to 87
mass <- 87

# load the definition
rimLimits <- ImportFameSettings(rim.file, mass = mass)

# sometimes you need to change the limits of a particular standard
rimLimits(rimLimits)[2,] <- c(410, 450)

# to change the mass value
rimMass(rimLimits) <- 85

# alternatively, you can pass a two column matrix as limits and pass the
# rest as parameters separately (toy example).
limits <- cbind(c(10, 20, 30), c(15, 25, 35))
mass <- 100
standard  <- c(100, 200, 300)
rimLimits <- ImportFameSettings(limits, mass, standard)

# or combine all into a matrix (note the column order!)
def <- cbind(limits, standard, mass)
rimLimits <- ImportFameSettings(def)

}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{RIcorrect}}, \code{\linkS4class{tsRim}},
  \code{\link[utils:read.table]{read.delim}}}
