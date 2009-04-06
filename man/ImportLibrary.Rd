\name{ImportLibrary}
\alias{ImportLibrary}
\title{ Library import }
\description{
    This function imports a metabolite library file that will be used to processed
    the GC-MS data.
}
\usage{
ImportLibrary(libfile, RI_dev = c(2000, 1000, 200), SelMasses = 5,
              TopMasses = 15, ExcludeMasses)
}
\arguments{
  \item{libfile}{ A character string naming a library file. See details. }
  \item{RI_dev}{ A three component vector with RI windows. }
  \item{SelMasses}{ The number of selective masses that will be used. }
  \item{TopMasses}{ The number of most intensive masses that will be taken from the spectrum,
    if no \code{TOP_MASSES} is provided.}
  \item{ExcludeMasses}{ Optional. A vector containing a list of masses that will be excluded. }
}
\details{

    The library file is a tab delimited text file with the following column names.
    \itemize{

    \item \code{Name} - The metabolite name.
    \item \code{RI} - The expected RI.
    \item \code{SEL_MASSES} - A list of selective masses separated with semicolon.
    \item \code{TOP_MASSES} - A list of the most abundant masses to be searched, separated
        with semicolons.
    \item \code{Win_k} - The RI windows, k = 1,2,3. Mass search is perfomed in three
        steps. A RI window required for each one of them.
    \item \code{SPECTRUM} - The metabolite spectrum. m/z and intensity are separated by
        spaces and colons.

    }

    The columns \code{Name} and \code{RI} are mandatory. At least one of columns \code{SEL_MASSES},
    \code{TOP_MASSES} and \code{SPECTRUM} must be given as well. By using the
    parameters \code{SelMasses} or \code{TopMasses} it is possible to set the selective
    masses or the top masses from the spectra. The parameter \code{ExcludeMasses} is
    used only when masses are obtained from the spectra.
    The parameter \code{RI_dev} can be used to set the RI windows.
    Note that in this case, all metabolites would have the same RI windows.

}    

\value{
  A \code{tsLib} object.
}
\examples{
# get the reference library file
cdfpath <- file.path(.find.package("TargetSearchData"), "gc-ms-data")
lib.file  <- file.path(cdfpath, "library.txt")

# Import the reference library
refLibrary <- ImportLibrary(lib.file) 

# set new names for the first 3 metabolites
libName(refLibrary)[1:3] <- c("Metab01", "Metab02", "Metab03")

# change the retention time deviations of Metabolite 3
RIdev(refLibrary)[3,] <- c(3000,1500,150)

}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{ImportSamples}}, \code{\linkS4class{tsLib}} }

