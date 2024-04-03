#' check if files exist and warn if not
#'
#' @param path the file path
#' @param files the files to check
#' @param quiet logical, be quiet
#' @return logical. whether the files exist

`.check.file.exists` <-
function(files, quiet=FALSE)
{
    z <- file.exists(files)
    if(quiet | all(z))
        return(z)
    cat("Note: Some files were not found in specified path", sep="\n")
    cat(files[!z], fill=TRUE, labels=" ->")
    warning("Some files are missing!")
    invisible(z)
}

#' check file extensions
#'
#' If a file has a proper file extension, we leave it. If not, then tries
#' to find an appropriate extension. If found, then returns, if not, return
#' the file unchanged.
#' @param files the files to check
#' @return a character vector with file names
`.detect.files` <-
function(files)
{
    fun <- function(x, m) {
        if(m != -1)
            return(x)
        for(e in c('.nc4', '.cdf', '.CDF')) {
            f <- paste0(x, e)
            if(file.exists(f))
                return(f)
       }
       return(x)
    }
    m <- regexpr("\\.(cdf|nc4)$", files, perl=TRUE, ignore.case=TRUE)
    mapply(fun, files, m, USE.NAMES=FALSE)
}

`ImportSamples` <-
function(sampfile, CDFpath, RIpath, ftype=c("binary", "text"), ...)
{
    ftype <- match.arg(ftype)
    as_chr <- as.character

    Samples <- if(is.data.frame(sampfile)) sampfile else
                   read.delim(sampfile, ...)

    # look for column 'CDF_FILE' and 'MEASUREMENT_DAY' and 'SAMPLE_NAME'
    cdf <- 'CDF_FILE'
    day <- 'MEASUREMENT_DAY'
    snm <- 'SAMPLE_NAME'

    if(!cdf %in% colnames(Samples)) {
        # greps for 'CDF'
        k <- grep('cdf', colnames(Samples), ignore.case=TRUE)
        if(length(k) == 0)
            stop("Column 'CDF_FILE' is missing!!")
        cdf <- colnames(Samples)[k][1]
        message("Note: Using '", cdf, "' as CDF column")
    }

    CDFfiles <- as_chr(Samples[, cdf])

    if(!day %in% colnames(Samples)) {
        # greps for 'DAY'
        day <- grep('day', colnames(Samples), ignore.case=TRUE, value=TRUE)
        if(length(day) == 0) {
            days <- "1"
            warning("Column 'MEASUREMENT_DAY' not found. Using default setting.")
        } else {
            days <- Samples[, day[1]]
            message("Note: Using '", day[1], "' as MEASUREMENT_DAY column")
        }
    } else {
        days <- Samples[, day]
    }

    if(!snm %in% colnames(Samples)) {
        k <- FALSE
        patterns <- c('^sample.?id', '^sample.?name', 'name')
        for(p in patterns) {
            k <- grepl(p, colnames(Samples), ignore.case=TRUE)
            if(any(k))
                break
        }
        if(any(k)) {
            col <- colnames(Samples)[k][1]
            message("Note: Using '", col, "' as sample names")
            Names <- as_chr(Samples[, col])
        } else {
            Names <- .trim_file_ext(basename(CDFfiles), c('nc4','cdf'))
        }
    } else {
        Names <- as_chr(Samples[, snm])
    }

    if(any(duplicated(Names)))
        stop('Duplicated sample names detected.')

    if(!missing(CDFpath))
        CDFfiles <- file.path(CDFpath, CDFfiles)

    # warn and print missing samples.
    CDFfiles <- .detect.files(CDFfiles)
    nul <- .check.file.exists(CDFfiles)

    if(missing(RIpath))
        RIpath <- dirname(CDFfiles)

    new("tsSample", Names = Names, CDFfiles = CDFfiles, days = days,
        RIpath = RIpath, data = Samples)
}

# function to extract the /measurement day/ of a set of names.

getDays <- function(x) {
    y <- regexpr('[0-9]+', x)

     # check if one element doesn't have digits. If so, just return ones.
    if(any(y== -1)) {
        return(rep('1', length(x)))
    }

    d <- substring(x, y, y+attr(y, "match.length")-1)
    # count number of elements per day
    n <- sapply(unique(d), function(z) sum(z == d))
    if(any(n == 1)) {
        # warning TODO
        return(rep('1', length(x)))
    }
    return(d)
}

`ImportSamplesFromDir` <-
function(CDFpath=".", RIfiles=FALSE, ftype=c("binary", "text"), verbose=FALSE, ...)
{
    if(RIfiles == TRUE) {
        ftype <- match.arg(ftype)
        ext <- switch(ftype, binary='dat', text='txt')
        pattern <- sprintf('^RI_(.*)\\.%s$', ext)
        rifiles <- dir(path=CDFpath, pattern=pattern, full.names=TRUE, ...)

        if(verbose) {
            cat("Searching for RI files\n======================", sep="\n")
            cat(paste("File Path      :", CDFpath), sep="\n")
            cat(paste("File Extension :", ext), sep="\n")
            cat(paste("File Pattern   :", pattern), sep="\n")
            cat("Detected Files :", sep="\n")
            cat(rifiles, fill=TRUE, labels=" ->")
        }

        if(length(rifiles) == 0)
            stop('Error: No RI files were found in the directory')

        path    <- dirname(rifiles)
        cdffiles <- sub(pattern, "\\1.nc4", basename(rifiles), perl=TRUE, ignore.case=TRUE)
        cdffiles <- file.path(path, cdffiles)
        ret <- new('tsSample', CDFfiles=cdffiles, RIfiles=rifiles, ftype=ftype)
    } else {
        pattern <- '\\.(cdf|nc4)$'
        cdffiles <- dir(path=CDFpath, pattern=pattern, full.names=TRUE, ...)

        if(verbose) {
            cat("Searching for CDF files\n======================", sep="\n")
            cat(paste("File Path      :", CDFpath), sep="\n")
            cat(paste("File Extension : cdf, nc4"), sep="\n")
            cat(paste("File Pattern   :", pattern), sep="\n")
            cat("Detected Files :", sep="\n")
            cat(cdffiles, fill=TRUE, labels=" ->")
        }

        if(length(cdffiles) == 0)
            stop('Error: No CDF files were find in the directory')

        cdffiles <- .find_uniq_files(cdffiles, c('nc4', 'cdf'))

        ret <- new('tsSample', CDFfiles=cdffiles, ftype=ftype)
    }
    ret
}

# vim: set ts=4 sw=4 et:
