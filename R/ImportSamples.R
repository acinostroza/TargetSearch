#' check if files exist and warn if not
#'
#' @param path the file path
#' @param files the files to check
#' @param quiet logical, be quiet
#' @return logical. whether the files exist

`.check.file.exists` <-
function(path, files, quiet=FALSE)
{
    z <- file.exists(file.path(path, files))
    if(quiet | all(z))
        return(z)
    cat(sprintf("=> Notice: Some files were not found in path `%s`", path), sep="\n")
    cat(sprintf("   - %s", files[!z]), sep="\n")
    warning("some files are missing")
    invisible(z)
}

#' detect file extension if not present
#'
#' tries extensions cdf or nc4 if missing. Returns the files if they exist
#' or just assumes cdf extension
#' @param path the file path
#' @param files the files to check
#' @return a character vector with file names
`.detect.file.ext` <-
function(path, files)
{
    ret <- regexpr("\\.(cdf|nc4)", files)
    if(all(ret != -1)) # files have extension
        return(files)
    a <- file.path(path, sprintf("%s.cdf", files))
    b <- file.path(path, sprintf("%s.nc4", files))
    ret <- mapply(function(a, b, f) {
                      if(file.exists(a)) return(a)
                      if(file.exists(b)) return(b)
                      f}, a, b, files)
    basename(unname(ret))
}

`ImportSamples` <-
function(sampfile, CDFpath = ".", RIpath = ".", ftype=c("binary", "text"), ...)
{
    ftype <- pmatch(ftype[1], c("binary", "text"))
    stopifnot(!is.na(ftype))

    Samples <- if(is.data.frame(sampfile)) sampfile else
                   read.delim(sampfile, ...)

    # look for column 'CDF_FILE' and 'MEASUREMENT_DAY'
    cdf <- 'CDF_FILE'
    day <- 'MEASUREMENT_DAY'

    if(!cdf %in% colnames(Samples)) {
        # greps for 'CDF'
        k <- grep('cdf', colnames(Samples), ignore.case=TRUE)
        if(length(k) == 0)
            stop("Column 'CDF_FILE' is missing!!")
        cdf <- colnames(Samples)[k][1]
        message("Note: Using '", cdf, "' as CDF column")
    }

    Samples[[ cdf ]] <- as.character(Samples[, cdf])

    if(!day %in% colnames(Samples)) {
        # greps for 'DAY'
        k <- grep('day', colnames(Samples), ignore.case=TRUE)
        if(length(k) == 0) {
            Samples[[ day ]] <- "1"
            warning("Column 'MEASUREMENT_DAY' not found. Using default setting.")
        } else {
            day <- colnames(Samples)[k][1]
        }
    }
    CDFfiles <- .detect.file.ext(CDFpath, Samples[[ cdf ]])
    chk     <- .check.file.exists(CDFpath, CDFfiles)
    ext     <- c("dat","txt")[ftype]
    RIfiles <- paste0("RI_", CDFfiles, ".", ext)
    RIfiles <- sub("\\.(cdf|nc4)\\.", ".", RIfiles, ignore.case = T)
    Names   <- gsub(".(cdf|nc4)$", "", Samples[[ cdf ]], ignore.case = T)
    new("tsSample", Names = Names, CDFfiles = CDFfiles, days = as.character(Samples[[ day ]]),
        RIfiles = RIfiles, CDFpath = CDFpath, RIpath = RIpath, data = Samples)
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
function(CDFpath=".", RIfiles=FALSE, ignore.case=TRUE, ftype=c("binary", "text")) {

    if(RIfiles == TRUE) {
        cdffiles <- dir(path=CDFpath, pattern='^RI_', ignore.case=ignore.case)
        ext <- unique(substring(tolower(cdffiles), nchar(cdffiles)-2))
        if(!all(ext %in% c("dat","txt")))
            stop("RI file extension error. It must be either 'dat' or 'txt'")

        if(length(ext) == 1) {
            cdffiles <- gsub("RI_", "", gsub(ext, "cdf", cdffiles, ignore.case=ignore.case))
            ftype <- if(ext=="dat") "binary" else "text"
        } else {
            ft <- pmatch(ftype, c("binary", "text"))[1]
            if(is.null(ft) | is.na(ft) | length(ft)==0)
                stop(paste("Binary and text RI files found. Please select one file type",
                           "with the parameter 'ftype'."))
            pattern <- c("\\.dat$","\\.txt$")[ft]
            cdffiles <- gsub("RI_", "", gsub(pattern, ".cdf", cdffiles, ignore.case=ignore.case))
            cdffiles <- grep("cdf", cdffiles, value=TRUE, ignore.case=ignore.case)
        }
    } else {
        cdffiles <- dir(path=CDFpath, pattern='\\.cdf$', ignore.case=ignore.case)
    }
    if(length(cdffiles) == 0)
        stop('Error: No CDF files were find in the directory')

    new('tsSample', CDFfiles=cdffiles,RIpath=CDFpath,CDFpath=CDFpath, ftype=ftype)
}
# vim: set ts=4 sw=4 et:
