# common functions

is_nullOrNA <- function(x)
{
    if(is.null(x)) return(TRUE)
    if(any(is.na(x))) return(TRUE)
    FALSE
}

#' like `dirname`, but checks if all dirs are the same
#'
#' @param x a vector of file paths
#' @return the dirname for each member of x. If all names are equal, return
#'         the unique value
.dirname <- function(x)
{
    ret <- dirname(x)
    unq <- unique(ret)
    if(length(unq) == 1)
        return(unq)
    return(ret)
}

#' update file path
#'
#' Updates the path of new files. If the length of new is 1, then apply the
#' the same path to all files. if NULL or NA, return old paths
#' Assumes that both length are equal or the paths can be recycled. Raises
#' and error if it is not possible
#'
#' @param old a vector of file paths
#' @param new a vector of new paths
#' @return the updated file paths
.setpath <- function(old, new)
{
    if(is_nullOrNA(new))
        return(old)
    if(length(old) %% length(new) != 0)
        stop("lengths must be multiples of each other")
    if(length(new) > length(old)) {
        warning('length of `new` is longer than `old`. Extra elements will be ignored')
        new <- new[seq(length(old))]
    }
    new <- rep(new, length(old)/length(new))
    file.path(new, basename(old))
}

#' parse file names
#'
#' A simple parser for file names
.parse_file_names <- function(files, exts)
{
    e <- paste(exts, collapse="|")
    m <- regexpr(sprintf("^(.+)\\.(%s)$", e), files, perl=TRUE, ignore.case=TRUE)
    f <- function(x, m, len) if(m == -1) NA_character_ else substring(x, m, m + len - 1)
    cstart <- attr(m, 'capture.start')
    clen <- attr(m, 'capture.length')
    name <- mapply(f, files, cstart[, 1], clen[, 1], USE.NAMES=FALSE)
    extension <- mapply(f, files, cstart[, 2], clen[, 2], USE.NAMES=FALSE)
    name[ is.na(name) ] <- files[ is.na(name) ]
    cbind(name, extension)
}

#' trim file extension
#'
#' remove known file extensions from file path. if none found, then does
#' not do anything.
#'
#' @param files file names
#' @param exts  file extensions to trim. the extensions are converted to
#'              to regular expressions
#' @return files with extension trimmed.
.trim_file_ext <- function(files, exts)
{
    .parse_file_names(files, exts)[,1]
}

#' add file extension
#'
#' add file extensions to extension trimmed file names and checks if
#' the file exists. In case of ambiguity, takes the first
#' @param files a list of file paths
#' @param exts a list of file extensions
#' @return take the first path/extension combo that exists, if none exists
#'         take the first.
.add_file_ext <- function(files, exts)
{
    f <- lapply(exts, function(e) paste(files, e, sep="."))
    f <- lapply(seq(files), function(i) sapply(f, getElement, i))
    sapply(f, function(z) {
               y <- file.exists(z)
               if(all(!y)) z[1] else z[y][1] })
}

#' make RI files
#'
#' Takes a list of CDF files and create a list of RI files by appending
#' RI_ to the file. It takes care of paths.
#' @param files a list of file paths
#' @param exts a list of file extensions
#' @return take the first path/extension combo that exists, if none exists
.make_RI_files <- function(files, type, exts)
{
    bs <- basename(files)
    dr <- dirname(files)
    mapply(function(d, b) {
              z <- sprintf("RI_%s.%s", .trim_file_ext(b, exts), type)
              if(d == ".") z else file.path(d, z)
        }, dr, bs, USE.NAMES=FALSE)
}

`.find_uniq_files` <- function(files, exts)
{
    if(length(files) <= 1)
        return(files)
    p <- .parse_file_names(files, exts)
    p <- cbind(files, p)
    k <- match(tolower(p[, 'extension']), tolower(exts))
    k[is.na(k)] <- 1
    p <- p[order(k, p[, 'name']),,drop=FALSE]
    p[ !duplicated(p[, 'name']), 'files']
}

# vim: set ts=4 sw=4 et:
