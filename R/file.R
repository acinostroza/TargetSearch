# miscellaneous file functions

`get.columns.name.default` <- function()
{
	c(spectrum="SPECTRUM",retIndex="RETENTION_TIME_INDEX",retTime="RETENTION_TIME")
}

# function to get the default RI text columns
`get.columns.name` <- function(cols)
{
	if(missing(cols) || is.null(cols) || any(is.na(cols)))
		cols <- getOption('TS_RI_columns', get.columns.name.default())
	assert_that(length(cols) == 3, msg='Option `TS_RI_columns` must have length=3')
	if(is.character(cols) || is.integer(cols))
		return(cols)
	if(is.numeric(cols))
		return(as.integer(cols))
	stop('Option `TS_RI_columns` must be integer or character')
}

# return file header as string to write TXT files
`get.file.header` <- function() {
	cols <- get.columns.name.default()
	paste(cols['retTime'], cols['spectrum'], cols['retIndex'], sep="\t")
}

# Convert RI files from text to binary format and viceversa.
`.convert_ri_file` <-
function(in.file, ...)
{
	if((ret <- .Call(c_convert_ri_file, in.file, ...)) != 0)
		stop("An error occurred converting file '", in.file, "'")
	ret
}

`bin2text` <-
function(in.files, out.files=NULL)
{
	if(is.null(out.files))
		out.files <- paste(sub("\\.\\w+$", "", in.files), ".txt", sep="")

	assert_that(length(in.files) == length(out.files))
	header <- get.file.header()
	ret <- vapply(seq(length(in.files)), function(k)
				  .convert_ri_file(in.files[k], out.files[k], 0L, NULL, header), 0L)
	invisible(out.files)
}

`text2bin` <-
function(in.files, out.files=NULL, columns=NULL)
{
	if(is.null(out.files))
		out.files <- paste(sub("\\.\\w+$", "", in.files), ".dat", sep="")

	assert_that(length(in.files) == length(out.files))
	columns <- get.columns.name(columns)
	ret <- vapply(seq(length(in.files)), function(k)
				.convert_ri_file(in.files[k], out.files[k], 1L, columns, NULL), 0L)
	invisible(out.files)
}

# Read peak list in binary format using the R (not C) interface.
# Args:
#  - f: RI binary file name
# Value:
#   a list with unformatted raw data with components
#    - retIndex: retention time index (length: n)
#    - retTime:  retention time (length: n)
#    - N:        number of mass/intensity pairs per scan (length: n)
#    - Values:   a vector of masses and intensities (length: sum(n)*2)
#   where 'n' is the number of scans.

`readRIBin` <-
function(f) {
	z <- file(f, "rb")
	sig <- readBin(z, what="int", n=2, endian="little")
	stopifnot(all(sig == c(169603882, 84919)))
	n   <- readBin(z, what="int", n=1, endian="little")
	m   <- readBin(z, what="int", n=1, endian="little")
	RI  <- readBin(z, what="numeric", n=n, endian="little")
	RT  <- readBin(z, what="numeric", n=n, endian="little")
	N   <- readBin(z, what="int", n=n, endian="little")
	tmp <- readBin(z, what="int", n=2*m, endian="little")
	close(z)
	list(retIndex=RI,retTime=RT,N=N,Values=tmp)
}

# Write peak list in binary format using the R (not C) interface.
# Args:
#  - x: a list with components retIndex, retTime, N, Values (see readRIBin)
#  - f: a file.

`writeRIBin` <-
function(x, f) {
	z <- file(f, "wb")
	sig <- c(169603882, 84919)
	writeBin(as.integer(sig), z, endian="little")
	writeBin(as.integer(length(x$N)), z, endian="little")
	writeBin(as.integer(sum(x$N)), z, endian="little")
	writeBin(x$retIndex, z, endian="little")
	writeBin(x$retTime, z, endian="little")
	writeBin(as.integer(x$N), z, endian="little")
	writeBin(as.integer(x$Values), z, endian="little")
	close(z)
	invisible()
}

# Read retention times from a binary file
# f = RI binary file
`readRetTimes` <-
function(f)
{
	z <- file(f, "rb")
	sig <- readBin(z, what="int", n=2, endian="little")
	if(!all(sig == c(169603882, 84919)))
		stop(sprintf("Incorrect binary format of file %s", f))
	n   <- readBin(z, what="int", n=2, endian="little")
	RI  <- readBin(z, what="numeric", n=n[1], endian="little")
	RT  <- readBin(z, what="numeric", n=n[1], endian="little")
	close(z)
	cbind(retIndex=RI,retTime=RT)
}

#' Wrapper to C function guess_file_type
#'
#' This is calls the C function `guess_file_type` which in turns calls
#' `file_type`. This function returns 1 for text files, 0 for binary
#' files, and -1 for errors (for example if file does not exist).
#'
#' @param filename (string) path to file
#' @return integer. The (guessed) type of the file
`file_type` <-
function(filename)
{
	assert_that(is.string(filename))
	.Call(c_guess_file_type, filename)
}

# vim: set ts=4 sw=4 noet:
