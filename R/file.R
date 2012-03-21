# miscellaneous file functions

# Convert RI files from text to binary format and viceversa.

`bin2text` <-
function(in.files, out.files=NULL)
{
	if(is.null(out.files))
		out.files <- paste(sub("\\.\\w+$", "", in.files), ".txt", sep="")

	swap <- pmatch(.Platform$endian, c("little", "big")) - 1
	header <- "RETENTION_TIME\tSPECTRUM\tRETENTION_TIME_INDEX"

	stopifnot(length(in.files) == length(out.files))
	for(i in 1:length(in.files)) {
		opt <- get.file.format.opt(in.files[i], "none")
		if(opt[1] != 1) stop("Incorrect file format")
		res <- .C("dat2txt", as.character(in.files[i]), as.character(out.files[i]),
				as.integer(swap), as.character(header), PACKAGE="TargetSearch")
	}
	invisible(out.files)
}

`text2bin` <-
function(in.files, out.files=NULL,
	columns=c("SPECTRUM","RETENTION_TIME_INDEX","RETENTION_TIME"))
{
	if(is.null(out.files))
		out.files <- paste(sub("\\.\\w+$", "", in.files), ".dat", sep="")

	swap <- pmatch(.Platform$endian, c("little", "big")) - 1
	header <- "RETENTION_TIME\tSPECTRUM\tRETENTION_TIME_INDEX"

	stopifnot(length(in.files) == length(out.files))
	for(i in 1:length(in.files)) {
		opt <- get.file.format.opt(in.files[i], columns)
		if(opt[1] != 0) stop("Incorrect file format")
		res <- .C("txt2dat", as.character(in.files[i]), as.character(out.files[i]),
				as.integer(swap), as.integer(opt[3:5]), PACKAGE="TargetSearch")
	}
	invisible(out.files)
}

# Function to guess the file type (binary or text)
# get column indices (text file only)
`get.columns` <-
function(my.file, columns)
{
	if(length(columns) != 3)
		stop("Incorrect length of 'columns' argument. Should be exactly 3.")

	if(is.character(columns)) {
		header <- scan(my.file, what = "character", nlines = 1, quiet = TRUE)
		tmp <- sapply(columns, function(x) which( header == x ))
		if(length(unlist(tmp)) != length(columns))
			stop("Column name not found. Check your RI file.")
		columns <- unlist(tmp) - 1
	}
	return(columns)
}

# Checks automatically the file format and return a numeric vector
#  - file type: 0 = TXT; 1 = DAT
#  - swap: 0, 1 in little, big endian platforms
#  - SPECTRUM column number (*)
#  - RETENTION_TIME_INDEX column number (*)
#  - RETENTION_TIME column number (*)
#    (*) numbering starts from 0
`get.file.format.opt` <-
function(my.file, columns)
{
	x <- readBin(my.file, what="int", n=2, endian="little")
	if(all(x == c(169603882,84919))) { # bin file signature
		opt <- c(1, pmatch(.Platform$endian, c("little", "big")) - 1,0,0,0)
	} else { # text format
		opt <- c(0,0, get.columns(my.file, columns))
	}
	opt
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

