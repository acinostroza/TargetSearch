# function to import a library. a general function for
# MSP or tab-delimited format.
# Arguments:
#  - x: a library file or a data.frame
#  - type: the file type: options are "auto" (autodetection), "tab" and "msp"

ImportLibrary <- function(x, type = c("auto", "tab", "msp"), ...) {

    # detects file type: either tab-delimited or NIST MSP format
    if(is.data.frame(x)) {
        return(ImportLibrary.tab(libdata=x, ...))
    }
    type <- match.arg(type)

    if(type == "auto") {
        line <- readLines(x, n = 1)
        if(grepl("\t", line)) {
            type <- "tab"
        } else if(grepl("^Name:", line)) {
            type <- "msp"
        } else {
            stop("Error: library format not recognized")
        }
    }

    if(type == "tab")
        ImportLibrary.tab(x, ...)
    else
        ImportLibrary.msp(x, ...)
}

.optfun <- function(file, opt)
{
	stopifnot(is.list(opt) | is.null(opt))
	# default options passed to read.table
	def <- list(file=file, header=TRUE, sep="\t", quote="", dec=".",
			fill=TRUE, comment.char="")
	for(n in names(opt))
		def[[n]] <- opt[[n]]
	def
}

# function to import a library (from a tab delimited file):
# Arguments
#  - libfile: The library file
#  - fields: meaningless parameter. Included to be compatible with ImportLibrary.txt
#  - RI_dev: RI deviation
#  - SelMasses: How many masses should be used as selective if no SelMasses are found.
#               Masses are taken from the TopMasses
#  - TopMasses: How many masses from the spectra should be uses as top Masses.
#  - ExcludeMasses: Don't take this masses as TopMasses automatically.
#  - libdata: A data frame containing the data from a library file (optional)
#  - file.opt: A list. Further options passed to read.delim

ImportLibrary.tab <- function(libfile, fields = NULL, RI_dev = c(2000,1000,200),
    SelMasses = 5, TopMasses = 15, ExcludeMasses = NULL, libdata, file.opt=NULL)
{
	if(missing(libfile)) {
		if(missing(libdata)) {
			stop("argument \"libfile\" and \"libdata\" are missing, with no default")
		} else if(is.data.frame(libdata)) {
			M <- libdata
		} else {
			stop("argument \"libdata\" should be a 'data.frame'")
		}
	} else {
		file.opt <- .optfun(libfile, file.opt)
		M <- do.call("read.table", file.opt)
	}
	M <- .check.data.frame(M)

    if(is.null(M[[ 'Name' ]]))
        stop("Column 'Name' is missing!!")

    if(is.null(M[[ 'RI' ]]))
        stop("Column 'RI' is missing!!")

    cols <- c('SEL_MASS', 'TOP_MASS', 'SPECTRUM')
    if(all(!cols %in% colnames(M)))
        stop("One of the following columns must be given: ",
             paste(cols, collapse=", "), ".")

    # force numeric
    M$RI <- as.numeric(M$RI)
    if(any(is.na(M$RI)))
        stop("Missing values in RI definition found. Please check input data.")

    if(is.null(M[[ 'SEL_MASS' ]]))
        M$SEL_MASS <- ""
    if(is.null(M[[ 'TOP_MASS' ]]))
        M$TOP_MASS <- ""

    q_mass <- M[[ 'QUANT_MASS' ]]

    v <- Vectorize(function(x) {
                       x <- as.numeric(x)
                       x[!is.na(x)]
    }, SIMPLIFY=FALSE, USE.NAMES=FALSE)

    sel_mass  <- v(strsplit(as.character(M$SEL_MASS), "[;\\|:, ]+"))
    top_mass  <- v(strsplit(as.character(M$TOP_MASS), "[;\\|:, ]+"))

    if(!is.null(spectra <- M[[ 'SPECTRUM']])) {
        spectra <- Spectra(as.character(M$SPECTRUM))
        temp <- lapply(spectra, Top.Masses, TopMasses, ExcludeMasses)
        top_mass <- mapply(function(a, b) if(length(a) > 0) a else b, top_mass,
                       temp, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    }

    # consolidate selective and top masses
    f <- function(a, b) if(length(a) > 0) a else b
    sel_mass <- mapply(f, sel_mass, top_mass, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    top_mass <- mapply(union, sel_mass, top_mass, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    if(any(k <- (vapply(sel_mass, length, 0) == 0)))
        stop("Invalid entry detected (# ", which(k)[1], ")")

    w1 <- if(is.null(M$Win_1)) rep(RI_dev[1], length(M$RI)) else M$Win_1
    w2 <- if(is.null(M$Win_2)) rep(RI_dev[2], length(M$RI)) else M$Win_2
    w3 <- if(is.null(M$Win_3)) rep(RI_dev[3], length(M$RI)) else M$Win_3

    new("tsLib", Name = as.character(M$Name), RI = M$RI, medRI = M$RI,
        RIdev = cbind(w1,w2,w3), selMass = sel_mass, topMass = top_mass,
        quantMass=q_mass, spectra = spectra, libData = M)
}

Top.Masses <- function(sp, TopMasses, ExcludeMasses = NULL) {
    if(is.null(sp) || length(sp) == 0)
        return(numeric(0))
	m <- sp[order(sp[,2], decreasing = TRUE),1]
	if(!missing(ExcludeMasses) & !is.null(ExcludeMasses)) {
			if(!is.numeric(ExcludeMasses))
				stop("'ExcludeMasses' must be a numeric vector")
			m <- m[!m %in% ExcludeMasses]
	}

	m[1:min(length(m), TopMasses)]
}

Spectra <- function(spectrum)
{
    # check spectra format
    assert_that(is.character(spectrum))
    if(any(k <- grepl("[^0-9 :;]", spectrum)))
        stop("Error: Invalid characters in entry #", which(k)[1], ": '", spectrum[k][1], "'")

    spectrum[is.na(spectrum)] <- ''
    spectrum <- gsub("[:;]", " ", spectrum)
    y <- lapply(strsplit(spectrum, " +"), function(z) as.numeric(z[z != ""]))
    makemat <- function(x) {
        n <- length(x)
        n <- n - (n %% 2)
        if(n > 0) matrix(x[seq(n)], ncol=2, byrow=TRUE) else numeric(0)
    }
    z <- sapply(y, makemat, simplify=FALSE)
    return(z)
}

# function to read a msp file. Returns a list where every component represent
# a metabolite/analyte entry in the MSP. Every entry is a list with
# following components:
#   Name: The metabolite/analite name:
#   formula: the chemical formula (NA if no "Formula:" field is found)
#   mw: Molecular weight (NA if no "MW:" field is found)
#   n_peaks: Number of peaks.
#   spectrum: A two column matrix with m/z and intensities.
#   Synon: A character vector with the metabolite synonyms (Synon: fields)

read.msp <- function(msp.file) {
    lines   <- readLines(msp.file)
    msp.id.end   <- which(lines == "")
    if(length(msp.id.end) == 1) {
        msp.id.end   <- msp.id.end - 1
        msp.id.start <- 1
    } else if (length(msp.id.end) > 1) {
        msp.id.start <- c(0, msp.id.end[1:(length(msp.id.end)-1)]) + 1
        msp.id.end   <- msp.id.end - 1
    } else {
        stop("Empty line separators not found. Make sure this is a proper MSP file.")
    }

    msp <- lapply( 1:length(msp.id.start), function(x) {
        tmp  <- lines[msp.id.start[x]:msp.id.end[x]]

        if(grep("^Name: ?", tmp, ignore.case = TRUE) != 1)
            stop("Invalid MSP format. Expecting: 'Name: <metab name>', found:", tmp[1])

        name <- sub("^Name: ?", "", tmp[1], ignore.case = TRUE)
        mw       <- sub("^MW: ?", "", tmp[ grep("^MW: ", tmp) ])
        if(length(mw) == 0) mw <- NA
        formula  <- sub("^Formula: ?", "", tmp[ grep("^Formula: ?", tmp) ])
        if(length(formula) == 0) formula <- NA
        n_peaks.id    <- grep("^Num Peaks: ?", tmp)
        if(length(n_peaks.id) == 0) stop("Error in MSP file: Metab = '",name, "'\nField 'Num peaks: ' doesn't exist")
        n_peaks  <- as.numeric(sub("^Num Peaks: ?", "", tmp[ n_peaks.id ]))
        spectrum <- as.numeric(unlist(strsplit( paste(tmp[(n_peaks.id+1):length(tmp)], collapse = " "), "[: ;\t]+", perl = TRUE)))
        spectrum <- spectrum[!is.na(spectrum)]
        if(length(spectrum) != 2*n_peaks)
            stop("Invalid number of peaks in MSP file: ", name)
        spectrum <- matrix(spectrum, ncol = 2, byrow = TRUE, dimnames = list(1:n_peaks, c("mz", "intensity")))
        syn <- grep("Synon", tmp, value = TRUE)
        list( Name = name, mw = mw, formula = formula, n_peaks = n_peaks, spectrum = spectrum, Synon = syn)
    })
    rm(lines)
    msp
}

# import a library from a MSP file. It has the same arguments as ImportLibrary.tab
# except fields.

ImportLibrary.msp <- function(libfile, fields = NULL, RI_dev = c(2000,1000,200),
   SelMasses = 5, TopMasses = 15, ExcludeMasses = NULL) {

    if(length(RI_dev) != 3)
        stop("'RI_dev' must be a vector of length 3")

    msp <- read.msp(libfile)
    topMass <- lapply(msp, function(x) Top.Masses(x$spectrum, TopMasses, ExcludeMasses))

    if(is.null(fields)) {
        ri <- rep(0, length(msp))
        selMass <- lapply(topMass, function(x) x[1:SelMasses])
    } else if(is.list(fields) == FALSE) {
        stop("'fields' must be a list")
    } else {
        n <- length(fields)
        if(n >= 1) {
            ri <- sapply(msp, function(x) {
                re <- paste("^Synon: ?", fields[[1]][1], " *(.+)$", sep = "")
                id <- grep(re, x$Synon)
                if(length(id) > 0) {
                    as.numeric( sub(re, "\\1", x$Synon[id[1]] ) )
                } else 0
            })

        } else stop("'fields' must have at list one member")

        if(n >= 2) {
            selMass <- lapply(msp, function(x) {
                re <- paste("^Synon: ?", fields[[2]][1], " *(.+)$", sep = "")
                id <- grep(re, x$Synon)
                if(length(id) > 0) {
                    strsplit(sub(re, "\\1", x$Synon[id[1]]), "[;\\|:, ]+" )[[1]]
                } else "NA"
            })
        } else
            selMass <- as.list(rep("NA", length(msp)))

        selMass <- mapply(function(x,y) if(any(x == "NA")) y[1:SelMasses] else x, selMass, topMass, SIMPLIFY = FALSE)
        selMass <- lapply(selMass, as.numeric)
    }

    libData <- data.frame(Name = sapply(msp, function(x) x$Name), libRI = ri,
        MW = as.numeric(sapply(msp, function(x) x$mw)), Formula = sapply(msp, function(x) x$formula),
        stringsAsFactors = FALSE)

    new("tsLib", Name = libData$Name, RI = ri, medRI = ri,
        RIdev = matrix(rep(RI_dev, length(msp)), ncol = 3, byrow = TRUE),
        selMass = selMass, topMass = topMass,
        spectra = lapply(msp, function(x) x$spectrum), libData = libData)
}

# function to remove unexpected quotation marks (usually generated by
# external software such as Excel). Also converts all factors to character

.check.data.frame <- function(x) {
	if(!is.data.frame(x))
		return(x)
	for(i in 1:ncol(x)) {
		if(is.numeric(x[,i]))
			next
		if(is.factor(x[,i]))
			x[,i] <- as.character(x[,i])
		if(any(grepl("^['\"]", x[,i])) | any(grepl("['\"]$", x[,i]))) {
			warning("Unexpected single/double quotation characters found in column '", names(x)[i],
				"'.\n  These characters will be removed. Please check input file or data.frame.")
			x[,i] <- sub("^['\"]", "", x[,i])
			x[,i] <- sub("['\"]$", "", x[,i])
		}
	}
	return(x)
}

# vim: set ts=4 sw=4:
