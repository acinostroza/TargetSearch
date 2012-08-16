# function to import a library. a general function for
# MSP or tab-delimited format.
# Arguments:
#  - libfile: The library file
#  - type: the file type: options are "auto" (autodetection), "tab" and "msp"

ImportLibrary <- function(libfile, type = "auto", ...) {

    # detects file type: either tab-delimited or NIST MSP format
    type <- pmatch(type[1], c("auto", "tab", "msp"))
    if(is.na(type))
        stop("Invalided 'type' option.")
        
    if(type == 1) {
        line <- readLines(libfile, n = 1)
        if(grepl("\t", line)) {
            type <- 2
        } else if(grepl("^Name:", line)) {
            type <- 3
        } else {
            stop("Error: library format not recognized")
        }
    }

    if(type == 2)
        ImportLibrary.tab(libfile, ...)
    else
        ImportLibrary.msp(libfile, ...)
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

ImportLibrary.tab <- function(libfile, fields = NULL, RI_dev = c(2000,1000,200),
	SelMasses = 5, TopMasses = 15, ExcludeMasses = NULL, libdata) {

	op <- options(stringsAsFactors=FALSE)

	if(missing(libfile)) {
		if(missing(libdata)) {
			stop("argument \"libfile\" and \"libdata\" are missing, with no default")
		} else if(class(libdata) == 'data.frame'){
			M <- libdata
		} else {
			stop("argument \"libdata\" should be a 'data.frame'")
		}
	} else {
		M <- read.delim(libfile, as.is = T, quote="")
	}
	M <- .check.data.frame(M)

    if(is.null(M$Name))     stop("Column 'Name' is missing!!")
    if(is.null(M$RI)) {
        stop("Column 'RI' is missing!!")
    } else { # force numeric
        M$RI <- as.numeric(M$RI)
        if(any(is.na(M$RI))) {
            stop("Missing values in RI definition found. Please check input data.")
        }
    }

	has.spectra <- TRUE
	spectra     <- list()
	
	if(is.null(M$SPECTRUM)) has.spectra <- FALSE
	if(is.null(M$SEL_MASS)) M$SEL_MASS <- NA
	if(is.null(M$TOP_MASS)) M$TOP_MASS <- NA 

	sel.mass  <- strsplit(as.character(M$SEL_MASS), "[;\\|:, ]+")
	top.mass  <- strsplit(as.character(M$TOP_MASS), "[;\\|:, ]+")

	if(has.spectra)	spectra <- Spectra(as.character(M$SPECTRUM))

	selMass <- list()
	topMass <- list()

	for(i in 1:length(sel.mass)) {
		sm <- as.numeric(sel.mass[[i]])
		tm <- as.numeric(top.mass[[i]])

		if(all(is.na(tm))) {
			if(has.spectra) {
				tm <- Top.Masses(spectra[[i]], TopMasses, ExcludeMasses)
			} else {
				if(all(is.na(sm))) stop ("Error importing library: line ", i+1, ".\nNo mass was found. Please check Library file\n.")
				tm <- sm
			}
		}

		tm <- tm[!is.na(tm)]

		if(all(is.na(sm))) {
			sm <- tm[1:min(SelMasses,length(tm))]
		} else {
			sm <- sm[!is.na(sm)]
		}

		selMass[[i]] <- sm
		topMass[[i]] <- unique(c(sm,tm))
	}

	if(is.null(M$Win_1)) w1 <- rep(RI_dev[1], length(M$RI)) else w1 <- M$Win_1
	if(is.null(M$Win_2)) w2 <- rep(RI_dev[2], length(M$RI)) else w2 <- M$Win_2
	if(is.null(M$Win_3)) w3 <- rep(RI_dev[3], length(M$RI)) else w3 <- M$Win_3
	qM <- if(is.null(M$QUANT_MASS)) numeric(nrow(M)) else as.numeric(M$QUANT_MASS)

	options(op)

	new("tsLib", Name = M$Name, RI = M$RI, medRI = M$RI, RIdev = cbind(w1,w2,w3),
		selMass = selMass, topMass = topMass, quantMass=qM, spectra = spectra, libData = M)
	
}

Top.Masses <- function(sp, TopMasses, ExcludeMasses = NULL) {
	m <- sp[order(sp[,2], decreasing = T),1]
	if(!missing(ExcludeMasses) & !is.null(ExcludeMasses)) {
			if(!is.numeric(ExcludeMasses))
				stop("'ExcludeMasses' must be a numeric vector")
			m <- m[!m %in% ExcludeMasses]
	}

	m[1:min(length(m), TopMasses)]
}

Spectra <- function(x) {
	# check spectra format
	if(any(is.na(x))) stop("Error: Invalid spectrum format.")
	if(length(grep("[^0-9 :]", x))) stop("Error: Invalid spectrum format.")

	# remove extra spaces
	x <- gsub(" +", " ", x)
	x <- sub("^ ","", x)
	x <- sub(" $","", x)

	y <- strsplit(as.character(x), " +")
	z <- sapply(y, function(x) unlist(strsplit(x, ":")), simplify=FALSE)
	z <- sapply(z, function(x) matrix(as.numeric(x), ncol=2, byrow=TRUE), simplify=FALSE)
	return(z)
}

# function to read a msp file. Returns a list where every component represent
# a metabolite/analyte entry in the MSP. Every entry is a list with
# following components:
#   Name: The metabolite/analite name:
#   formula: the chemical formula (NA if no "Formula:" field is found)
#   mw: Molecular weigth (NA if no "MW:" field is found)
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
        RIdev = matrix(rep(RI_dev, length(msp)), ncol = 3, byrow = T),
        selMass = selMass, topMass = topMass,
        spectra = sapply(msp, function(x) x$spectrum), libData = libData)
}

# function to remove unexpected quotation marks (usually generated by
# external softwares such as excel. Also converts all factors to character

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
