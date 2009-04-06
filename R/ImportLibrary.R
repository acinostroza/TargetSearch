# function to import a library:
# Arguments
#  - libfile: The library file
#  - RI_dev: RI deviation
#  - SelMasses: How many masses should be used as selective if no SelMasses are found.
#               Masses are taken from the TopMasses
#  - TopMasses: How many masses from the spectra should be uses as top Masses.                                                                                  
#  - ExcludeMasses: Don't take this masses as TopMasses automatically.

ImportLibrary <- function(libfile, RI_dev = c(2000,1000,200),
	SelMasses = 5, TopMasses = 15, ExcludeMasses) {

	if(missing(libfile))
		stop("argument \"libfile\" is missing, with no default")
		
	M <- read.delim(libfile, as.is = T)

  if(is.null(M$Name))     stop("Column 'Name' is missing!!")
	if(is.null(M$RI))       stop("Column 'RI' is missing!!")

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
				tm <- TopMasses(spectra[[i]], TopMasses, ExcludeMasses)
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

	new("tsLib", Name = M$Name, RI = M$RI, medRI = M$RI, RIdev = cbind(w1,w2,w3),
		selMass = selMass, topMass = topMass, spectra = spectra, libData = M)
	
}

TopMasses <- function(sp, TopMasses, ExcludeMasses) {
	m <- sp[order(sp[,2], decreasing = T),1]
	if(!missing(ExcludeMasses)) {
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
	z <- sapply(y, function(x) unlist(strsplit(x, ":")))
	z <- sapply(z, function(x) matrix(as.numeric(x), ncol = 2, byrow = T))
	return(z)
}

