`ImportFameSettings` <-
function(file, mass=NULL, standard=NULL, colnames=FALSE, ...)
{
    assert_that(is.flag(colnames))

  ## next line specifies the correct RI-value's of your marker metabolites
  ## the values are for a LECO machine with fatty acid methyl ester (FAME)
  ## standard
  rim.perfect <- c(262320, 323120, 381020, 487220, 582620, 668720, 747420,
                   819620, 886620, 948820, 1006900, 1061700, 1113100)

    ## default m/z value for FAME
    rim.mass <- 87

    ## default retention time limits
    rim.limits <- rbind(c( 230,  280),   # lower/upper/RI_correct_value limit of first Marker
                        c( 290,  340),   # lower/upper limit of second Marker
                        c( 350,  400),   # ...
                        c( 450,  500),
                        c( 540,  590),
                        c( 630,  665),
                        c( 705,  745),
                        c( 775,  820),
                        c( 845,  885),
                        c( 905,  940),
                        c( 965, 1015),
                        c(1020, 1060),
                        c(1070, 1110))

    ## expected column names
    rim.cols <- c('LowerLimit', 'UpperLimit', 'RIstandard', 'mass')

    ## here you set up the search frames for these markers (in seconds)
    ## if you specify a path for 'file' you can load settings from a tab-delimited file
    if(!missing(file)) {
        if(is.string(file))
            dat <- read.delim(file, ...)
        else if(is.data.frame(file) || is.matrix(file))
            dat <- file
        else
            stop('Invalid object `file`: it must be either a file path or a data.frame',
                 ' or a matrix.')

        ## match column names of fame file
        if(colnames) {
            cols <- match(tolower(rim.cols), tolower(colnames(dat)))

            if(any(is.na(cols[1:2])))
                stop('Could not find columns `', rim.cols[1], '` and `', rim.cols[2], '`.')

            rim.limits <- as.matrix(dat[, cols[1:2]])
            if(!is.na(cols[3]))
                rim.perfect <- dat[, cols[3]]

            if(!is.na(cols[4]))
                rim.mass <- dat[, cols[4]]
        }

        ## use column positions
        else if(ncol(dat) >= 2) {
            rim.limits <- as.matrix(dat[, 1:2])

            if(ncol(dat) >= 3)
                rim.perfect <- dat[,3]
            if(ncol(dat) >= 4)
                rim.mass <- dat[, 4]
        }
        else {
            stop("Error reading FAME file. The file doesn't have 2, 3 or 4 columns",
                 " (LowerLimits, UpperLimit, RIstandard, mass)")
        }
    }

    if(any(rim.limits[,1] > rim.limits[,2]))
        stop("Error: LowerLimits are greater than UpperLimits. Please check your file (rows ",
              paste(which(rim.limits[,1] > rim.limits[,2]), collapse = " "), ")")

    colnames(rim.limits) <- rim.cols[1:2]

    ## override mass and standard values
    if(is.null(mass)) mass <- rim.mass
    if(is.null(standard)) standard <- rim.perfect

    new("tsRim", limits = rim.limits, standard = standard, mass = mass)
}
