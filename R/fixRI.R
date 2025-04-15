# A function to manually correct the RIs
# Deprecated! see fixRI
fixRIcorrection <- function(...) {

   .Defunct("fixRI")
}

# function to fix the RI of the files indexed by 'idx'
.fixRIfile <- function(ri.files, RImatrix, standard, idx, quiet)
{
   for(i in idx) {
       fameTimes <- RImatrix[,i]
       if(!quiet)
           message(sprintf("Correcting File %s", ri.files[i]))

       ftype <- file_type(ri.files[i])

       if(ftype == 1) {
           if(is.integer(cols <- get.columns.name()))
               cols <- cols + 1
           tmp  <- read.delim(ri.files[i], as.is = TRUE)
           ri_col <- cols[2] ; rt_col <- cols[3]
           tmp[, ri_col] <- rt2ri(tmp[, rt_col], fameTimes, standard)
           write.table(tmp, file = ri.files[i], row.names = FALSE, sep="\t", quote=FALSE)
       } else if(ftype == 0) {
           z <- readRIBin(ri.files[i])
           z$retIndex <- rt2ri(z$retTime, fameTimes, standard)
           writeRIBin(z, ri.files[i])
       } else {
           stop(sprintf("incorrect file format: %s", ri.files[i]))
       }
   }
}

# validate an RI matrix. returns the same matrix
.validateRImatrix <- function(smp, rim, RImat)
{
	if(length(rimStandard(rim)) != nrow(RImat))
		stop("Invalid RI matrix: number of rows don't match rimLimit object.")
	if(is.null(colnames(RImat)))
		stop("Invalid RI matrix: missing column names")
	if(!all(sampleNames(smp) %in% colnames(RImat)))
		stop("Sample names and columns of the RI matrix do not match")
	RImat[, sampleNames(smp), drop=FALSE]
}

# returns the index of the sampl
.sampNameIndex <- function(smp, snames)
{
	if(is.character(snames))
		which(sampleNames(smp) %in% snames)
	else if(is.numeric(snames))
		snames
	else
		stop("Invalid argument")
}

riMatrix <- function(samples, rim)
{
	ri.files <- RIfiles(samples)
	mass     <- rimMass(rim)
	std      <- rimStandard(rim)
	rLimits  <- rimLimits(rim)
	RImat    <- matrix(nrow=dim(rLimits)[1], ncol=length(ri.files))
	colnames(RImat) <- sampleNames(samples)
	rownames(RImat) <- rownames(rLimits)
	if(length(mass) == 1)
		mass <- rep(mass, dim(rLimits)[1])
	RIint    <- RImat

	for (i in seq_along(ri.files)) {
		out  <- .c_find_peaks(
						as.character(ri.files[i]), # MyFile
						as.integer(mass),          # Mass
						NULL,                      # RI_exp
						as.numeric(rLimits[,1]),   # RI_min
						as.numeric(rLimits[,2]),   # RI_max
						TRUE,                      # useRT
						"maxInt",                  # max intensity
						NULL)
		assert_that(!is.null(out), msg=sprintf("Error processing `%s`", ri.files[i]))
		RImat[out[[4]]+1,i] <- out[[3]]
		RIint[out[[4]]+1,i] <- out[[1]]
	}
	structure(RImat, intensity=RIint, mass=rimMass(rim))
}

fixRI <- function(samples, rimLimits, RImatrix=NULL, sampleNames=NULL, quiet=TRUE)
{
	assert_that(is.flag(quiet))
	if(!is.null(RImatrix))
		RImatrix <- .validateRImatrix(samples, rimLimits, RImatrix)
	else
		RImatrix <- riMatrix(samples, rimLimits)

	idx <- seq(length(samples))
	if(!is.null(sampleNames))
		idx <- .sampNameIndex(samples, sampleNames)

	ri.files <- RIfiles(samples)
	standard  <- rimStandard(rimLimits)
	.fixRIfile(ri.files, RImatrix, standard, idx, quiet)
	invisible()
}

# vim: set ts=4 sw=4:
