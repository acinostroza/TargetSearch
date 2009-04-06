#

FindPeaks <-
function(my.files, refLib, columns = c("SPECTRUM", "RETENTION_TIME_INDEX"), showProgressBar = FALSE) {

	get.columns <- function(my.file, columns) {
		if(is.character(columns)) {
			header <- scan(my.file, what = "character", nlines = 1, quiet = TRUE)
			tmp <- sapply(columns, function(x) which( header == x ))
			if(length(unlist(tmp)) != length(columns))
				stop("Column name not found. Check your RI file.")
			columns <- unlist(tmp) - 1
		}
		return(columns)
	}
	
	my.names <- basename(my.files)
  resInt   <- matrix(nrow = nrow(refLib), ncol = length(my.files))
  colnames(resInt) <- my.names
  rownames(resInt) <- rownames(refLib)
  resRI <- resInt

	if(showProgressBar)  
	  pb <- ProgressBar(title="Finding Peaks...", label="File in processing...")
  
  for (i in 1:length(my.files)) {
  	if(showProgressBar)
    setProgressBar(pb, value=i/length(my.files),
			title=paste("Findind Peaks (", round(100*i/length(my.files)), "%)"),
			label=sprintf("Reading File %s", basename(my.files[i])))
  
  	cols <- get.columns(my.files[i], columns)
		out <- .Call("FindPeaks", as.character(my.files[i]), as.integer(refLib[,1]),
			as.integer(refLib[,2]),	as.integer(refLib[,3]), as.integer(cols), PACKAGE="TargetSearch")
	  resInt[, i] <- out[[1]]
    resRI[, i]  <- out[[2]]
  }
  if(showProgressBar)
		close(pb)
  return(new("tsMSdata", RI = resRI, Intensity = resInt))
}
