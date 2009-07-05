# A function to manually correct the RIs 
fixRIcorrection <- function(samples, rimLimits, RImatrix, sampleNames) {

   if( all( colnames(RImatrix) == sampleNames(samples) ) == FALSE) {
        stop("Sample names and columns of the RI matrix do not match")
   }

   if(is.character(sampleNames)) {
        idx <- which(sampleNames(samples) %in% sampleNames)
		if(length(idx) <= 0) {
		    stop("sampleNames not found")
		}
   } else if(is.numeric(sampleNames)) {
        idx <- sampleNames
   } else {
        stop("Invalid sampleNames argument")
   }
   
   ri.files <- RIfiles(samples)
   standard  <- rimStandard(rimLimits)
   
   for(i in idx) {
       fameTimes <- RImatrix[,i]
       message(sprintf("Correcting File %s", ri.files[i]))
       tmp    <- read.delim(ri.files[i], as.is = TRUE)
  	   tmp$RETENTION_TIME_INDEX <- rt2ri(tmp$RETENTION_TIME, fameTimes, standard)
       write.table(tmp, file = ri.files[i], row.names = FALSE, sep="\t", quote=FALSE)
   }
}
