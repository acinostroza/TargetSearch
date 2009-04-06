# Normalisation methods

`dayNorm` <-
function(samples, resInt) {
  run_days <- sampleDays(samples)
  Days <- unique(run_days)
	res  <- matrix(ncol = ncol(resInt), nrow = nrow(resInt))
	for(i in 1:length(Days)) {
    tmp <- resInt[,run_days == Days[i]]
    day.med <- apply(tmp,1, median, na.rm=T)
    tmp <- sweep(tmp,1,day.med, FUN="/")
    res[,run_days == Days[i]] <- tmp 
  }
	return(res)
}

`medianNorm` <-
function(my.files, resInt) {
	res.med <- apply(resInt,1, median, na.rm=T)
	return(sweep(resInt,1,res.med, FUN="/"))
}
