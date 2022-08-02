# Normalisation methods

`dayNorm` <-
function(samples, resInt) {
    run_days <- sampleDays(samples)
    lapply(resInt, function(z) {
        Days <- unique(run_days)
    	res  <- z
        for(i in 1:length(Days)) {
            tmp <- z[,run_days == Days[i],drop=FALSE]
            day.med <- apply(tmp,1, median, na.rm=TRUE)
            tmp <- sweep(tmp,1,day.med, FUN="/")
            res[,run_days == Days[i]] <- tmp
        }
	   return(res)
    })
}

`medianNorm` <-
function(my.files, resInt) {
    lapply(resInt, function(z) {
        res.med <- apply(z,1, median, na.rm=TRUE)
    	sweep(z,1,res.med, FUN="/")
    })
}
