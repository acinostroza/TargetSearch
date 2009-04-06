`findRetentionTime` <-
function(rTime, Intensity, rimLimits) {
  fameTimes <- rep(0, nrow(rimLimits))
  for(i in 1:nrow(rimLimits)) {
		window <- which(rTime > rimLimits[i,1] & rTime < rimLimits[i,2])
    if(length(window) == 0) {
      fameTimes[i] <- NA
      next
    }
    fameTimes[i] <- rTime[window[which.max(Intensity[window])]]
  }
  return(fameTimes)
}

