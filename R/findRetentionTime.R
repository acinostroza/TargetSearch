findRetentionTime <- function(rTime, Intensity, rimLimits) {
        nc <- if(is.null(dim(Intensity))) 1 else dim(Intensity)[2]
        if(nc != 1 & nc != nrow(rimLimits))
                stop("Error in 'findRetentionTime': # of columns of Intensity must be",
                " the same as number of rows of rimLimits")

        fameTimes <- numeric(nrow(rimLimits))
        for(i in 1:nrow(rimLimits)) {
                window <- which(rTime > rimLimits[i,1] & rTime < rimLimits[i,2])
                if(length(window) == 0) {
                        fameTimes[i] <- NA
                        next
                }
                fameTimes[i] <- if(nc == 1) rTime[window[which.max(Intensity[window])]]
                        else rTime[window[which.max(Intensity[window,i])]]
        }
        return(fameTimes)
}

