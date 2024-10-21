`FAMEoutliers` <-
function(samples, RImatrix, pdffile=NA, startDay=NULL, endDay=NULL,
         threshold=3, group.threshold=0.05)
{
  days <- sampleDays(samples)
  days.unique <- unique(days)

  manyFiles <- CDFfiles(samples)

    # Algorithm to Identify DayGroups
    if(is.null_or_na(startDay) && is.null_or_na(endDay)) {
        # get the daily mean RI per R.T. Standard
        m <- vapply(days.unique, function(d) {
                        rowMeans(RImatrix[, d == days, drop=FALSE], na.rm=TRUE)
             }, numeric(nrow(RImatrix)))

		# Identify the day groups
		day.groups <- 1
		if(length(days.unique) >= 2) {
			m.dist <- dist(t(m))
			hc <- hclust(m.dist, "ave")
			day.groups <- cutree(hc, h = group.threshold * max(m.dist))
		}
    }
    else {
        day.groups <- numeric(length(days.unique))
        day.index  <- seq(length(days.unique))
        assert_that(!has_na(startDay), msg="`startDay` cannot contain NAs")
        assert_that(!has_na(endDay), msg="`endDay` cannot contain NAs")
        assert_that(length(startDay) == length(endDay), length(startDay) > 0)
        st <- match(startDay, days.unique)
        en <- match(endDay, days.unique)
        assert_that(!has_na(st), !has_na(en),
                  msg="Unmatched days found in `startDay` and `endDay`")
        for(k in seq(length(startDay)))
            day.groups[ day.index >= st[k] & day.index <= en[k] ] <- k

		if(!all(day.groups > 0))
			stop("Check arguments 'startDay' & 'endDay', at least one of the measurement days is out of range.")
    }

	# find the outliers.
	RI_out <- matrix(FALSE, ncol = ncol(RImatrix), nrow = nrow(RImatrix))
	colnames(RI_out) <- colnames(RImatrix)
	rownames(RI_out) <- rownames(RImatrix)
	for(i in 1:max(day.groups)) {
		RI_tmp <- RImatrix[, days %in% days.unique[day.groups == i], drop=FALSE]
		RI_tmp.mean <- rowMeans(RI_tmp, na.rm=TRUE)
		RI_tmp.sd <- apply(RI_tmp, 1, sd, na.rm=TRUE)
		RI_out[, days %in% days.unique[day.groups == i]] <-
		RI_tmp < RI_tmp.mean - threshold*RI_tmp.sd | RI_tmp > RI_tmp.mean + threshold*RI_tmp.sd
	}
    RI_out[is.na(RI_out)] <- FALSE

	# find missing markers
	missingMarkers <- missingIndex <- NULL
	if(any(is.na(RI_out))) {
		message("\n Missing Markers: \n ===============\n\n")
		missingIndex   <- which(is.na(RI_out))
		for(j in 1:nrow(RI_out)) {
			missingMarkers <- append(missingMarkers, sprintf(" RT standard: %2d | Sample: %s", j, manyFiles[is.na(RI_out[j,])]))
		}
		cat(missingMarkers, sep="\n")
		# set missing markers temporarily to FALSE
		RI_out[missingIndex] <- FALSE
	}

	# print the outliers:
	cat("\n Outliers Report: \n ===============\n\n")
	outliers <- NULL
	for(j in 1:nrow(RI_out)) {
		if(any(RI_out[j,])) {
			outliers <- append(outliers, sprintf(" RT standard: %2d | Sample: %s", j, manyFiles[RI_out[j,]]))
		}
	}
	cat(outliers, sep = "\n")
	if(!any(RI_out)) cat("No outliers were found.\n")

 	if(is.na(pdffile))
		pdffile <- paste("TargetSearch-", Sys.Date(), ".FAME-report.pdf", sep = "")

	new.page <- function(x="Outliers Report") {
	  plot.new()
	  title(x)
	}
	# print outliers in FAME report file
	pdf(pdffile, width = 8, height = 8)
	on.exit(dev.off())
	if(length(outliers) > 0) {
		npages <- ceiling(length(outliers)/45)
		op <- par(mar = c(3,3,3,2))
			for(j in 1:length(outliers)) {
				if(j %% 45 == 1) new.page()
				text(-0.03, 1.05-((j-1) %% 45 + 1)/45 , outliers[j], adj = 0, cex = 0.8)
			}
		par(op)
	}

	# print missing RI markers, if any.
	if(length(missingMarkers) > 0) {
		npages <- ceiling(length(missingMarkers)/45)
		op <- par(mar = c(3,3,3,2))
			for(j in 1:length(missingMarkers)) {
				if(j %% 45 == 1) new.page('Missing RI Markers')
				text(-0.03, 1.05-((j-1) %% 45 + 1)/45 , missingMarkers[j], adj = 0, cex = 0.8)
			}
		par(op)
	}

	for(i in 1:nrow(RImatrix)) {
		plotFAME(samples, RImatrix, i)
		points(which(RI_out[i,]), RImatrix[i,RI_out[i,]], cex = 2, col = "blue")
	}
	message("FAMEs were saved in ", pdffile)

	# set missing markers to NA
	if(length(missingIndex) > 0)
		RI_out[missingIndex] <- NA

	return(RI_out)
}

