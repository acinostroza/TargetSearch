`FAMEoutliers` <-
function(samples, RImatrix, pdffile = NA, startDay = NA,
	endDay = NA, threshold = 3, group.threshold = 0.05){
  days <- sampleDays(samples)
  days.unique <- unique(days)
  
  manyFiles <- CDFfiles(samples)

	# Algorithm to Identify DayGroups
	if(is.na(startDay) & is.na(endDay)) {
		m <- matrix(nrow = nrow(RImatrix), ncol = length(days.unique))
		# get the daily mean RI per R.T. Standard
		for(d in days.unique) {
			tmp <- RImatrix[, d == days]
			m[, which(days.unique %in% d)] <- apply(tmp, 1, mean)
		}

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
		stopifnot(length(startDay) == length(endDay))
		for(i in 1:length(startDay))
			day.groups[days.unique >= startDay[i] & days.unique <= endDay[i]] <- i
		if(!all(day.groups > 0))
			stop("Check arguments 'startDay' & 'endDay', at least one of the measurement days is out of range.")
	}
	
	# find the outliers.
	RI_out <- matrix(FALSE, ncol = ncol(RImatrix), nrow = nrow(RImatrix))
	colnames(RI_out) <- colnames(RImatrix)
	rownames(RI_out) <- rownames(RImatrix)
	for(i in 1:max(day.groups)) {
		RI_tmp <- RImatrix[, days %in% days.unique[day.groups == i]]
		RI_tmp.mean <- apply(RI_tmp, 1, mean)
	  RI_tmp.sd <- apply(RI_tmp, 1, sd)
		RI_out[, days %in% days.unique[day.groups == i]] <-
		RI_tmp < RI_tmp.mean - threshold*RI_tmp.sd | RI_tmp > RI_tmp.mean + threshold*RI_tmp.sd
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

	new.page <- function() {
	  plot.new()
	  title("Outliers Report")
	}	
	# print outliers in FAME report file
	pdf(pdffile, width = 8, height = 8)
	if(length(outliers) > 0) {	
		npages <- ceiling(length(outliers)/45)
		op <- par(mar = c(3,3,3,2))
			for(j in 1:length(outliers)) {
				if(j %% 45 == 1) new.page()
				text(-0.03, 1.05-((j-1) %% 45 + 1)/45 , outliers[j], adj = 0, cex = 0.8)
			}
		par(op)
	}
  		
	for(i in 1:nrow(RImatrix)) {
		plotFAME(samples, RImatrix, i)
		points(which(RI_out[i,]), RImatrix[i,RI_out[i,]], cex = 2, col = "blue")
	}
	dev.off()
	message("FAMEs were saved in ", pdffile)
	
	return(RI_out)
}

