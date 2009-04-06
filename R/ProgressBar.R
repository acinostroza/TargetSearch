# functions to display a progress bar.

# ChangeLog 16.01.2009.
#  - winProgressBar option removed. Reason: R CMD check doesn't work in linux.

ProgressBar <- function(title, label) {
#	pb <- NULL
#	if(exists("winProgressBar")) {
#	  pb <- winProgressBar(title, label)
#	} else if (exists("txtProgressBar")) {
	  pb <- txtProgressBar(style = 3, width = 72)
	  message("  ", title)
#	}
	return(pb)
}

setProgressBar <-
function(pb, value, title, label) {
#	if(exists("setWinProgressBar")) {
#	  setWinProgressBar(pb, value, title, label)
#	} else if (exists("setTxtProgressBar")) {
	  setTxtProgressBar(pb, value)
#	}
}
