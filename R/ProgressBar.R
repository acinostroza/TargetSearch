# functions to display a progress bar.
# ChangeLog 14.11.2011
#  - rewritten
# ChangeLog 13.11.2011
#  - fixed progress bar bug
# ChangeLog 16.01.2009.
#  - winProgressBar option removed. Reason: R CMD check doesn't work in linux.

ProgressBar <- function(title, label) {
	message("  ", title)
	pb <- txtProgressBar(style = 3, width = 72)
	return(pb)
}

setProgressBar <-
function(pb, value, title, label) {
	setTxtProgressBar(pb, value)
}
