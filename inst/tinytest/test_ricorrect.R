## load common functions
source('mock_functions.R')

#' simulate the profile a RI marker
#' @param t the time vector
#' @param x the time of the markers
#' @param e error term to introduce random variation
sim_rim <- function(t, x, e=2) {
    n <- length(x)
    y <- x + runif(n, -e, e)
    A <- 10^runif(n, 4, 6)
    s <- runif(n, 0.5, 1.5)
    z <- mapply(function(xx, aa, ss) gauss_peak(t, xx, aa, ss), y, A, s)
    z <- round(rowSums(z))
    i <- sapply(x, function(u) {
                    j <- which(abs(t - u) < 3*e)
                    j[ which.max(z[j]) ]
         })
    list(z=z, rim=t[i], rim_int=z[i])
}

## generate random CDF files with expected RI markers

nfiles <- 7
nfames <- 5
mz_range <- c(85, 90)
time_range <- c(200, 300)
rim_time <- seq(210, 290, length=nfames)
RImat <- RIint <- matrix(nrow=nfames, ncol=nfiles)
cdf <- character(nfiles)

for(i in seq(nfiles)) {
    z <- mock_ncdf_data(mz_range, time_range, n_peaks=20, n_scans=501, ri_factor=1004)
    w <- sim_rim(z$Time, rim_time)
    z$Peaks[, 3] <- as.integer(w$z) # m/z 87
    z$Index <- NULL ; z$baselineCorrected <- FALSE
    cdf[i] <- mock_write_cdf(z, fileext='.nc4')
    RImat[, i] <- w$rim
    RIint[, i] <- w$rim_int
}

## create random sample and rim objects

expect_silent(smp <- new('tsSample', CDFfiles=cdf, days=1) )

limit <- cbind(rim_time - 5, rim_time + 5)
std   <- 100 * seq(nfames)
expect_silent(rim <- new('tsRim', limits=limit, standard=std, mass=87))

## performs RT correction and check we have the expected values

RI <- RIcorrect(smp, rim, Window=5, IntThreshold=100)
expect_equivalent(RImat, RI)
expect_equivalent(RIint, attr(RI, 'intensity'))

## check riMatrix function

RI2 <- riMatrix(smp, rim)
expect_equal(RI2, RI)

unlink(CDFfiles(smp))
unlink(RIfiles(smp))
