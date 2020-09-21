
mock_peak_data <- function(mz_range, time_range, n_peaks, n_scans, ri_factor=1004)
{
    stopifnot(length(mz_range) == 2)
    stopifnot(length(time_range) == 2)
    P <- matrix(0, nrow=n_scans, ncol=diff(mz_range) + 1)
    j <- sample(length(P), n_peaks)
    P[j] <- round( 10^runif(n_peaks, 2, 6))
    RT <- seq(time_range[1], time_range[2], length=n_scans)
    RI <- RT * ri_factor
    list(Time=RT, Index=RI, massRange=mz_range, Peaks=P)
}

mock_write_file <- function(data)
{
    outfile <- tempfile()
    TargetSearch:::writeRIFile(outfile, data, data$Index, data$massRange, 'binary')
    outfile
}

gauss_peak <- function(t, t0, A, s) A * exp(-(t-t0)^2/2/s^2)

gen_rand_trace <- function(t, n) {
    s <- runif(n, 0.5, 1.5)
    A <- 10^runif(n, 2, 6)
    t0 <- runif(n, min(t), max(t))
    as.integer(round(sapply(t, function(t) sum(gauss_peak(t, t0, A, s)))))
}

mock_ncdf_data <- function(mz_range, time_range, n_peaks, n_scans, ri_factor=1004)
{
    stopifnot(length(mz_range) == 2)
    stopifnot(length(time_range) == 2)
    t <- seq(time_range[1], time_range[2], length=n_scans)
    n <- diff(mz_range) + 1
    p <- sapply(sample(n_peaks, n, replace=T), function(k) gen_rand_trace(t, k))
    list(Time=t, Peaks=p, massRange=mz_range, Index=t * ri_factor, baselineCorrected=TRUE)
}

mock_write_cdf <- function(peaks)
{
    outfile <- tempfile()
    TargetSearch:::ncdf4_write(outfile, peaks)
    outfile
}
