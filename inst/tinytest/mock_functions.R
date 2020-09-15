
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
