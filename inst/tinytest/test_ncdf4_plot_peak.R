# test mock data
source('mock_functions.R')

repeats <- 5
nfiles <- 2
mz_range <- c(100, 110)
time_range <- c(300, 400)

search_data <- function(data, mz_values, time_range) {
    sapply(data, function(x) {
                    i <- which( x$Time > time_range[1] & x$Time < time_range[2] )
                    j <- match(mz_values, x$massRange[1]:x$massRange[2])
                    p <- x$Peaks[i, j]
                    colnames(p) <- (x$massRange[1]:x$massRange[2])[j]
                    list(Time=as.array(x$Time[i]), Index=as.array(x$Index[i]), Intensity=p,
                         massRange=as.array(as.integer(x$massRange)))
                }, simplify=FALSE)
}


data <- lapply(seq(nfiles), function(i) mock_ncdf_data(mz_range, time_range, n_peaks=20, n_scans=1001, ri_factor=1004))
files <- sapply(data, mock_write_cdf)
names(data) <- basename(files)

for(i in seq(repeats)) {
    mz <- sort( sample(seq(mz_range[1], mz_range[2]), 2) )
    tm <- runif(1,time_range[1], time_range[2]) + c(-5, 5)
    fake <- search_data(data, mz, tm)
    real <- ncdf4_plot_peak(files, massValues=mz, timeRange=tm, plot=FALSE)
    expect_identical(fake, real)
}

expect_warning(real <- ncdf4_plot_peak(files, massValues=mz, timeRange=tm + 200, plot=FALSE))
expect_null(real)
unlink(files)
