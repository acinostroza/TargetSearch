# test mock data
source('mock_functions.R')

# functions that re-implement ri_plot_peak and ri_data_extract
format_data <- function(x) {
    p <- which(x$Peaks > 0)
    p <- arrayInd(p, dim(x$Peaks))
    m <- x$massRange[1]:x$massRange[2]
    cbind(RI=x$Index[p[,1]], RT=x$Time[p[,1]], Intensity=x$Peaks[p], mz=m[p[,2]])
}

search_data <- function(data, mz_values, time_range) {
    z <- lapply(data, function(x)  x[ x[, 'RT'] > time_range[1] & x[, 'RT'] < time_range[2] & x[, 'mz'] %in% mz_values,,drop=FALSE])
    n <- rep(seq(data), sapply(z, nrow))
    cbind(do.call(rbind, z), sample=n)
}

repeats <- 10
nfiles <- 10
mz_range <- c(100, 110)
time_range <- c(300, 400)

data <- lapply(seq(nfiles), function(i) mock_peak_data(mz_range, time_range, n_peaks=200, n_scans=1001, ri_factor=1004))
files <- sapply(data, mock_write_file)
data <- lapply(data, format_data)

for(i in seq(repeats)) {
    mz <- sort( sample(seq(mz_range[1], mz_range[2]), 3) )
    tm <- runif(1,time_range[1], time_range[2]) + c(-5, 5)

    # run ri - plot function
    fake <- search_data(data, mz, tm)
    if(nrow(fake) == 0) # do not check if we do not get any peak
        next
    expect_silent(real <- ri_plot_peak(files, massValues=mz, timeRange=tm))
    expect_identical(real, fake)
}

# shoud raise a warning and return null
expect_warning(ret <- ri_plot_peak(files, massValues=c(300, 400), timeRange=tm, plot=FALSE))
expect_null(ret)
unlink(files)
