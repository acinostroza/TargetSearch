require(TargetSearchData)

smp <- new('tsSample', RIfiles=tsd_rifiles(), CDFfiles=tsd_cdffiles())
lib <- ImportLibrary(file.path(tsd_data_path(), 'library.txt'))

# check for invalid usage
expect_error(ret <- FindAllPeaks(smp))
expect_error(ret <- FindAllPeaks(smp, lib))

# check for leucine (both are equivalent)
m <- c(102, 158, 232, 260)
x <- FindAllPeaks(smp, dev=1500, RI=306800, mz=m)
y <- FindAllPeaks(smp, lib, 'GC.5')
expect_equal(x, y)

# check for valine
m <- c(100, 144, 156, 218, 246)
x <- FindAllPeaks(smp, dev=2000, RI=271500, mz=m)
y <- FindAllPeaks(smp, lib, 'GC.3')
expect_equal(x, y)

# check searching by RT
source('mock_functions.R')
peak_search <- function(data, dev, RT, mz)
{
    search <- function(x, fid) {
        k <- which(abs(x$Time - RT) < dev)
        y <- x$Peaks[k,  mz - x$massRange[1] + 1]
        j <- arrayInd(which(y > 0), dim(y))
        f <- if(length(j) > 0) fid else integer(0)
        cbind(Int=y[y > 0], RI=x$Index[k][j[,1]], RT=x$Time[k][j[,1]], mz=mz[j[,2]], fid=f)
    }
    res <- mapply(search, data, seq(data), USE.NAMES=FALSE, SIMPLIFY=FALSE)
    if(nrow(res <- do.call('rbind', res)) == 0)
        return(NULL)
    res
}

tmp <- tempdir()
ncfiles <- file.path(tmp, sprintf("file%d.nc4", 1:3))
smp <- new('tsSample', CDFfiles=ncfiles)

# repeat 10 times with random values
for(i in 1:20) {
    # write files and recover data
    data <- lapply(RIfiles(smp), mock_rifile, mz_range=c(100,200), time_range=c(200, 300), n_peaks=300, n_scans=50)

    RT <- runif(1, 200, 300)
    mz <- sample(100:200, 5)
    x <- FindAllPeaks(smp, dev=10, RT=RT, mz=mz)
    y <- peak_search(data, dev=10, RT=RT, mz=mz)
    expect_identical(x, y)
}
