# tests for quantMatrix function

# function to make a library
mk_lib <- function(nlib) {
    names <- sprintf('Metab_%d', seq(nlib))
    ri   <- seq(nlib) * 100
    sel  <- lapply(ri, function(x) sample(85:500, 4))
    top  <- lapply(ri, function(x) sample(85:500, 4))
    new("tsLib", Name = names, RI = ri, selMass = sel, topMass = top)
}

# function to simulate a matrix
mk_mat <- function(rnames, cnames)
{
    m <- matrix(0, length(rnames), length(cnames),
                dimnames=list(rnames, cnames))
    m[] <- as.integer(runif(length(m), 100, 10000))
    k <- round(length(m) / 3)
    m[ sample(length(m), k) ] <- NA
    m
}

# function to simulate correlation masses
mk_cormass <- function(x) {
    n <- sample(length(x), 1)
    z <- if(n == 1) x[1] else sample(x, n)
    paste(z, collapse=";")
}

# function to make a fake profile object
mk_prof <- function(lib, nsamp=7) {
    sname <- sprintf("samp_%d", seq(nsamp))
    Int <- lapply(topMass(lib), mk_mat, sname)
    pkdata <- new('tsMSdata', Intensity=Int, RI=Int, RT=Int)
    cormass <- vapply(topMass(lib), mk_cormass, "")
    nfo <- cbind(libData(lib), Masses=cormass)
    prof <- new('tsProfile', pkdata, info=nfo)
}

lib <- mk_lib(10)
prof <- mk_prof(lib)
expect_true(validObject(lib))
expect_true(validObject(prof))

# actual testing
Q <- quantMatrix(lib, prof)
E <- t(sapply(Intensity(prof), function(x) x[1,]))
expect_equivalent(Q, E)
expect_true(all(attr(Q, 'isSelMass')))

Q <- quantMatrix(lib, prof, 'maxint')
E <- t(sapply(Intensity(prof),
              function(x) {
                  z <- apply(x, 1, median, na.rm=TRUE)
                  x[ which.max(z), ]
              }))
expect_equivalent(Q, E)

Q <- quantMatrix(lib, prof, 'maxobs')
E <- t(sapply(Intensity(prof),
              function(x) {
                  z <- rowSums(!is.na(x))
                  x[ which.max(z), ]
              }))
expect_equivalent(Q, E)
