
# functions to create dummy matrices and profile info
make_matrix <- function(a, b) {
    z <- matrix(round(runif(a * b)*1000), a, b)
    rownames(z) <- sample(85:500, a)
    colnames(z) <- paste0('X', seq(b))
    z
}

make_info <- function(n) {
    a <- paste("Names ", seq(n))
    b <- round(runif(n)) * 1000
    data.frame(Names=a, Values=b)
}

#
# Test tsMSdata objects
#

N <- 5 # number of metabolites
S <- 6 # number of samples
M <- sample(3:10, N, replace=TRUE)

x <- lapply(M, make_matrix,  S)

expect_silent(dat <- new('tsMSdata', RI=x, RT=x, Intensity=x))
expect_true(validObject(dat))

y <- make_matrix(3, 4)
expect_error(retTime(dat)[[4]] <- y)
dat@RI[[4]] <- y
expect_error(validObject(dat))

expect_silent(dat <- new('tsMSdata', RI=x, RT=x, Intensity=x))
expect_error(retTime(dat)[[4]] <- dat@RT[, 1:3])
dat@RT[[1]] <- dat@RT[[1]][, 1:3]
expect_error(validObject(dat))

#
# Test tsProfile objects
#

expect_silent(dat <- new('tsMSdata', RI=x, RT=x, Intensity=x))
y <- make_matrix(N, S)
nfo <- make_info(N)

expect_silent(dat <- new('tsProfile', dat, info=nfo, profRT=y, profRI=y, profInt=y))
expect_true(validObject(dat))

dat2 <- dat
dat2@RT[[1]] <- dat2@RT[[1]][, 1:3]
expect_error(validObject(dat2))

dat2 <- dat
dat2@profRT <- dat2@profRT[, 1:3]
expect_error(validObject(dat2))
