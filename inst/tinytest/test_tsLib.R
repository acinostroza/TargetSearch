# tests for class tsSample

mkspec <- function(n) {
    f <- function(k)
        cbind(sort(sample(85:500, k)), round(runif(k, 1, 999)))
    k <- round(runif(n, 3, 20))
    lapply(k, f)
}

# make an library object
names <- sprintf('Metab_%d', seq(4))
ri   <- seq(4) * 100
sel  <- lapply(ri, function(x) sample(85:500, 4))
dev  <- matrix(rep(c(10,5,2), length(names)), ncol = 3, byrow = TRUE)
sp   <- mkspec(4)
lib  <- new("tsLib", Name = names, RI = ri, RIdev = dev, selMass = sel, spectra = sp)
expect_true(validObject(lib))
expect_equal(length(lib), 4L)

# split object
a <- lib[1:2]
b <- lib[3:4]
expect_equal(length(a), 2L)
expect_equal(length(b), 2L)

# test combination
x <- c(a,b)

expect_true(validObject(x))
expect_equal(length(x), 4L)
expect_equal(libName(x), libName(lib))
expect_equal(libRI(x), libRI(lib))
expect_equal(medRI(x), medRI(lib))
expect_equal(RIdev(x), RIdev(lib))
expect_equal(spectra(x), spectra(lib))
expect_equal(quantMass(x), quantMass(lib))
