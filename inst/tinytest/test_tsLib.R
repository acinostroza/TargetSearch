# tests for class tsSample

mkspec <- function(n) {
    f <- function(k)
        cbind(sort(sample(85:500, k)), round(runif(k, 1, 999)))
    k <- round(runif(n, 3, 20))
    lapply(k, f)
}

# make an library object
mk_lib <- function(prefix, size, data=NULL)
{
    names <- sprintf('%s_%d', prefix, seq(size))
    ri   <- seq(size) * 100
    sel  <- lapply(ri, function(x) sample(85:500, size))
    dev  <- matrix(rep(c(10,5,2), length(names)), ncol = 3, byrow = TRUE)
    sp   <- mkspec(size)
    lib  <- new("tsLib", Name = names, RI = ri, RIdev = dev, selMass = sel, spectra = sp, libData = data)
}

lib <- mk_lib('Metab', 4)
expect_true(validObject(lib))
expect_equal(length(lib), 4L)

# split object
a <- lib[1:2]
b <- lib[3:4]
expect_equal(length(a), 2L)
expect_equal(length(b), 2L)

# test combination
expect_silent(x <- c(a,b))

expect_true(validObject(x))
expect_equal(length(x), 4L)
expect_equal(libName(x), libName(lib))
expect_equal(libRI(x), libRI(lib))
expect_equal(medRI(x), medRI(lib))
expect_equal(RIdev(x), RIdev(lib))
expect_equal(spectra(x), spectra(lib))
expect_equal(quantMass(x), quantMass(lib))

# test spectra assignment methods
expect_error(spectra(lib) <- list(1,2, 3, 4))
expect_silent(spectra(lib) <- list())
expect_silent(spectra(lib) <- NULL)
expect_silent(spectra(lib) <- mkspec(4))
expect_silent(spectra(lib)[[2]] <- numeric(0))

# test combination with unequal data slot
da <- data.frame(libID=paste0('GC', 101:103),
                 Analyte=c('A', 'B', 'C'),
                 Formula=c('C1', 'C2', 'C3'))

db <- data.frame(libID=paste0('GC', 201:207),
                 Metabolite=paste0('MM_', 1:7))

expect_silent(la <- mk_lib('metab_A', 3, da))
expect_silent(lb <- mk_lib('metab_B', 7, db))
expect_silent(lib <- c(la, lb))
expect_equal(nrow(da), length(la))
expect_equal(nrow(db), length(lb))
expect_equal(length(lib), length(la) + length(lb))
