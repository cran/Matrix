## These are tests related to the centralization (since r~3454) of various
## methods for symmetrizing the (possibly asymmetric) 'Dimnames' of symmetric
## matrices.

library(Matrix)

if (interactive()) {
    options(Matrix.verbose = TRUE, warn = 1, error = recover)
} else {
    options(Matrix.verbose = TRUE, warn = 1)
}

## For getting and setting '[dD]imnames' on '[mM]atrix'
DN <- function(x) {
    if (is(x, "Matrix")) {
        x@Dimnames
    } else {
        dimnames(x)
    }
}
`DN<-` <- function(x, value) {
    if (is(x, "Matrix")) {
        x@Dimnames <- value
    } else {
        dimnames(x) <- value
    }
    x
}

## SDN1(dn) is documented to behave as SDN2(dn, NULL)
SDN1 <- Matrix:::symmDN
SDN2 <- function(dn, uplo = NULL) {
    J <-
        if (is.null(uplo)) {
            if (!is.null(dn[[1L]]) && is.null(dn[[2L]])) 1L else 2L
        } else {
            if (uplo == "U") 2L else 1L
        }
    rep(dn[J], 2L)
}

## isSDN1(dn) is documented to behave as isSDN2(dn)
isSDN1 <- Matrix:::isSymmetricDN
isSDN2 <- function(dn) {
    (is.null(ndn <- names(dn)) || !all(nzchar(ndn)) || ndn[1L] == ndn[2L]) &&
        (is.null(rn <- dn[[1L]]) || is.null(cn <- dn[[2L]]) ||
         isTRUE(all(rn == cn | (is.na(rn) & is.na(cn)))))
}

## Various possible (a)symmetries of 'Dimnames'
n <- 4L
rn <- letters[seq_len(n)]
cn <- LETTERS[seq_len(n)]
ldn <- list(list(rn, rn),
            list(rn, cn),
            list(rn, NULL),
            list(NULL, cn),
            list(NULL, NULL),
            list(x = rn, rn),
            list(x = rn, cn),
            list(x = rn, NULL),
            list(x = NULL, cn),
            list(x = NULL, NULL),
            list(rn, y = rn),
            list(rn, y = cn),
            list(rn, y = NULL),
            list(NULL, y = cn),
            list(NULL, y = NULL),
            list(x = rn, y = rn),
            list(x = rn, y = cn),
            list(x = rn, y = NULL),
            list(x = NULL, y = cn),
            list(x = NULL, y = NULL))

## 'matrix' and _most_ 'd..Matrix' ...
## zero matrices are fine for the purpose of testing handling of 'Dimnames'
lM <- c(list(matrix(0, n, n),
             new("ddiMatrix", x = double(n), Dim = c(n, n)),
             new("dgeMatrix", x = double(n * n), Dim = c(n, n))),
        .mapply(new,
                expand.grid(Class = c("dsyMatrix", "dtrMatrix"),
                            uplo = c("U", "L"),
                            stringsAsFactors = FALSE),
                list(x = double(n * n), Dim = c(n, n))),
        .mapply(new,
                expand.grid(Class = c("dspMatrix", "dtpMatrix"),
                            uplo = c("U", "L"),
                            stringsAsFactors = FALSE),
                list(x = double((n * (n + 1L)) %/% 2L), Dim = c(n, n))),
        list(new("dgCMatrix", x = double(0L), Dim = c(n, n),
                 i = integer(0L), p = rep.int(0L, n + 1L))),
        .mapply(new,
                expand.grid(Class = c("dsCMatrix", "dtCMatrix"),
                            uplo = c("U", "L"),
                            stringsAsFactors = FALSE),
                list(x = double(0L), Dim = c(n, n),
                     i = integer(0L), p = rep.int(0L, n + 1L))))

## A few dense symmetric matrices, which are _not_ symmetricMatrix
## and whose symmetry (in the sense of 'isSymmetric') should depend
## only on their 'Dimnames' slot
.d <- diag(n)
.lM <- list(new("dgeMatrix",
                x = as.vector(.d), Dim = c(n, n)),
            new("ltrMatrix",
                x = as.vector(.d != 0), Dim = c(n, n), uplo = "U"),
            new("ntpMatrix",
                x = .d[upper.tri(.d, TRUE)] != 0, Dim = c(n, n), uplo = "U"))
.iS <- function(M, dn) {
    M@Dimnames <- dn
    isSymmetric(M, tol = 0, checkDN = TRUE)
}

for (dn in ldn) {
    stopifnot(identical(sdn <- SDN1(dn), SDN2(dn)),
              (isdn <- isSDN1(dn)) == isSDN2(dn),
              vapply(.lM, .iS, NA, dn = dn) == isdn)

    for (M in lM) {
        DN(M) <- dn
        if (is.s <- is(M, "symmetricMatrix")) {
            ## 'dimnames' should symmetrize
            stopifnot(identical(dimnames(M), sdn))
        }

        if (is.s && !identical(dn[1L], dn[2L])) {
            ## Methods for 'symmetricMatrix' assume symmetric 'Dimnames'
            ## for efficiency ... should they?
            next
        }
        stopifnot(identical(DN(forceSymmetric(M)), sdn),
                  identical(DN(symmpart(M)),       sdn),
                  identical(DN(skewpart(M)),       sdn))
        ## others?
    }
}

## r3459: allowing initialization with typeof(Dimnames[[i]]) != "character"
## ... nothing to do with symmetry, but here for now ...
stopifnot(identical(new("dgeMatrix", x = as.double(1:4), Dim = c(2L, 2L),
                        Dimnames = list(1:2, as.factor(3:4))),
                    new("dgeMatrix", x = as.double(1:4), Dim = c(2L, 2L),
                        Dimnames = list(c("1", "2"), c("3", "4")))))

stopifnot(vapply(ldn, isSDN1, NA) == vapply(ldn, isSDN2, NA))

cat("Time elapsed:", proc.time(), "\n") # "stats"
