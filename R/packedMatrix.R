## METHODS FOR CLASS: packedMatrix (virtual)
## dense triangular or symmetric matrices with packed storage
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.pM.subclasses <- names(getClassDef("packedMatrix")@subclasses)

setMethod("unpack", signature(x = "packedMatrix"),
          function(x, ...) .Call(R_dense_as_unpacked, x))

setMethod("pack", signature(x = "packedMatrix"),
          function(x, ...) x)

setMethod("forceSymmetric", signature(x = "packedMatrix", uplo = "missing"),
          function(x, uplo) .Call(packedMatrix_force_symmetric, x, NULL))
setMethod("forceSymmetric", signature(x = "packedMatrix", uplo = "character"),
          function(x, uplo) .Call(packedMatrix_force_symmetric, x, uplo))

## Not all of these .pM.is.* are used, because all packedMatrix inherit
## from symmetricMatrix or triangularMatrix, and those classes have
## their own methods.  They are retained here somewhat for completeness ...

.pM.is.sy <- function(object, checkDN = TRUE, ...) {
    ## backwards compatibility: don't check DN if check.attributes=FALSE
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    ## requiring exact symmetry (fast):
    .Call(packedMatrix_is_symmetric, object, checkDN)
}

.pM.is.sy.dz <- function(object, tol = 100 * .Machine$double.eps,
                         tol1 = 8 * tol, checkDN = TRUE, ...) {
    if (tol <= 0)
        .
    else {
        ## going via all.equal (slow):
        isSymmetric(unpack(object), tol = tol, tol1 = tol1,
                    checkDN = checkDN, ...)
    }
}
body(.pM.is.sy.dz) <-
    do.call(substitute, list(body(.pM.is.sy.dz), list(. = body(.pM.is.sy))))

.pM.is.tr <- function(object, upper = NA, ...)
    .Call(packedMatrix_is_triangular, object, upper)

.pM.is.di <- function(object) .Call(packedMatrix_is_diagonal, object)

## method for     .spMatrix in ./symmetricMatrix.R
## method for [lni]tpMatrix in ./triangularMatrix.R
for (.cl in grep("^[dz]tpMatrix$", .pM.subclasses, value = TRUE))
    setMethod("isSymmetric", signature(object = .cl), .pM.is.sy.dz)

setMethod("isDiagonal", signature(object = "packedMatrix"), .pM.is.di)

rm(.pM.is.sy, .pM.is.sy.dz, .pM.is.tr, .pM.is.di, .cl)

setMethod("t", signature(x = "packedMatrix"),
          function(x)
              .Call(packedMatrix_transpose, x))
setMethod("diag", signature(x = "packedMatrix"),
          function(x, nrow, ncol, names = TRUE)
              .Call(packedMatrix_diag_get, x, names))
setMethod("diag<-", signature(x = "packedMatrix"),
          function(x, value)
              .Call(packedMatrix_diag_set, x, value))

setMethod("symmpart", signature(x = "packedMatrix"),
          function(x) .Call(packedMatrix_symmpart, x))
setMethod("skewpart", signature(x = "packedMatrix"),
          function(x) .Call(packedMatrix_skewpart, x))

setMethod("cov2cor", signature(V = "packedMatrix"),
          function(V) as(forceSymmetric(V), "pcorMatrix"))

rm(.pM.subclasses)
