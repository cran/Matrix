## METHODS FOR CLASS: unpackedMatrix (virtual) ... and many for base matrices
## dense matrices with unpacked storage
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.upM.subclasses <- names(getClassDef("unpackedMatrix")@subclasses)

.upM.pack <- function(x, ...)
    .Call(R_dense_as_packed, x, NULL, NULL)

.upM.pack.ge <- .m.pack <- function(x, symmetric = NA, upperTri = NA, ...) {
    if(((sna <- is.na(symmetric)) || symmetric) && isSymmetric(x, ...))
        .Call(R_dense_as_packed, x, "U", "")
    else if((sna || !symmetric) &&
            (it <- isTriangular(x, upper = upperTri))) {
        uplo <-
            if(is.na(upperTri))
                attr(it, "kind")
            else if(upperTri)
                "U"
            else "L"
        .Call(R_dense_as_packed, x, uplo, "N")
    } else {
        if(sna)
            stop("matrix is not symmetric or triangular")
        else if(symmetric)
            stop("matrix is not symmetric")
        else stop("matrix is not triangular")
    }
}
body(.m.pack)[[2L]][[3L]] <-
    quote(.Call(R_matrix_as_dense, x, ".sp", "U", NULL))
body(.m.pack)[[2L]][[4L]][[3L]][[3L]] <-
    quote(.Call(R_matrix_as_dense, x, ".tp", uplo, "N"))

setMethod("unpack", signature(x = "unpackedMatrix"),
          function(x, ...) x)
setMethod("unpack", signature(x = "matrix"),
          function(x, ...) .m2dense.checking(x, "."))

setMethod("pack", signature(x = "unpackedMatrix"), .upM.pack)
for(.cl in grep("^.geMatrix$", .upM.subclasses, value = TRUE))
setMethod("pack", signature(x = .cl), .upM.pack.ge)
setMethod("pack", signature(x = "matrix"), .m.pack)

rm(.cl, .upM.pack, .upM.pack.ge, .m.pack)

setMethod("forceSymmetric", signature(x = "unpackedMatrix", uplo = "missing"),
          function(x, uplo) .Call(unpackedMatrix_force_symmetric, x, NULL))
setMethod("forceSymmetric", signature(x = "unpackedMatrix", uplo = "character"),
          function(x, uplo) .Call(unpackedMatrix_force_symmetric, x, uplo))

setMethod("forceSymmetric", signature(x = "matrix", uplo = "missing"),
          function(x, uplo) .Call(R_matrix_as_dense, x, ".sy",  "U", NULL))
setMethod("forceSymmetric", signature(x = "matrix", uplo = "character"),
          function(x, uplo) .Call(R_matrix_as_dense, x, ".sy", uplo, NULL))

.upM.is.sy <- function(object, checkDN = TRUE, ...) {
    ## backwards compatibility: don't check DN if check.attributes=FALSE
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    ## requiring exact symmetry:
    .Call(unpackedMatrix_is_symmetric, object, checkDN)
}

.upM.is.sy.dz <- function(object, tol = 100 * .Machine$double.eps,
                          tol1 = 8 * tol, checkDN = TRUE, ...) {
    ## backwards compatibility: don't check DN if check.attributes=FALSE
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    ## be very fast when requiring exact symmetry
    if(tol <= 0)
        return(.Call(unpackedMatrix_is_symmetric, object, checkDN))
    ## pretest: is it square?
    d <- object@Dim
    if((n <- d[1L]) != d[2L])
        return(FALSE)
    ## pretest: are DN symmetric in the sense of validObject(<symmetricMatrix>)?
    if(checkDN && !isSymmetricDN(object@Dimnames))
        return(FALSE)
    if(n <= 1L)
        return(TRUE)
    object <- .M2gen(object)
    ## now handling n-by-n [dz]geMatrix, n >= 2:

    Cj <- if(is.complex(object@x)) Conj else identity
    ae <- function(check.attributes, ...) {
        ## discarding possible user-supplied check.attributes
        all.equal(..., check.attributes = FALSE)
    }

    ## pretest: outermost rows ~= outermost columns? (fast for large asymmetric)
    ## FIXME: quite inefficient, though, if subsetting must go through "matrix"
    if(length(tol1)) {
        i. <- if (n <= 4L) 1:n else c(1L, 2L, n-1L, n)
        for(i in i.)
            if(!isTRUE(ae(target = object[i, ], current = Cj(object[, i]),
                          tolerance = tol1, ...)))
                return(FALSE)
    }
    ## followed by slower test using 't'
    isTRUE(ae(target = object@x, current = Cj(t(object))@x,
              tolerance = tol, ...))
}

.upM.is.tr <- function(object, upper = NA, ...)
    .Call(unpackedMatrix_is_triangular, object, upper)

.upM.is.di <- function(object)
    .Call(unpackedMatrix_is_diagonal, object)

.m.is.sy <- function(object, tol = 100 * .Machine$double.eps,
                     tol1 = 8 * tol, checkDN = TRUE, ...) {
    ## backwards compatibility: don't check DN if check.attributes=FALSE
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    if(is.logical(object) || is.integer(object) || tol <= 0)
        ## requiring exact symmetry:
        return(.Call(matrix_is_symmetric, object, checkDN))
    if(checkDN && !is.null(dn <- dimnames(object)) && !isSymmetricDN(dn))
        return(FALSE)
    ## discarding possible user-supplied check.attributes:
    iS.m <- function(check.attributes, ...) {
        isSymmetric.matrix(..., check.attributes = FALSE)
    }
    iS.m(object = object, tol = tol, tol1 = tol1, ...)
}

.m.is.tr <- function(object, upper = NA, ...)
    .Call(matrix_is_triangular, object, upper)

.m.is.di <- function(object)
    .Call(matrix_is_diagonal, object)

## method for     .syMatrix in ./symmetricMatrix.R
## method for [lni]trMatrix in ./triangularMatrix.R
for (.cl in grep("^[lni]geMatrix$", .upM.subclasses, value = TRUE))
    setMethod("isSymmetric", signature(object = .cl), .upM.is.sy)
for (.cl in grep("^[dz](ge|tr)Matrix$", .upM.subclasses, value = TRUE))
    setMethod("isSymmetric", signature(object = .cl), .upM.is.sy.dz)

## method for .syMatrix in ./symmetricMatrix.R
## method for .trMatrix in ./triangularMatrix.R
for (.cl in grep("^.geMatrix$", .upM.subclasses, value = TRUE))
    setMethod("isTriangular", signature(object = .cl), .upM.is.tr)

setMethod("isDiagonal", signature(object = "unpackedMatrix"), .upM.is.di)

if(FALSE) {
## Would override isSymmetric.matrix and be faster in the logical and integer
## cases and in the tol<=0 case, but use a looser notion of symmetric 'dimnames'
## and so probably break too much ...
setMethod("isSymmetric", signature(object = "matrix"), .m.is.sy)
}
setMethod("isTriangular", signature(object = "matrix"), .m.is.tr)
setMethod("isDiagonal", signature(object = "matrix"), .m.is.di)

rm(.upM.is.sy, .upM.is.sy.dz, .upM.is.tr, .upM.is.di,
   .m.is.sy, .m.is.tr, .m.is.di, .cl)

setMethod("t", signature(x = "unpackedMatrix"),
          function(x)
              .Call(unpackedMatrix_transpose, x))
setMethod("diag", signature(x = "unpackedMatrix"),
          function(x, nrow, ncol, names = TRUE)
              .Call(unpackedMatrix_diag_get, x, names))
setMethod("diag<-", signature(x = "unpackedMatrix"),
          function(x, value)
              .Call(unpackedMatrix_diag_set, x, value))

setMethod("symmpart", signature(x = "unpackedMatrix"),
          function(x) .Call(unpackedMatrix_symmpart, x))
setMethod("symmpart", signature(x = "matrix"),
          ## function(x) .Call(matrix_symmpart, x)) # returning .syMatrix
          function(x) 0.5 * symmetrizeDimnames(x + t(x))) # returning matrix

setMethod("skewpart", signature(x = "unpackedMatrix"),
          function(x) .Call(unpackedMatrix_skewpart, x))
setMethod("skewpart", signature(x = "matrix"),
          ## function(x) .Call(matrix_skewpart, x)) # returning .geMatrix
          function(x) 0.5 * symmetrizeDimnames(x - t(x))) # returning matrix

setMethod("cov2cor", signature(V = "unpackedMatrix"),
          function(V) as(forceSymmetric(V), "corMatrix"))

rm(.upM.subclasses)
