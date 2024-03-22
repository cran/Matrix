## METHODS FOR CLASS: denseMatrix (virtual)
## dense matrices with unpacked _or_ packed storage
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.dense.band <- function(x, k1, k2, ...)
    .Call(R_dense_band, x, k1, k2)
.dense.triu <- function(x, k = 0L, ...)
    .Call(R_dense_band, x, k, NULL)
.dense.tril <- function(x, k = 0L, ...)
    .Call(R_dense_band, x, NULL, k)
.dense.diag.get <- function(x = 1, nrow, ncol, names = TRUE)
    .Call(R_dense_diag_get, x, names)
.dense.diag.set <- function(x, value)
    .Call(R_dense_diag_set, x, value)
.dense.t <- function(x)
    .Call(R_dense_transpose, x)
.dense.fS1  <- function(x, uplo)
    .Call(R_dense_force_symmetric, x, NULL)
.dense.fS2  <- function(x, uplo)
    .Call(R_dense_force_symmetric, x, uplo)
.dense.symmpart <- function(x)
    .Call(R_dense_symmpart, x)
.dense.skewpart <- function(x)
    .Call(R_dense_skewpart, x)
.dense.is.di <- function(object)
    .Call(R_dense_is_diagonal, object)
.dense.is.tr <- function(object, upper = NA, ...)
    .Call(R_dense_is_triangular, object, upper)
.dense.is.sy <- function(object, checkDN = TRUE, ...) {
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    .Call(R_dense_is_symmetric, object, checkDN)
}
.dense.is.sy.dz <- function(object, checkDN = TRUE,
                            tol = 100 * .Machine$double.eps,
                            tol1 = 8 * tol, ...) {
    ## backwards compatibility: don't check DN if check.attributes=FALSE
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    ## be very fast when requiring exact symmetry
    if(tol <= 0)
        return(.Call(R_dense_is_symmetric, object, checkDN))
    ## pretest: is it square?
    d <- object@Dim
    if((n <- d[2L]) != d[1L])
        return(FALSE)
    ## pretest: are DN symmetric in the sense of validObject(<symmetricMatrix>)?
    if(checkDN && !isSymmetricDN(object@Dimnames))
        return(FALSE)
    if(n == 0L)
        return(TRUE)
    object <- .M2gen(object)

    ## now handling n-by-n [dz]geMatrix, n >= 1:

    Cj <- if(is.complex(object@x)) Conj else identity
    ae <- function(check.attributes, ...) {
        ## discarding possible user-supplied check.attributes
        all.equal.numeric(..., check.attributes = FALSE)
    }

    ## pretest: outermost rows ~= outermost columns?
    ## (fast for large asymmetric)
    if(length(tol1)) {
        i. <- if(n <= 4L) 1L:n else c(1L, 2L, n - 1L, n)
        for(i in i.)
            if(!isTRUE(ae(target = object[i, ], current = Cj(object[, i]),
                          tolerance = tol1, ...)))
                return(FALSE)
    }
    isTRUE(ae(target    =      object  @x,
              current   = Cj(t(object))@x,
              tolerance = tol, ...))
}

setMethod("diff", c(x = "denseMatrix"),
          ## Mostly cut and paste of base::diff.default :
          function(x, lag = 1L, differences = 1L, ...) {
              if(length(lag) != 1L || length(differences) != 1L ||
                  lag < 1L || differences < 1L)
                  stop(gettextf("'%s' and '%s' must be positive integers",
                                "lag", "differences"),
                       domain = NA)
              if(lag * differences >= x@Dim[1L])
                  return(x[0L])
              i1 <- -seq_len(lag)
              for(i in seq_len(differences)) {
                  m <- x@Dim[1L]
                  x <- x[i1, , drop = FALSE] -
                      x[-m:-(m - lag + 1L), , drop = FALSE]
              }
              x
          })

setMethod("mean", c(x = "denseMatrix"),
          function(x, ...) mean.default(.M2v(x), ...))

setMethod("rep", c(x = "denseMatrix"),
          function(x, ...)          rep(.M2v(x), ...))

setMethod("band"  , c(x = "denseMatrix"), .dense.band)

setMethod("triu"  , c(x = "denseMatrix"), .dense.triu)

setMethod("tril"  , c(x = "denseMatrix"), .dense.tril)

setMethod("diag"  , c(x = "denseMatrix"), .dense.diag.get)

setMethod("diag<-", c(x = "denseMatrix"), .dense.diag.set)

setMethod("t"     , c(x = "denseMatrix"), .dense.t)

setMethod("forceSymmetric", c(x = "denseMatrix", uplo =   "missing"), .dense.fS1)

setMethod("forceSymmetric", c(x = "denseMatrix", uplo = "character"), .dense.fS2)

setMethod("symmpart", c(x = "denseMatrix"), .dense.symmpart)

setMethod("skewpart", c(x = "denseMatrix"), .dense.skewpart)

setMethod("isSymmetric" , c(object = "denseMatrix"), .dense.is.sy)

setMethod("isTriangular", c(object = "denseMatrix"), .dense.is.tr)

setMethod("isDiagonal"  , c(object = "denseMatrix"), .dense.is.di)

.dense.subclasses <- names(getClassDef("denseMatrix")@subclasses)
for (.cl in grep("^[dz](ge|tr|tp)Matrix$", .dense.subclasses, value = TRUE))
setMethod("isSymmetric" , c(object = .cl), .dense.is.sy.dz)
rm(.cl, .dense.subclasses)


## METHODS FOR CLASS: unpackedMatrix (virtual)
## dense matrices with unpacked storage
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("unpack", c(x = "packedMatrix"),
          function(x, ...) .Call(R_dense_as_unpacked, x))

setMethod("pack", c(x = "packedMatrix"),
          function(x, ...) x)


## METHODS FOR CLASS: packedMatrix (virtual)
## dense matrices with packed storage
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.uM.pack <-
function(x, ...) .Call(R_dense_as_packed, x, NULL, NULL)

.uM.pack.ge <-
function(x, symmetric = NA, upperTri = NA, ...) {
    if(((sna <- is.na(symmetric)) || symmetric) && isSymmetric(x, ...))
        .Call(R_dense_as_packed, x, "U", NULL)
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

setMethod("unpack", c(x = "unpackedMatrix"),
          function(x, ...) x)

setMethod("pack", c(x = "unpackedMatrix"), .uM.pack)

.uM.subclasses <- names(getClassDef("unpackedMatrix")@subclasses)
for(.cl in grep("^.geMatrix$", .uM.subclasses, value = TRUE))
setMethod("pack", c(x = .cl), .uM.pack.ge)
rm(.cl, .uM.subclasses)


## METHODS FOR CLASS: matrix
## traditional matrices, which really are "dense"
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.m.pack <- .uM.pack.ge
body(.m.pack)[[2L]][[3L]]             <-  quote(.m2dense(x, ".sp",  "U"))
body(.m.pack)[[2L]][[4L]][[3L]][[3L]] <-  quote(.m2dense(x, ".tp", uplo))

setMethod("unpack", c(x = "matrix"),
          function(x, ...) .m2dense.checking(x, "."))
setMethod("pack", c(x = "matrix"), .m.pack)
setMethod("band", c(x = "matrix"), .dense.band)
setMethod("triu", c(x = "matrix"), .dense.triu)
setMethod("tril", c(x = "matrix"), .dense.tril)
setMethod("forceSymmetric", c(x = "matrix", uplo = "missing"),
          function(x, uplo) .m2dense(x, ".sy",  "U"))
setMethod("forceSymmetric", c(x = "matrix", uplo = "character"),
          function(x, uplo) .m2dense(x, ".sy", uplo))
setMethod("symmpart", c(x = "matrix"),
          function(x) symmetrizeDN(0.5 * (x + t(x))))
setMethod("skewpart", c(x = "matrix"),
          function(x) symmetrizeDN(0.5 * (x - t(x))))
setMethod("isTriangular", c(object = "matrix"), .dense.is.tr)
setMethod("isDiagonal"  , c(object = "matrix"), .dense.is.di)

rm(.uM.pack, .uM.pack.ge, .m.pack,
   list = c(grep("^[.]dense[.](band|tri[ul]|diag[.](get|set)|t|fS[12]|symmpart|skewpart|is[.](sy|tr|di)([.]dz)?)$",
                 ls(all.names = TRUE), value = TRUE)))
