## METHODS FOR CLASS: sparseMatrix, [CRT]sparseMatrix (virtual)
## sparse matrices, in some cases restricted to CSC, CSR, triplet
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.sparse.band <- function(x, k1, k2, ...)
    .Call(R_sparse_band, x, k1, k2)
.sparse.triu <- function(x, k = 0L, ...)
    .Call(R_sparse_band, x, k, NULL)
.sparse.tril <- function(x, k = 0L, ...)
    .Call(R_sparse_band, x, NULL, k)
.sparse.diag.get <- function(x, nrow, ncol, names = TRUE)
    .Call(R_sparse_diag_get, x, names)
.sparse.diag.set <- function(x, value)
    .Call(R_sparse_diag_set, x, value)
.sparse.t <- function(x)
    .Call(R_sparse_transpose, x, FALSE)
.sparse.fS1  <- function(x, uplo)
    .Call(R_sparse_force_symmetric, x, NULL)
.sparse.fS2  <- function(x, uplo)
    .Call(R_sparse_force_symmetric, x, uplo)
.sparse.symmpart <- function(x)
    .Call(R_sparse_symmpart, x)
.sparse.skewpart <- function(x)
    .Call(R_sparse_skewpart, x)
.sparse.is.di <- function(object)
    .Call(R_sparse_is_diagonal, object)
.sparse.is.tr <- function(object, upper = NA, ...)
    .Call(R_sparse_is_triangular, object, upper)
.sparse.is.sy <- function(object, checkDN = TRUE, ...) {
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    .Call(R_sparse_is_symmetric, object, checkDN)
}
.sparse.is.sy.dz <- function(object, checkDN = TRUE,
                             tol = 100 * .Machine$double.eps, ...) {
    ## backwards compatibility: don't check DN if check.attributes=FALSE
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    ## be very fast when requiring exact symmetry
    if(tol <= 0)
        return(.Call(R_sparse_is_symmetric, object, checkDN))
    ## pretest: is it square?
    d <- object@Dim
    if((n <- d[2L]) != d[1L])
        return(FALSE)
    ## pretest: are DN symmetric in the sense of validObject(<symmetricMatrix>)?
    if(checkDN && !isSymmetricDN(object@Dimnames))
        return(FALSE)
    if(n == 0L)
        return(TRUE)

    ## now handling an n-by-n [dz]g[CRT]Matrix, n >= 1:

    Cj <- if(is.complex(object@x)) Conj else identity
    ae <- function(check.attributes, ...) {
        ## discarding possible user-supplied check.attributes
        all.equal(..., check.attributes = FALSE)
    }

    isTRUE(ae(target    = .M2V(     object  ),
              current   = .M2V(Cj(t(object))),
              tolerance = tol, ...))
}

setMethod("diff", signature(x = "sparseMatrix"),
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

setMethod("mean", signature(x = "sparseMatrix"),
          function(x, ...) mean(as(x, "sparseVector"), ...))

setMethod("rep", "sparseMatrix",
          function(x, ...)  rep(as(x, "sparseVector"), ...))

for(.cl in paste0(c("C", "R", "T"), "sparseMatrix")) {
setMethod("band"  , signature(x = .cl), .sparse.band)
setMethod("triu"  , signature(x = .cl), .sparse.triu)
setMethod("tril"  , signature(x = .cl), .sparse.tril)
setMethod("diag"  , signature(x = .cl), .sparse.diag.get)
setMethod("diag<-", signature(x = .cl), .sparse.diag.set)
setMethod("t"     , signature(x = .cl), .sparse.t)
setMethod("forceSymmetric", signature(x = .cl, uplo =   "missing"), .sparse.fS1)
setMethod("forceSymmetric", signature(x = .cl, uplo = "character"), .sparse.fS2)
setMethod("symmpart", signature(x = .cl), .sparse.symmpart)
setMethod("skewpart", signature(x = .cl), .sparse.skewpart)
setMethod("isSymmetric" , signature(object = .cl), .sparse.is.sy)
setMethod("isTriangular", signature(object = .cl), .sparse.is.tr)
setMethod("isDiagonal"  , signature(object = .cl), .sparse.is.di)
}

.sparse.subclasses <- names(getClassDef("sparseMatrix")@subclasses)
for(.cl in grep("^[dz][gt][CRT]Matrix$", .sparse.subclasses, value = TRUE))
setMethod("isSymmetric" , signature(object = .cl), .sparse.is.sy.dz)
rm(.cl, .sparse.subclasses)

rm(list = c(grep("^[.]sparse[.](band|tri[ul]|diag[.](get|set)|t|fS[12]|symmpart|skewpart|is[.](sy|tr|di)([.]dz)?)$",
                 ls(all.names = TRUE), value = TRUE)))


if(FALSE) ### FIXME: This would *NOT* be needed, if    as.matrix(<sparseMatrix>) was a no-op ;
         ### -----  and then,  base::scale() -> base::scale.default() would work "magically" already..
scale.sparseMatrix <- function(x, center = FALSE, scale = TRUE) {
    if(center) warning("a sparseMatrix should rarely be centered: will not be sparse anymore")
    ## x <- as.matrix(x)

    ## This rest is *identically*  == base :: scale.default :
    nc <- ncol(x)
    if (is.logical(center)) {
        if (center) {
            center <- colMeans(x, na.rm=TRUE)
            x <- sweep(x, 2L, center, check.margin=FALSE)
        }
    }
    else if (is.numeric(center) && (length(center) == nc))
        x <- sweep(x, 2L, center, check.margin=FALSE)
    else
        stop("length of 'center' must equal the number of columns of 'x'")
    if (is.logical(scale)) {
        if (scale) {
            f <- function(v) {
                v <- v[!is.na(v)]
                sqrt(sum(v^2) / max(1, length(v) - 1L))
            }
            scale <- apply(x, 2L, f)
            x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
        }
    }
    else if (is.numeric(scale) && length(scale) == nc)
        x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
    else
        stop("length of 'scale' must equal the number of columns of 'x'")
    if(is.numeric(center)) attr(x, "scaled:center") <- center
    if(is.numeric(scale)) attr(x, "scaled:scale") <- scale
    x
}
