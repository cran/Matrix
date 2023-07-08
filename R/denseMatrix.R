## METHODS FOR CLASS: denseMatrix (virtual)
## dense matrices with unpacked _or_ packed storage
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## README: Many of these methods are for _subclasses_ of denseMatrix,
##         for base matrices, or for base vectors.  They have been
##         centralized here quite on purpose, for easier maintenance,
##         to prevent infelicities among groups of similar methods,
##         and to avoid accidental gaps in implementation.


## ~~~~ COERCIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

..dense2dsparse <- function(from)
    .Call(R_dense_as_sparse, from, "d.C", NULL, NULL)
..dense2lsparse <- function(from)
    .Call(R_dense_as_sparse, from, "l.C", NULL, NULL)
..dense2nsparse <- function(from)
    .Call(R_dense_as_sparse, from, "n.C", NULL, NULL)
..dense2Csparse <- function(from)
    .Call(R_dense_as_sparse, from, "..C", NULL, NULL)
..dense2Rsparse <- function(from)
    .Call(R_dense_as_sparse, from, "..R", NULL, NULL)
..dense2Tsparse <- function(from)
    .Call(R_dense_as_sparse, from, "..T", NULL, NULL)

..dense2dge <- function(from) .Call(R_dense_as_general, from, "d")
..dense2lge <- function(from) .Call(R_dense_as_general, from, "l")
..dense2nge <- function(from) .Call(R_dense_as_general, from, "n")
..dense2g   <- function(from) .Call(R_dense_as_general, from, ".")

..m2ddense <- function(from) .m2dense.checking(from, "d")
..m2ldense <- function(from) .m2dense.checking(from, "l")
..m2ndense <- function(from) .m2dense.checking(from, "n")
..m2dense  <- function(from) .m2dense.checking(from, ".")

..m2dsparse <- function(from) .m2sparse.checking(from, "d", "C")
..m2lsparse <- function(from) .m2sparse.checking(from, "l", "C")
..m2nsparse <- function(from) .m2sparse.checking(from, "n", "C")
..m2Csparse <- function(from) .m2sparse.checking(from, ".", "C")
..m2Rsparse <- function(from) .m2sparse.checking(from, ".", "R")
..m2Tsparse <- function(from) .m2sparse.checking(from, ".", "T")

..m2dge <- function(from)
    .Call(R_matrix_as_dense, from, "dge", NULL, NULL)
..m2lge <- function(from)
    .Call(R_matrix_as_dense, from, "lge", NULL, NULL)
..m2nge <- function(from)
    .Call(R_matrix_as_dense, from, "nge", NULL, NULL)
..m2ge  <- function(from)
    .Call(R_matrix_as_dense, from, ".ge", NULL, NULL)

..m2dgC <- function(from)
    .Call(R_dense_as_sparse, from, "dgC", NULL, NULL)
..m2lgC <- function(from)
    .Call(R_dense_as_sparse, from, "lgC", NULL, NULL)
..m2ngC <- function(from)
    .Call(R_dense_as_sparse, from, "ngC", NULL, NULL)
..m2gC  <- function(from)
    .Call(R_dense_as_sparse, from, ".gC", NULL, NULL)
..m2gR  <- function(from)
    .Call(R_dense_as_sparse, from, ".gR", NULL, NULL)
..m2gT  <- function(from)
    .Call(R_dense_as_sparse, from, ".gT", NULL, NULL)

..pack   <- function(from)   pack(from)
..unpack <- function(from) unpack(from)

..ge2m  <- function(from) .Call(R_geMatrix_as_matrix, from, FALSE)
..nge2m <- function(from) .Call(R_geMatrix_as_matrix, from,  TRUE)

..ge2v  <- function(from) .Call(R_geMatrix_as_vector, from, FALSE)
..nge2v <- function(from) .Call(R_geMatrix_as_vector, from,  TRUE)


## To denseMatrix ..........................................

setAs("ANY", "denseMatrix",
      function(from) Matrix(from, sparse = FALSE, doDiag = FALSE))

setAs("matrix", "denseMatrix", ..m2dense)

setAs("numLike", "denseMatrix", ..m2ge)

## To sparseMatrix .........................................

setAs("denseMatrix",  "sparseMatrix", ..dense2Csparse)
setAs("denseMatrix", "CsparseMatrix", ..dense2Csparse)
setAs("denseMatrix", "RsparseMatrix", ..dense2Rsparse)
setAs("denseMatrix", "TsparseMatrix", ..dense2Tsparse)

setAs("matrix",  "sparseMatrix", ..m2Csparse)
setAs("matrix", "CsparseMatrix", ..m2Csparse)
setAs("matrix", "RsparseMatrix", ..m2Rsparse)
setAs("matrix", "TsparseMatrix", ..m2Tsparse)

setAs("numLike",  "sparseMatrix", ..m2gC)
setAs("numLike", "CsparseMatrix", ..m2gC)
setAs("numLike", "RsparseMatrix", ..m2gR)
setAs("numLike", "TsparseMatrix", ..m2gT)

## As many people have depended on this being fast:
setAs("matrix", "dgCMatrix", ..m2dgC)

## To base matrix, base vector .............................

setAs("denseMatrix", "matrix", .dense2m)
setAs("denseMatrix", "vector", .dense2v)

setMethod("as.vector", signature(x = "denseMatrix"),
          function(x, mode = "any") as.vector(.dense2v(x), mode))
setMethod("as.numeric", signature(x = "denseMatrix"),
          function(x, ...) as.double(.dense2v(x)))
setMethod("as.logical", signature(x = "denseMatrix"),
          function(x, ...) as.logical(.dense2v(x)))

## Faster:
for (.kind in c("d", "l", "n")) {
    .from <- paste0(.kind, "geMatrix")
    setAs(.from, "matrix", if(.kind != "n") ..ge2m else ..nge2m)
    setAs(.from, "vector", if(.kind != "n") ..ge2v else ..nge2v)

    setMethod("as.vector", signature(x = .from),
              if(.kind != "n")
                  function(x, mode) as.vector(.ge2v(x, FALSE), mode)
              else
                  function(x, mode) as.vector(.ge2v(x,  TRUE), mode))
}
rm(.kind, .from)

## To "kind" ...............................................

setAs("denseMatrix", "dMatrix", ..dense2d)
setAs("denseMatrix", "lMatrix", ..dense2l)
setAs("denseMatrix", "nMatrix", ..dense2n)

setAs("denseMatrix", "ddenseMatrix", ..dense2d)
setAs("denseMatrix", "ldenseMatrix", ..dense2l)
setAs("denseMatrix", "ndenseMatrix", ..dense2n)

setAs("denseMatrix", "dsparseMatrix", ..dense2dsparse)
setAs("denseMatrix", "lsparseMatrix", ..dense2lsparse)
setAs("denseMatrix", "nsparseMatrix", ..dense2nsparse)

setAs("matrix", "dMatrix",
      function(from) {
          storage.mode(from) <- "double"
          if(isDiagonal(from))
              forceDiagonal(from)
          else if(sparseDefault(from))
              .m2sparse.checking(from, "d", "C")
          else .m2dense.checking(from, "d")
      })
setAs("matrix", "lMatrix",
      function(from) {
          storage.mode(from) <- "logical"
          if(isDiagonal(from))
              forceDiagonal(from)
          else if(sparseDefault(from))
              .m2sparse.checking(from, "l", "C")
          else .m2dense.checking(from, "l")
      })
setAs("matrix", "nMatrix",
      function(from) {
          storage.mode(from) <- "logical"
          if(sparseDefault(from))
              .m2sparse.checking(from, "n", "C")
          else .m2dense.checking(from, "n")
      })

setAs("matrix", "ddenseMatrix", ..m2ddense)
setAs("matrix", "ldenseMatrix", ..m2ldense)
setAs("matrix", "ndenseMatrix", ..m2ndense)

setAs("matrix", "dsparseMatrix", ..m2dsparse)
setAs("matrix", "lsparseMatrix", ..m2lsparse)
setAs("matrix", "nsparseMatrix", ..m2nsparse)

setAs("numLike", "dMatrix",
      function(from) {
          if(sparseDefault(from))
              .m2sparse(from, "dgC", NULL, NULL)
          else .m2dense(from, "dge", NULL, NULL)
      })

setAs("numLike", "lMatrix",
      function(from) {
          if(sparseDefault(from))
              .m2sparse(from, "lgC", NULL, NULL)
          else .m2dense(from, "lge", NULL, NULL)
      })

setAs("numLike", "nMatrix",
      function(from) {
          if(sparseDefault(from))
              .m2sparse(from, "ngC", NULL, NULL)
          else .m2dense(from, "nge", NULL, NULL)
      })

setAs("numLike", "ddenseMatrix", ..m2dge)
setAs("numLike", "ldenseMatrix", ..m2lge)
setAs("numLike", "ndenseMatrix", ..m2nge)

setAs("numLike", "dsparseMatrix", ..m2dgC)
setAs("numLike", "lsparseMatrix", ..m2lgC)
setAs("numLike", "nsparseMatrix", ..m2ngC)

## To general ..............................................

setAs("denseMatrix", "generalMatrix", ..dense2g)
setAs(     "matrix", "generalMatrix", ..m2ge)
setAs(    "numLike", "generalMatrix", ..m2ge)

## To symmetric ............................................

## setAs("denseMatrix", "symmetricMatrix", .) # inherited from Matrix
## setAs(     "matrix", "symmetricMatrix", .) # in ./symmetricMatrix.R

## To triangular ...........................................

## setAs("denseMatrix", "triangularMatrix", .) # inherited from Matrix
## setAs(     "matrix", "triangularMatrix", .) # in ./triangularMatrix.R

## To unpacked .............................................

setAs("denseMatrix", "unpackedMatrix", ..unpack)
setAs(     "matrix", "unpackedMatrix", ..m2dense)
setAs(    "numLike", "unpackedMatrix", ..m2ge)

## To packed ...............................................

setAs("denseMatrix",   "packedMatrix", ..pack)
setAs(     "matrix",   "packedMatrix", ..pack)

rm(..dense2dsparse, ..dense2lsparse, ..dense2nsparse,
   ..dense2Csparse, ..dense2Rsparse, ..dense2Tsparse,
   ..dense2dge, ..dense2lge, ..dense2nge, ..dense2g,
   ..m2ddense, ..m2ldense, ..m2ndense, ..m2dense,
   ..m2dsparse, ..m2lsparse, ..m2nsparse,
   ..m2Csparse, ..m2Rsparse, ..m2Tsparse,
   ..m2dge, ..m2lge, ..m2nge, ..m2ge,
   ..m2dgC, ..m2lgC, ..m2ngC, ..m2gC, ..m2gR, ..m2gT,
   ..pack, ..unpack,
   ..ge2m, ..nge2m,
   ..ge2v, ..nge2v)


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dim<-", signature(x = "denseMatrix"),
          function(x, value) {
              if(!is.numeric(value) || length(value) != 2L)
                  stop("dimensions must be numeric of length 2")
              if(anyNA(value))
                  stop("dimensions cannot contain NA")
              if(any(value < 0))
                  stop("dimensions cannot contain negative values")
              if(!is.integer(value)) {
                  if(any(value > .Machine$integer.max))
                      stop("dimensions cannot exceed 2^31-1")
                  value <- as.integer(value)
              }
              if(all(value == (d <- x@Dim)))
                  return(x)
              if((pv <- prod(value)) != (pd <- prod(d)))
                  stop(gettextf("assigned dimensions [product %.0f] do not match Matrix length [%.0f]",
                                pv, pd, domain = NA))
              r <- .dense2g(x)
              r@Dim <- value
              r@factors <- list()
              r
          })

setMethod("mean", signature(x = "denseMatrix"),
          function(x, trim = 0, na.rm = FALSE, ...) {
              if(is.numeric(trim) && length(trim) == 1L && !is.na(trim) &&
                 trim == 0) {
                  ## Be fast in this special case :
                  if(isTRUE(na.rm))
                      x <- x[!is.na(x)]
                  sum(x) / length(x)
              } else mean.default(.dense2v(x), trim = trim, na.rm = na.rm, ...)
          })

setMethod("rep", "denseMatrix",
          function(x, ...) rep(.dense2v(x), ...))

setMethod("show", "denseMatrix",
          function(object) prMatrix(object))

.dense.band <- function(x, k1, k2, ...) .Call(R_dense_band, x, k1, k2)
.dense.triu <- function(x, k = 0L,  ...) .Call(R_dense_band, x, k, NULL)
.dense.tril <- function(x, k = 0L,  ...) .Call(R_dense_band, x, NULL, k)
for (.cl in c("denseMatrix", "matrix")) {
    setMethod("band", signature(x = .cl), .dense.band)
    setMethod("triu", signature(x = .cl), .dense.triu)
    setMethod("tril", signature(x = .cl), .dense.tril)
}
rm(.dense.band, .dense.triu, .dense.tril, .cl)

## x[] <- value :
setReplaceMethod("[", signature(x = "denseMatrix", i = "missing", j = "missing",
                                value = "ANY"),## double/logical/...
                 function (x, value) {
                     x <- .dense2g(x)
                     x@x[] <- value
                     validObject(x)# check if type and lengths above match
                     x
                 })

## FIXME: 1) These are far from efficient
## -----
setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "missing",
                                value = "replValue"),
                 function (x, i, j, ..., value) {
                     r <- as(x, "matrix")
                     ## message("`[<-` with nargs()= ",nargs())
                     if((na <- nargs()) == 3)
                         r[i] <- value
                     else if(na == 4)
                         r[i, ] <- value
                     else stop(gettextf("invalid nargs()= %d", na), domain=NA)
                     .m2ge(r, .M.kind(x))
                 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "missing", j = "index",
                                value = "replValue"),
                 function (x, i, j, ..., value) {
                     r <- as(x, "matrix")
                     r[, j] <- value
                     .m2ge(r, .M.kind(x))
                 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "index",
                                value = "replValue"),
                 function (x, i, j, ..., value) {
                     r <- as(x, "matrix")
                     r[i, j] <- value
                     as_denseClass(r, class(x)) ## was as(r, class(x))
                 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "matrix",  # 2-col.matrix
                                j = "missing", value = "replValue"),
                 function(x, i, j, ..., value) {
                     r <- as(x, "matrix")
                     r[ i ] <- value
                     .m2ge(r, .M.kind(x))
                 })
