## METHODS FOR GENERIC: qr
## QR factorization of dense and sparse matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## FIXME? We could have methods for generalMatrix, symmetricMatrix,
##        triangularMatrix, and diagonalMatrix instead?  We could
##        construct the result "directly" in the upper triangular
##        (incl. diagonal) cases, and cache it in the general and
##        symmetric cases.  See the methods for chol() in ./chol.R ...

## base::qr.default(x) does x <- as.matrix(x) and storage.mode(x) <- "double",
## hence if we deprecate as.matrix(<Matrix>), then we would need this:
if(.Matrix.avoiding.as.matrix) {
setMethod("qr", signature(x = "denseMatrix"),
	  function(x, ...)
              qr.default(as(x, "matrix"), ...))
}

setMethod("qr", signature(x = "sparseMatrix"),
	  function(x, ...)
              qr(..sparse2d(.sparse2g(as(x, "CsparseMatrix"))), ...))

setMethod("qr", signature(x = "dgCMatrix"),
          function(x,
                   ## also had
                   ## > tol = 1e-07, LAPACK = FALSE,
                   ## from base::qr.default() but these are unused
                   ## _and_ not needed for generic consistency ...
                   keep.dimnames = TRUE,
                   verbose = getOption("Matrix.verbose", FALSE), ...)
              .Call(dgCMatrix_QR, x, if(verbose) -1L else 1L, keep.dimnames))
