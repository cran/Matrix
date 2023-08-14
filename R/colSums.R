## METHODS FOR GENERIC: colSums, rowSums, colMeans, rowMeans
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## FIXME? *Sums(<logical matrix>) is currently always of type "double";
##        should *Sums(<([nl]|ind)Matrix>) behave the same?  We are not
##        consistent.  Currently:
##
##        double result:           integer result:
##        * [nl]denseMatrix        * [nl]sparseMatrix
##        * ldiMatrix              * indMatrix
##
##        hence we might consider changing to always give double ...


## ==== denseMatrix ====================================================

setMethod("colSums",  signature(x = "denseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...)
              .Call(R_dense_colSums, x, na.rm, FALSE))
setMethod("colMeans", signature(x = "denseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...)
              .Call(R_dense_colSums, x, na.rm, TRUE))
setMethod("rowSums",  signature(x = "denseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...)
              .Call(R_dense_rowSums, x, na.rm, FALSE))
setMethod("rowMeans", signature(x = "denseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...)
              .Call(R_dense_rowSums, x, na.rm, TRUE))


## ==== sparseMatrix ===================================================

## ---- diagonalMatrix -------------------------------------------------

.diag.cS <- .diag.rS <- function(x, na.rm = FALSE, dims = 1L, ...) {
    if((n <- x@Dim[1L]) == 0L)
        return(double(0L))
    else if(x@diag != "N")
        r <- rep.int(1, n)
    else {
        r <- as.double(x@x)
        if(na.rm)
            r[is.na(r)] <- 0
    }
    if(!is.null(nms <- x@Dimnames[[.MARGIN]]))
        names(r) <- nms
    r
}
body(.diag.cS) <- do.call(substitute, list(body(.diag.cS), list(.MARGIN = 2L)))
body(.diag.rS) <- do.call(substitute, list(body(.diag.rS), list(.MARGIN = 1L)))

.diag.cM <- .diag.rM <- function(x, na.rm = FALSE, dims = 1L, ...) {
    if((n <- x@Dim[1L]) == 0L)
        return(double(0L))
    else if(x@diag != "N")
        r <- rep.int(1 / n, n)
    else {
        r <- as.double(x@x) / n
        if(na.rm)
            r[is.na(r)] <- if(n == 1L) NaN else 0
    }
    if(!is.null(nms <- x@Dimnames[[.MARGIN]]))
        names(r) <- nms
    r
}
body(.diag.cM) <- do.call(substitute, list(body(.diag.cM), list(.MARGIN = 2L)))
body(.diag.rM) <- do.call(substitute, list(body(.diag.rM), list(.MARGIN = 1L)))

setMethod("colSums",  signature(x = "diagonalMatrix"), .diag.cS)
setMethod("colMeans", signature(x = "diagonalMatrix"), .diag.cM)
setMethod("rowSums",  signature(x = "diagonalMatrix"), .diag.rS)
setMethod("rowMeans", signature(x = "diagonalMatrix"), .diag.rM)

rm(.diag.cS, .diag.cM, .diag.rS, .diag.rM)


## ---- indMatrix (incl. pMatrix) --------------------------------------

setMethod("colSums",  signature(x = "indMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...) {
              n <- x@Dim[2L]
              r <- if(x@margin == 1L)
                       tabulate(x@perm, n)
                   else rep.int(1L, n)
              if(!is.null(nms <- x@Dimnames[[2L]]))
                  names(r) <- nms
              r
          })
setMethod("colMeans",  signature(x = "indMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...) {
              n <- (d <- x@Dim)[2L]
              r <- if(x@margin == 1L)
                       tabulate(x@perm, n) / d[1L]
                   else rep.int(1 / d[1L], n)
              if(!is.null(nms <- x@Dimnames[[2L]]))
                  names(r) <- nms
              r
          })
setMethod("rowSums",  signature(x = "indMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...) {
              m <- x@Dim[1L]
              r <- if(x@margin == 1L)
                       rep.int(1L, m)
                   else tabulate(x@perm, m)
              if(!is.null(nms <- x@Dimnames[[1L]]))
                  names(r) <- nms
              r
          })
setMethod("rowMeans",  signature(x = "indMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...) {
              m <- (d <- x@Dim)[1L]
              r <- if(x@margin == 1L)
                       rep.int(1 / d[2L], m)
                   else tabulate(x@perm, m) / d[2L]
              if(!is.null(nms <- x@Dimnames[[1L]]))
                  names(r) <- nms
              r
          })


## ---- CsparseMatrix --------------------------------------------------

setMethod("colSums",  signature(x = "CsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(CRsparse_colSums, x, na.rm, FALSE, sparseResult))
setMethod("colMeans", signature(x = "CsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(CRsparse_colSums, x, na.rm, TRUE, sparseResult))
setMethod("rowSums",  signature(x = "CsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(CRsparse_rowSums, x, na.rm, FALSE, sparseResult))
setMethod("rowMeans", signature(x = "CsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(CRsparse_rowSums, x, na.rm, TRUE, sparseResult))


## ---- RsparseMatrix --------------------------------------------------

setMethod("colSums",  signature(x = "RsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(CRsparse_rowSums, x, na.rm, FALSE, sparseResult))
setMethod("colMeans", signature(x = "RsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(CRsparse_rowSums, x, na.rm, TRUE, sparseResult))
setMethod("rowSums",  signature(x = "RsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(CRsparse_colSums, x, na.rm, FALSE, sparseResult))
setMethod("rowMeans", signature(x = "RsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(CRsparse_colSums, x, na.rm, TRUE, sparseResult))


## ---- TsparseMatrix --------------------------------------------------

setMethod("colSums",  signature(x = "TsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(CRsparse_colSums, .M2C(x), na.rm, FALSE, sparseResult))
setMethod("colMeans",  signature(x = "TsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(CRsparse_colSums, .M2C(x), na.rm, TRUE, sparseResult))
setMethod("rowSums",  signature(x = "TsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(CRsparse_colSums, .M2R(x), na.rm, FALSE, sparseResult))
setMethod("rowMeans",  signature(x = "TsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(CRsparse_colSums, .M2R(x), na.rm, TRUE, sparseResult))
