#### Collect methods for  colSums(), rowSums(), colMeans(), rowMeans() here.
####			  =======    -------    --------    --------

### Dense Matrices: -------------------------------------------------

setMethod("colSums", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call(dgeMatrix_colsums, x, na.rm, TRUE, FALSE))

setMethod("colMeans", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call(dgeMatrix_colsums, x, na.rm, TRUE, TRUE))

setMethod("rowSums", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call(dgeMatrix_colsums, x, na.rm, FALSE, FALSE))

setMethod("rowMeans", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call(dgeMatrix_colsums, x, na.rm, FALSE, TRUE))

## FIXME: "works" but not optimally for triangular/symmetric (esp. packed)
.recall.as.dge <- function(x, na.rm = FALSE, dims = 1) {
    x <- .dense2g(x, "d")
    callGeneric()
}
setMethod("colSums",  signature(x = "denseMatrix"), .recall.as.dge)
setMethod("colMeans", signature(x = "denseMatrix"), .recall.as.dge)
setMethod("rowSums",  signature(x = "denseMatrix"), .recall.as.dge)
setMethod("rowMeans", signature(x = "denseMatrix"), .recall.as.dge)
rm(.recall.as.dge)


### Sparse Matrices: -------------------------------------------------

## Diagonal ones:
.diag.Sum <- function(x, na.rm = FALSE, dims = 1)
    if(x@diag == "U") rep(1, x@Dim[1]) else as.numeric(x@x)
.diag.Mean <- function(x, na.rm = FALSE, dims = 1) {
    n <- x@Dim[1L]
    if(x@diag == "U") rep(1/n, n) else as.numeric(x@x)/n
}

setMethod("colSums",  signature(x = "diagonalMatrix"), .diag.Sum)
setMethod("rowSums",  signature(x = "diagonalMatrix"), .diag.Sum)
setMethod("colMeans", signature(x = "diagonalMatrix"), .diag.Mean)
setMethod("rowMeans", signature(x = "diagonalMatrix"), .diag.Mean)

rm(.diag.Sum, .diag.Mean)

### Csparse --- the fast workhorse ones

### 1) those with .Call(.), {d, i, l, n} gCMatrix  x  {col|row}{Sums|Means} :

## the last two arguments to .gCMatrix_(col|col)(Sums|Means)  are 'trans' and 'means'
setMethod("colSums", signature(x = "dgCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(dgCMatrix_colSums, x, na.rm, sparseResult, FALSE, FALSE))

setMethod("rowSums", signature(x = "dgCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(dgCMatrix_colSums, x, na.rm, sparseResult, TRUE, FALSE))

setMethod("colMeans", signature(x = "dgCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(dgCMatrix_colSums, x, na.rm, sparseResult, FALSE, TRUE))

setMethod("rowMeans", signature(x = "dgCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(dgCMatrix_colSums, x, na.rm, sparseResult, TRUE, TRUE))

## not yet
if(FALSE) {
setMethod("colSums", signature(x = "igCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(igCMatrix_colSums, x, na.rm, sparseResult, FALSE, FALSE))

setMethod("rowSums", signature(x = "igCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(igCMatrix_colSums, x, na.rm, sparseResult, TRUE, FALSE))

setMethod("colMeans", signature(x = "igCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(igCMatrix_colSums, x, na.rm, sparseResult, FALSE, TRUE))

setMethod("rowMeans", signature(x = "igCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(igCMatrix_colSums, x, na.rm, sparseResult, TRUE, TRUE))
}

setMethod("colSums", signature(x = "lgCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(lgCMatrix_colSums, x, na.rm, sparseResult, FALSE, FALSE))

setMethod("rowSums", signature(x = "lgCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(lgCMatrix_colSums, x, na.rm, sparseResult, TRUE, FALSE))

setMethod("colMeans", signature(x = "lgCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(lgCMatrix_colSums, x, na.rm, sparseResult, FALSE, TRUE))

setMethod("rowMeans", signature(x = "lgCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(lgCMatrix_colSums, x, na.rm, sparseResult, TRUE, TRUE))

setMethod("colSums", signature(x = "ngCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(ngCMatrix_colSums, x, na.rm, sparseResult, FALSE, FALSE))

setMethod("rowSums", signature(x = "ngCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(ngCMatrix_colSums, x, na.rm, sparseResult, TRUE, FALSE))

setMethod("colMeans", signature(x = "ngCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(ngCMatrix_colSums, x, na.rm, sparseResult, FALSE, TRUE))

setMethod("rowMeans", signature(x = "ngCMatrix"),
	  function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
          .Call(ngCMatrix_colSums, x, na.rm, sparseResult, TRUE, TRUE))

### 2) the other Csparse ones are "just" coerced to a *gCMatrix :
.recall.as.g <- function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE) {
    x <- .sparse2g(x)
    callGeneric()
}
setMethod("colSums",  signature(x = "CsparseMatrix"), .recall.as.g)
setMethod("colMeans", signature(x = "CsparseMatrix"), .recall.as.g)
setMethod("rowSums",  signature(x = "CsparseMatrix"), .recall.as.g)
setMethod("rowMeans", signature(x = "CsparseMatrix"), .recall.as.g)
rm(.recall.as.g)

## --- Tsparse ----

## .as.C.Fun -- since there's now  C code for dgCMatrix_colSums
##     Note: in the past, these went (quite inefficiently)
##           via dgTMatrix, using sparsapply() in ./Auxiliaries.R
.recall.as.C <- function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE) {
    x <- .T2C(x)
    callGeneric()
}
setMethod("colSums",  signature(x = "TsparseMatrix"), .recall.as.C)
setMethod("colMeans", signature(x = "TsparseMatrix"), .recall.as.C)
setMethod("rowSums",  signature(x = "TsparseMatrix"), .recall.as.C)
setMethod("rowMeans", signature(x = "TsparseMatrix"), .recall.as.C)
rm(.recall.as.C)

## --- Rsparse ----

## row <-> col of the "transposed, seen as C" :
setMethod("colSums", signature(x = "RsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
              rowSums(.tCR2RC(x), na.rm = na.rm, dims = dims,
                      sparseResult = sparseResult))
setMethod("colMeans", signature(x = "RsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
              rowMeans(.tCR2RC(x), na.rm = na.rm, dims = dims,
                       sparseResult = sparseResult))
setMethod("rowSums", signature(x = "RsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
              colSums(.tCR2RC(x), na.rm = na.rm, dims = dims,
                      sparseResult = sparseResult))
setMethod("rowMeans", signature(x = "RsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
              colMeans(.tCR2RC(x), na.rm = na.rm, dims = dims,
                       sparseResult = sparseResult))

## --- indMatrix [incl pMatrix ] ---

setMethod("colSums",  signature(x = "indMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
	  tabulate(x@perm, nbins = x@Dim[2L]))
setMethod("colMeans",  signature(x = "indMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
	  tabulate(x@perm, nbins = x@Dim[2L]) / x@Dim[1L])
## for completeness:
setMethod("rowSums",  signature(x = "indMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          rep.int(1, x@Dim[1L]))
setMethod("rowMeans",  signature(x = "indMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          rep.int(1 / x@Dim[2L], x@Dim[1L]))
