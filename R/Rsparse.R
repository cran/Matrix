#### Sparse Matrices in Compressed row-oriented format
####                               --- "R"

### ``mainly for completeness'' --- we *do* favour Csparse
##    - - - - - - - - - - - -   hence only "minimal" methods here !
##  see also ./SparseM-conv.R

### contains = "dMatrix"

.R.2.T <- function(from) .Call(compressed_to_TMatrix, from, FALSE)
.R.2.C <- function(from) .Call(R_to_CMatrix, from)

setAs("RsparseMatrix", "TsparseMatrix", .R.2.T)
setAs("RsparseMatrix", "CsparseMatrix", .R.2.C)
## for printing etc:
setAs("RsparseMatrix", "dgeMatrix",
      function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "dgeMatrix"))
setAs("RsparseMatrix", "matrix",
      function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "matrix"))

setAs("RsparseMatrix", "dsparseMatrix",
      function(from) as(.Call(R_to_CMatrix, from), "dsparseMatrix"))
setAs("RsparseMatrix", "lsparseMatrix",
      function(from) as(.Call(R_to_CMatrix, from), "lsparseMatrix"))
setAs("RsparseMatrix", "nsparseMatrix",
      function(from) as(.Call(R_to_CMatrix, from), "nsparseMatrix"))

setAs("RsparseMatrix", "dMatrix",
      function(from) as(.Call(R_to_CMatrix, from), "dMatrix"))
setAs("RsparseMatrix", "lMatrix",
      function(from) as(.Call(R_to_CMatrix, from), "lMatrix"))
setAs("RsparseMatrix", "nMatrix",
      function(from) as(.Call(R_to_CMatrix, from), "nMatrix"))

##--- and all these are just "the essential low-level coercions" : ----------

## setAs("dgRMatrix", "matrix",
##       function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "matrix"))
## setAs("lgRMatrix", "matrix",
##       function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "matrix"))
## setAs("ngRMatrix", "matrix",
##       function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "matrix"))

setAs("dgRMatrix", "dgeMatrix",
      function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "dgeMatrix"))
setAs("lgRMatrix", "lgeMatrix",
      function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "lgeMatrix"))
setAs("ngRMatrix", "ngeMatrix",
      function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "ngeMatrix"))

setAs("dgRMatrix", "dgCMatrix", .R.2.C)
setAs("lgRMatrix", "lgCMatrix", .R.2.C)
setAs("ngRMatrix", "ngCMatrix", .R.2.C)
## really needed? :
setAs("dgRMatrix", "CsparseMatrix", function(from) as(from, "dgCMatrix"))


setAs("dgRMatrix", "dgTMatrix", .R.2.T)
setAs("lgRMatrix", "lgTMatrix", .R.2.T)
setAs("ngRMatrix", "ngTMatrix", .R.2.T)

##=== Now the same stories for the "s" (symmetric) and "t" (triangular) ones ===

setAs("dsRMatrix", "dsCMatrix", .R.2.C)
setAs("lsRMatrix", "lsCMatrix", .R.2.C)
setAs("nsRMatrix", "nsCMatrix", .R.2.C)

setAs("dsRMatrix", "dsTMatrix", .R.2.T)
setAs("lsRMatrix", "lsTMatrix", .R.2.T)
setAs("nsRMatrix", "nsTMatrix", .R.2.T)

setAs("dsRMatrix", "dsyMatrix",
      function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "dsyMatrix"))
setAs("lsRMatrix", "lsyMatrix",
      function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "lsyMatrix"))
setAs("nsRMatrix", "nsyMatrix",
      function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "nsyMatrix"))

setAs("dtRMatrix", "dtCMatrix", .R.2.C)
setAs("ltRMatrix", "ltCMatrix", .R.2.C)
setAs("ntRMatrix", "ntCMatrix", .R.2.C)

setAs("dtRMatrix", "dtTMatrix", .R.2.T)
setAs("ltRMatrix", "ltTMatrix", .R.2.T)
setAs("ntRMatrix", "ntTMatrix", .R.2.T)

setAs("dtRMatrix", "dtrMatrix",
      function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "dtrMatrix"))
setAs("ltRMatrix", "ltrMatrix",
      function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "ltrMatrix"))
setAs("ntRMatrix", "ntrMatrix",
      function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "ntrMatrix"))

##setAs("matrix", "dgRMatrix",
##      function(from) {
##          storage.mode(from) <- "double"
##          .Call(matrix_to_csc, from)
##      })

## **VERY** cheap substitutes:  work via dgC and t(.)
.viaC.to.dgR <- function(from) {
    m <- as(t(from), "dgCMatrix")
    new("dgRMatrix", Dim = dim(from), Dimnames = .M.DN(from),
	p = m@p, j = m@i, x = m@x)
}

setAs("matrix",    "dgRMatrix", .viaC.to.dgR)
setAs("dgeMatrix", "dgRMatrix", .viaC.to.dgR)
setAs("dgCMatrix", "dgRMatrix", .viaC.to.dgR)
setAs("dgTMatrix", "dgRMatrix", .viaC.to.dgR)

## symmetric: can use same 'p' slot
setAs("dsCMatrix", "dsRMatrix",
      function(from) new("dsRMatrix", Dim = dim(from), Dimnames = .M.DN(from),
	      p = from@p, j = from@i, x = from@x,
	      uplo = if (from@uplo == "U") "L" else "U"))

setAs("dtCMatrix", "dtRMatrix", .viaC.to.dgR) # should work; can NOT use 'p'


##setAs("dgRMatrix", "dgeMatrix",
##      function(from) .Call(csc_to_dgeMatrix, from))

##setAs("matrix", "dgRMatrix",
##      function(from) {
##          storage.mode(from) <- "double"
##          .Call(matrix_to_csc, from)
##      })


##setMethod("diag", signature(x = "dgRMatrix"),
##          function(x = 1, nrow, ncol = n) .Call(csc_getDiag, x))

## try to define for "Matrix" -- once and for all -- but that fails -- why? __ FIXME __
## setMethod("dim", signature(x = "dgRMatrix"),
##           function(x) x@Dim, valueClass = "integer")

##setMethod("t", signature(x = "dgRMatrix"),
##          function(x) .Call(csc_transpose, x),
##          valueClass = "dgRMatrix")

setMethod("image", "dgRMatrix",
          function(x, ...) {
              x <- as(x, "dgTMatrix")
              callGeneric()
          })

setMethod("t", "RsparseMatrix", function(x) as_Rsparse(t(.R.2.T(x))))


## Want tril(), triu(), band() --- just as "indexing" ---
## return a "close" class:
setMethod("tril", "RsparseMatrix",
	  function(x, k = 0, ...) as_Rsparse(tril(.R.2.C(x), k = k, ...)))
setMethod("triu", "RsparseMatrix",
	  function(x, k = 0, ...) as_Rsparse(triu(.R.2.C(x), k = k, ...)))
setMethod("band", "RsparseMatrix",
	  function(x, k1, k2, ...)
	  as_Rsparse(band(.R.2.C(x), k1 = k1, k2 = k2, ...)))
