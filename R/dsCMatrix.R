#### Symmetric Sparse Matrices in compressed column-oriented format

##setAs("dgCMatrix", "dsCMatrix", ...

## Specific conversions, should they be necessary.  Better to convert as
## as(x, "TsparseMatrix") or as(x, "denseMatrix")

## Moved to ./Csparse.R
## setAs("dsCMatrix", "dsTMatrix",
##       function(from) .Call(Csparse_to_Tsparse, from, FALSE))

setAs("dsCMatrix", "dgTMatrix", # needed for show(), image()
      function(from)
      ## pre-Cholmod -- FIXME: get rid of
      .Call(dsCMatrix_to_dgTMatrix, from))

setAs("dsCMatrix", "dgeMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgeMatrix"))

setAs("dsCMatrix", "matrix",
      function(from) as(as(from, "generalMatrix"), "matrix"))
setAs("matrix", "dsCMatrix",
      function(from) as(as(from, "CsparseMatrix"), "symmetricMatrix"))

setAs("dsCMatrix", "lsCMatrix",
      function(from) new("lsCMatrix", i = from@i, p = from@p, uplo = from@uplo,
                         x = as.logical(from@x),
                         Dim = from@Dim, Dimnames = from@Dimnames))
setAs("dsCMatrix", "nsCMatrix",
      function(from) new("nsCMatrix", i = from@i, p = from@p, uplo = from@uplo,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("dsCMatrix", "dgCMatrix",
      function(from) .Call(Csparse_symmetric_to_general, from))

setAs("dsCMatrix", "dsyMatrix",
      function(from) as(from, "denseMatrix"))

## have rather tril() and triu() methods than
## setAs("dsCMatrix", "dtCMatrix", ....)
setMethod("tril", "dsCMatrix",
	  function(x, k = 0, ...) {
	      if(x@uplo == "L" && k == 0)
		  ## same internal structure (speedup potential !?)
		  new("dtCMatrix", uplo = x@uplo, i = x@i, p = x@p,
		      x = x@x, Dim = x@Dim, Dimnames = x@Dimnames)
	      else tril(as(x, "dgCMatrix"), k = k, ...)
	  })

setMethod("triu", "dsCMatrix",
	  function(x, k = 0, ...) {
	      if(x@uplo == "U" && k == 0)
		  ## same internal structure (speedup potential !?)
		  new("dtCMatrix", uplo = x@uplo, i = x@i, p = x@p,
		      x = x@x, Dim = x@Dim, Dimnames = x@Dimnames)
	      else triu(as(x, "dgCMatrix"), k = k, ...)
	  })

setMethod("solve", signature(a = "dsCMatrix", b = "ddenseMatrix"),
          function(a, b, ...) {
              if (class(b) != "dgeMatrix")
                  b <- .Call(dup_mMatrix_as_dgeMatrix, b)
              .Call(dsCMatrix_matrix_solve, a, b)
          },
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dsCMatrix", b = "matrix"),
          function(a, b, ...)
          .Call(dsCMatrix_matrix_solve, a,
                .Call(dup_mMatrix_as_dgeMatrix, b)),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dsCMatrix", b = "numeric"),
          function(a, b, ...)
          .Call(dsCMatrix_matrix_solve, a,
                .Call(dup_mMatrix_as_dgeMatrix, b)),
          valueClass = "dgeMatrix")

## `` Fully-sparse'' solve() :
setMethod("solve", signature(a = "dsCMatrix", b = "dsparseMatrix"),
	  function(a, b, ...) {
	      if (!is(b, "CsparseMatrix"))
		  b <- as(b, "CsparseMatrix")
	      if (is(b, "symmetricMatrix")) ## not supported (yet) by cholmod_spsolve
		  b <- as(b, "dgCMatrix")
	      .Call(dsCMatrix_Csparse_solve, a, b)
	  })


setMethod("chol", signature(x = "dsCMatrix", pivot = "missing"),
	  function(x, pivot, ...) .Call(dsCMatrix_chol, x, FALSE),
	  valueClass = "dtCMatrix")

setMethod("chol", signature(x = "dsCMatrix", pivot = "logical"),
	  function(x, pivot, ...) .Call(dsCMatrix_chol, x, pivot),
	  valueClass = "dtCMatrix")

setMethod("Cholesky", signature(A = "dsCMatrix"),
          function(A, perm = TRUE, LDL = TRUE, super = FALSE, ...)
          .Call(dsCMatrix_Cholesky, A, perm, LDL, super))


setMethod("t", signature(x = "dsCMatrix"),
          function(x) .Call(Csparse_transpose, x, FALSE),
          valueClass = "dsCMatrix")

setMethod("determinant", signature(x = "dsCMatrix", logarithm = "missing"),
          function(x, logarithm, ...) determinant(x, TRUE))

setMethod("determinant", signature(x = "dsCMatrix", logarithm = "logical"),
          function(x, logarithm, ...)
      {
          stop("Temporarily disabled until we work out the LDL factorization diagonal")
          ldet <- sum(log(chol(x)@D))
          modulus <- if (logarithm) ldet else exp(ldet)
          attr(modulus, "logarithm") <- logarithm
          val <- list(modulus = modulus, sign = as.integer(1))
          class(val) <- "det"
          val
      })

## setMethod("writeHB", signature(obj = "dsCMatrix"),
##           function(obj, file, ...) {
##               .Deprecated("writeMM")
##               .Call(Matrix_writeHarwellBoeing,
##                     if (obj@uplo == "U") t(obj) else obj,
##                     as.character(file), "DSC")
##           })
