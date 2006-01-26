#### Symmetric Sparse Matrices in compressed column-oriented format

### contains = "dgCMatrix"

setAs("dsCMatrix", "dsTMatrix",
      function(from) ## Cholmod:
      .Call("Csparse_to_Tsparse", from, PACKAGE = "Matrix"))

setAs("dsCMatrix", "dgTMatrix", # needed for image()
      function(from) ## pre-Cholmod:
      .Call("dsCMatrix_to_dgTMatrix", from, PACKAGE = "Matrix"))

setAs("dsCMatrix", "dgeMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgeMatrix"))
setAs("dgeMatrix", "dsCMatrix",
      function(from) as(as(from, "dsyMatrix"), "dsTMatrix"))

setAs("dsCMatrix", "matrix",
      function(from) as(as(from, "dgTMatrix"), "matrix"))
setAs("matrix", "dsCMatrix",
      function(from) as(as(from, "dgTMatrix"), "dsCMatrix"))

setAs("dsCMatrix", "lsCMatrix",
      function(from) new("lsCMatrix", i = from@i, p = from@p, uplo = from@uplo,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("dsCMatrix", "dgCMatrix",
      function(from) .Call("sCMatrix_to_gCMatrix", from, PACKAGE = "Matrix"))

if(FALSE) # have 'C' version above
setAs("dsCMatrix", "dsTMatrix",
      function(from)
      new("dsTMatrix", i = from@i,
          j = .Call("Matrix_expand_pointers", from@p, PACKAGE = "Matrix"),
          x = from@x, uplo = from@uplo,
          Dim= from@Dim, Dimnames = from@Dimnames)
      )

setAs("dsCMatrix", "dsyMatrix",
      function(from) as(as(from, "dsTMatrix"), "dsyMatrix"))

setMethod("solve", signature(a = "dsCMatrix", b = "dgeMatrix"),
          function(a, b, ...)
          .Call("dsCMatrix_matrix_solve", a, b, TRUE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dsCMatrix", b = "matrix"),
          function(a, b, ...)
          .Call("dsCMatrix_matrix_solve", a, b, FALSE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

##setMethod("solve", signature(a = "dsCMatrix", b = "numeric"),
##          function(a, b, ...) callGeneric(a, as.matrix(b)),
##          valueClass = "dgeMatrix")

setMethod("chol", signature(x = "dsCMatrix", pivot = "missing"),
          function(x, pivot, LINPACK) .Call("dsCMatrix_chol", x, TRUE, PACKAGE = "Matrix"))

setMethod("chol", signature(x = "dsCMatrix", pivot = "logical"),
          function(x, pivot, LINPACK) .Call("dsCMatrix_chol", x, pivot, PACKAGE = "Matrix"))

setMethod("t", signature(x = "dsCMatrix"),
          function(x) .Call("ssc_transpose", x, PACKAGE = "Matrix"),
          valueClass = "dsCMatrix")

setMethod("determinant", signature(x = "dsCMatrix", logarithm = "missing"),
          function(x, logarithm, ...) determinant(x, TRUE))

setMethod("determinant", signature(x = "dsCMatrix", logarithm = "logical"),
          function(x, logarithm, ...)
      {
          ldet <- sum(log(chol(x)@D))
          modulus <- if (logarithm) ldet else exp(ldet)
          attr(modulus, "logarithm") <- logarithm
          val <- list(modulus = modulus, sign = as.integer(1))
          class(val) <- "det"
          val
      })

setMethod("writeHB", signature(obj = "dsCMatrix"),
          function(obj, file, ...)
          .Call("Matrix_writeHarwellBoeing",
                if (obj@uplo == "U") t(obj) else obj,
                as.character(file), "DSC", PACKAGE = "Matrix"))

setMethod("writeMM", signature(obj = "dsCMatrix"),
          function(obj, file, ...)
          .Call("Matrix_writeMatrixMarket",
                if (obj@uplo == "U") t(obj) else obj,
                as.character(file), "DSC", PACKAGE = "Matrix"))
