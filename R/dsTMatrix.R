  ### Coercion and Methods for Symmetric Triplet Matrices

setAs("dsTMatrix", "dsCMatrix",
      function(from) .Call("dsTMatrix_as_dsCMatrix", from))

## Conversion to dense storage is first to a dsyMatrix
setAs("dsTMatrix", "dsyMatrix",
      function(from) .Call("dsTMatrix_as_dsyMatrix", from))

setAs("dsTMatrix", "dgeMatrix",
      function(from) as(as(from, "dsyMatrix"), "dgeMatrix"))

setAs("dsTMatrix", "matrix",
      function(from) as(as(from, "dsyMatrix"), "matrix"))

setMethod("t", signature(x = "dsTMatrix"),
          function(x)
          new("dsTMatrix", Dim = x@Dim,
              i = x@j, j = x@i, x = x@x,
              uplo = if (x@uplo == "U") "L" else "U"),
          valueClass = "dsTMatrix")

setMethod("writeHB", signature(obj = "dsTMatrix"),
          function(obj, file, ...)
          .Call("Matrix_writeHarwellBoeing",
                if (obj@uplo == "U") t(obj) else obj,
                as.character(file), "DST"))

setMethod("writeMM", signature(obj = "dsTMatrix"),
          function(obj, file, ...)
          .Call("Matrix_writeMatrixMarket",
                if (obj@uplo == "U") t(obj) else obj,
                as.character(file), "DST"))
