#### Sparse Matrices in Compressed column-oriented format

### contains = "dsparseMatrix"

setAs("dgCMatrix", "dgTMatrix",
      function(from) .Call("compressed_to_dgTMatrix", from, TRUE, PACKAGE = "Matrix"))

setAs("dgCMatrix", "matrix",
      function(from) .Call("csc_to_matrix", from, PACKAGE = "Matrix"))

setAs("dgCMatrix", "dgeMatrix",
      function(from) .Call("csc_to_dgeMatrix", from, PACKAGE = "Matrix"))

setAs("dgCMatrix", "lgCMatrix",
      function(from) new("lgCMatrix", i = from@i, p = from@p,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("matrix", "dgCMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .Call("matrix_to_csc", from, PACKAGE = "Matrix")
      })

setAs("dgeMatrix", "dgCMatrix",
      function(from) .Call("dgeMatrix_to_csc", from, PACKAGE = "Matrix"))


setMethod("crossprod", signature(x = "dgCMatrix", y = "missing"),
          function(x, y = NULL) .Call("csc_crossprod", x, PACKAGE = "Matrix"),
          valueClass = "dsCMatrix")

setMethod("crossprod", signature(x = "dgCMatrix", y = "dgeMatrix"),
          function(x, y = NULL)
          .Call("csc_matrix_crossprod", x, y, TRUE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dgCMatrix", y = "matrix"),
          function(x, y = NULL)
          .Call("csc_matrix_crossprod", x, y, FALSE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

##setMethod("crossprod", signature(x = "dgCMatrix", y = "numeric"),
##          function(x, y = NULL) callGeneric(x, as.matrix(y)),
##          valueClass = "dgeMatrix")

## setMethod("crossprod", signature(x = "dgCMatrix", y = "numeric"),
##           function(x, y = NULL) .Call("csc_matrix_crossprod", x, as.matrix(y)))

setMethod("tcrossprod", signature(x = "dgCMatrix", y = "missing"),
          function(x, y = NULL) .Call("csc_tcrossprod", x, PACKAGE = "Matrix"))

setMethod("diag", signature(x = "dgCMatrix"),
          function(x = 1, nrow, ncol = n)
          .Call("csc_getDiag", x, PACKAGE = "Matrix"))

## try to define for "Matrix" -- once and for all -- but that fails -- why?
setMethod("dim", signature(x = "dgCMatrix"),
          function(x) x@Dim, valueClass = "integer")

setMethod("t", signature(x = "dgCMatrix"),
          function(x) .Call("csc_transpose", x, PACKAGE = "Matrix"),
          valueClass = "dgCMatrix")

setMethod("image", "dgCMatrix",
          function(x, ...) {
              x <- as(x, "dgTMatrix")
              callGeneric()
          })

setMethod("%*%", signature(x = "dgCMatrix", y = "dgeMatrix"),
          function(x, y) .Call("csc_matrix_mm", x, y, TRUE, FALSE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgCMatrix", y = "matrix"),
          function(x, y) .Call("csc_matrix_mm", x, y, FALSE, FALSE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgeMatrix", y = "dgCMatrix"),
          function(x, y) .Call("csc_matrix_mm", y, x, TRUE, TRUE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "matrix", y = "dgCMatrix"),
          function(x, y) .Call("csc_matrix_mm", y, x, FALSE, TRUE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

## Group Methods, see ?Arith (e.g.)
## -----

### TODO:

if(FALSE) ## FIXME
setMethod("Arith", ##  "+", "-", "*", "^", "%%", "%/%", "/"
          signature(e1 = "dgCMatrix", e2 = "dgCMatrix"),
          function(e1, e2) {
              d <- dimCheck(e1, e2)
              dn <- dimNamesCheck(e1, e2)
              ij1 <- non0ind(e1)
              ij2 <- non0ind(e2)
              switch(.Generic,
                     "+" =, "-" =, "*" =
                     new("dgTMatrix", Dim = d, Dimnames = dn,
                         i = c(ij1[,1], ij2[,1]),
                         j = c(ij1[,2], ij2[,2]),
                         x = c(callGeneric(e1@x, 0), callGeneric(0, e2@x)))
                     ,
                     "^" = { ## X^0 |-> 1 (also for X=0)
                         r <- new("dgTMatrix", Dim = d, Dimnames = dn,
                                  i = c(ij1[,1], ij2[,1]),
                                  j = c(ij1[,2], ij2[,2]),
                                  x = c(rep.int(1, nrow(ij1)), 0 ^ e2@x))
                         ...
                     },
                     "%%" = , "%/%" = , "/" = {## 0 op 0  |-> NaN
                         ...
                     })
          })

setMethod("Arith",
	  signature(e1 = "dgCMatrix", e2 = "numeric"),
	  function(e1, e2) {
	      if(length(e2) == 1) {
		  f0 <- callGeneric(0, e2)
                  if(!is.na(f0) && f0 == 0.) {
		      e1@x <- callGeneric(e1@x, e2)
		      e1
		  } else {
		      ## FIXME: dgeMatrix [cbind(i,j)] <- .. is not yet possible
		      ##		  r <- as(e1, "dgeMatrix")
		      ##		  r[] <- f0
		      ##		  r[non0ind(e1)] <- callGeneric(e1@x, e2)
		      r <- as(e1, "matrix")
		      r[] <- f0
		      r[non0ind(e1)] <- callGeneric(e1@x, e2)
		      as(r, "dgeMatrix")
		  }
	      } else {
		  ## FIXME: maybe far from optimal:
		  warning("coercing sparse to dense matrix for arithmetic")
		  callGeneric(as(e1, "dgeMatrix"), e2)
	      }
	  })

setMethod("Arith",
	  signature(e1 = "numeric", e2 = "dgCMatrix"),
	  function(e1, e2) {
	      if(length(e1) == 1) {
		  f0 <- callGeneric(e1, 0)
                  if(!is.na(f0) && f0 == 0.) {
		      e2@x <- callGeneric(e1, e2@x)
		      e2
		  } else {
		      ## FIXME: dgeMatrix [cbind(i,j)] <- .. is not yet possible
		      r <- as(e2, "matrix")
		      r[] <- f0
		      r[non0ind(e2)] <- callGeneric(e1, e2@x)
		      as(r, "dgeMatrix")
		  }
	      } else {
		  ## FIXME: maybe far from optimal:
		  warning("coercing sparse to dense matrix for arithmetic")
		  callGeneric(e1, as(e2, "dgeMatrix"))
	      }
	  })


setMethod("Math",
	  signature(x = "dgCMatrix"),
	  function(x) {
              f0 <- callGeneric(0.)
	      if(!is.na(f0) && f0 == 0.) {
		  ## sparseness preserved
		  x@x <- callGeneric(x@x)
		  x
	      } else { ## no sparseness
		  callGeneric(as(x, "dgeMatrix"))
	      }
	  })

if(FALSE) ## unneeded with "Math2" in ./dMatrix.R
setMethod("Math2",
	  signature(x = "dgCMatrix", digits = "numeric"),
	  function(x, digits) {
	      f0 <- callGeneric(0., digits = digits)
	      if(!is.na(f0) && f0 == 0.) {
		  ## sparseness preserved
		  x@x <- callGeneric(x@x, digits = digits)
		  x
	      } else { ## no sparseness
		  callGeneric(as(x, "dgeMatrix"), digits = digits)
	      }
	  })

###---- end {Group Methods} -----------------



setMethod("writeHB", signature(obj = "dgCMatrix"),
          function(obj, file, ...)
          .Call("Matrix_writeHarwellBoeing", obj, as.character(file), "DGC", PACKAGE = "Matrix"))

setMethod("writeMM", signature(obj = "dgCMatrix"),
          function(obj, file, ...)
          .Call("Matrix_writeMatrixMarket", obj, as.character(file), "DGC", PACKAGE = "Matrix"))
