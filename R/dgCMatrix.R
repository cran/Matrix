#### Sparse Matrices in Compressed column-oriented format

### contains = "dsparseMatrix", "CsparseMatrix"

## Specific conversions, should they be necessary.  Better to convert as
## as(x, "TsparseMatrix") or as(x, "denseMatrix")

## Moved to ./Csparse.R :
## setAs("dgCMatrix", "dgTMatrix", ....
## setAs("dgCMatrix", "dgeMatrix", ....
## setAs("dgeMatrix", "dgCMatrix", ....

## rather use Csparse* to lsparse* in ./lsparseMatrix.R ,
## but this is for "back-compatibility" (have had tests for it..):

setAs("dgCMatrix", "ngCMatrix",
      function(from) .Call(Csparse_to_nz_pattern, from, FALSE))

setAs("dgCMatrix", "lgCMatrix",
      function(from) { ## FIXME use .Call() too!
	  r <- new("lgCMatrix")
	  r@x <- as.logical(from@x)
	  ## and copy the other slots
	  for(nm in c("i", "p", "Dim", "Dimnames"))
	      slot(r, nm) <- slot(from, nm)
          r
      })

setMethod("image", "dgCMatrix",
          function(x, ...) {
              x <- as(x, "dgTMatrix")
              callGeneric()
          })

## Group Methods, see ?Arith (e.g.)
## -----

setMethod("Arith", ##  "+", "-", "*", "^", "%%", "%/%", "/"
          signature(e1 = "dgCMatrix", e2 = "dgCMatrix"),
          function(e1, e2) {
              d <- dimCheck(e1, e2)
              dn <- dimNamesCheck(e1, e2)
              ij1 <- non0ind(e1)
              ij2 <- non0ind(e2)
              switch(.Generic,
                     "+" = , "-" =
                     ## special "T" convention: repeated entries are *summed*
                     as(new("dgTMatrix", Dim = d, Dimnames = dn,
                            i = c(ij1[,1], ij2[,1]),
                            j = c(ij1[,2], ij2[,2]),
                            x = c(callGeneric(e1@x, 0), callGeneric(0,e2@x))),
                        "dgCMatrix"),

                     "*" =
                 { ##  X * 0 == 0 * X == 0 --> keep common non-0
                     ii <- WhichintersectInd(ij1, ij2, nrow=d[1])
                     ij <- ij1[ii[[1]], , drop = FALSE]
                     as(new("dgTMatrix", Dim = d, Dimnames = dn,
                            i = ij[,1],
                            j = ij[,2],
                            x = e1@x[ii[[1]]] * e2@x[ii[[2]]]),
                        "dgCMatrix")
                 },

                     "^" =
                 {
                     ii <- WhichintersectInd(ij1, ij2, nrow=d[1])
                     ## 3 cases:
                     ## 1) X^0 := 1  (even for X=0) ==> dense
                     ## 2) 0^Y := 0  for Y != 0         =====
                     ## 3) x^y :

                     ## FIXME:  dgeM[cbind(i,j)] <- V  is not yet possible
                     ##     nor dgeM[ i_vect   ] <- V
                     ## r <- as(e2, "dgeMatrix")
                     ## ...
                     r <- as(e2, "matrix")
                     Yis0 <- is0(r)
                     r[complementInd(ij1, dim=d)] <- 0      ## 2)
                     r[1:1 + ij2[ii[[2]], , drop=FALSE]] <-
                         e1@x[ii[[1]]] ^ e2@x[ii[[2]]]      ## 3)
                     r[Yis0] <- 1                           ## 1)
                     as(r, "dgeMatrix")
                 },

                     "%%" = , "%/%" = , "/" = ## 0 op 0  |-> NaN => dense
                     callGeneric(as(e1, "dgeMatrix"), e2)
                     )
          })

setMethod("Arith",
	  signature(e1 = "dgCMatrix", e2 = "numeric"),
	  function(e1, e2) {
	      if(length(e2) == 1) { ## e.g.,  Mat ^ a
		  f0 <- callGeneric(0, e2)
		  if(is0(f0)) { # remain sparse
		      e1@x <- callGeneric(e1@x, e2)
		      e1
		  } else { ## non-sparse, since '0 o e2' is not 0

		      ## FIXME: dgeMatrix [cbind(i,j)] <- .. is not yet possible
		      ##		  r <- as(e1, "dgeMatrix")
		      ##		  r[] <- f0
		      ##		  r[non0ind(e1)] <- callGeneric(e1@x, e2)
		      r <- as(e1, "matrix")
		      r[] <- f0
		      r[non0ind(e1) + 1:1] <- callGeneric(e1@x, e2)
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
                  if(is0(f0)) {
		      e2@x <- callGeneric(e1, e2@x)
		      e2
		  } else {
		      ## FIXME: dgeMatrix [cbind(i,j)] <- .. is not yet possible
		      r <- as(e2, "matrix")
		      r[] <- f0
		      r[non0ind(e2) + 1:1] <- callGeneric(e1, e2@x)
		      as(r, "dgeMatrix")
		  }
	      } else {
		  ## FIXME: maybe far from optimal:
		  warning("coercing sparse to dense matrix for arithmetic")
		  callGeneric(e1, as(e2, "dgeMatrix"))
	      }
	  })


## "Math" is up in ./Csparse.R

## "Math2" is up in ./dMatrix.R


###---- end {Group Methods} -----------------


## "[<-" methods { setReplaceMethod()s }  are now in ./Csparse.R


setMethod("writeHB", signature(obj = "dgCMatrix"),
          function(obj, file, ...)
          .Call(Matrix_writeHarwellBoeing, obj, as.character(file), "DGC"))

setMethod("writeMM", signature(obj = "dgCMatrix"),
          function(obj, file, ...)
          .Call(Matrix_writeMatrixMarket, obj, as.character(file), "DGC"))


## TODO (in C):
## setMethod("colSums", signature(x = "dgCMatrix"),
## 	  function(x, na.rm = FALSE, dims = 1)
##           .Call(dgCMatrix_colsums, x, na.rm, TRUE, FALSE),
## 	  valueClass = "numeric")

## setMethod("colMeans", signature(x = "dgCMatrix"),
## 	  function(x, na.rm = FALSE, dims = 1)
##           .Call(dgCMatrix_colsums, x, na.rm, TRUE, TRUE),
## 	  valueClass = "numeric")

setMethod("colSums",  signature(x = "dgCMatrix"), .as.dgT.Fun)
setMethod("colMeans", signature(x = "dgCMatrix"), .as.dgT.Fun)

setMethod("rowSums", signature(x = "dgCMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          tapply1(x@x, factor(x@i, 0:(x@Dim[1]-1)), sum, na.rm = na.rm),
	  valueClass = "numeric")

setMethod("rowMeans", signature(x = "dgCMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          tapply1(x@x, factor(x@i, 0:(x@Dim[1]-1)), mean, na.rm = na.rm),
	  valueClass = "numeric")

setMethod("qr", signature(x = "dgCMatrix"),
          function(x, tol = 1e-07, LAPACK = FALSE)
          .Call(dgCMatrix_QR, x, TRUE))

setMethod("lu", signature(x = "dgCMatrix"),
          function(x, ...) .Call(dgCMatrix_LU, x, TRUE, 1))

setMethod("solve", signature(a = "dgCMatrix", b = "matrix"),
          function(a, b, ...) .Call(dgCMatrix_matrix_solve, a, b),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgCMatrix", b = "ddenseMatrix"),
          function(a, b, ...) .Call(dgCMatrix_matrix_solve, a, b),
          valueClass = "dgeMatrix")
