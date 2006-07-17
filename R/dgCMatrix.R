#### Sparse Matrices in Compressed column-oriented format

### contains = "dsparseMatrix"

setAs("dgCMatrix", "dgTMatrix",
      function(from) .Call(compressed_to_dgTMatrix, from, TRUE))

setAs("dgCMatrix", "matrix",
      function(from) .Call(csc_to_matrix, from))

setAs("dgCMatrix", "dgeMatrix",
      function(from) .Call(csc_to_dgeMatrix, from))

setAs("dgCMatrix", "lgCMatrix",
      function(from) new("lgCMatrix", i = from@i, p = from@p,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("matrix", "dgCMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .Call(matrix_to_csc, from)
      })

setAs("dgeMatrix", "dgCMatrix",
      function(from) .Call(dgeMatrix_to_csc, from))


setMethod("crossprod", signature(x = "dgCMatrix", y = "missing"),
          function(x, y = NULL) .Call(csc_crossprod, x),
          valueClass = "dsCMatrix")

setMethod("crossprod", signature(x = "dgCMatrix", y = "dgeMatrix"),
          function(x, y = NULL)
          .Call(csc_matrix_crossprod, x, y, TRUE),
          valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dgCMatrix", y = "matrix"),
          function(x, y = NULL) {
	      storage.mode(y) <- "double"
              .Call(csc_matrix_crossprod, x, y, FALSE)
          }, valueClass = "dgeMatrix")

##setMethod("crossprod", signature(x = "dgCMatrix", y = "numeric"),
##          function(x, y = NULL) callGeneric(x, as.matrix(y)),
##          valueClass = "dgeMatrix")

## setMethod("crossprod", signature(x = "dgCMatrix", y = "numeric"),
##           function(x, y = NULL) .Call(csc_matrix_crossprod, x, as.matrix(y)))

setMethod("tcrossprod", signature(x = "dgCMatrix", y = "missing"),
          function(x, y = NULL) .Call(csc_tcrossprod, x))

setMethod("diag", signature(x = "dgCMatrix"),
	  function(x, nrow, ncol = n) .Call(csc_getDiag, x))

## try to define for "Matrix" -- once and for all -- but that fails -- why? __ FIXME __
## setMethod("dim", signature(x = "dgCMatrix"),
##           function(x) x@Dim, valueClass = "integer")

setMethod("t", signature(x = "dgCMatrix"),
          function(x) .Call(csc_transpose, x),
          valueClass = "dgCMatrix")

setMethod("image", "dgCMatrix",
          function(x, ...) {
              x <- as(x, "dgTMatrix")
              callGeneric()
          })

setMethod("%*%", signature(x = "dgCMatrix", y = "dgeMatrix"),
          function(x, y) .Call(csc_matrix_mm, x, y, TRUE, FALSE),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgCMatrix", y = "matrix"),
          function(x, y) {
	      storage.mode(y) <- "double"
              .Call(csc_matrix_mm, x, y, FALSE, FALSE)
          }, valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgeMatrix", y = "dgCMatrix"),
          function(x, y) .Call(csc_matrix_mm, y, x, TRUE, TRUE),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "matrix", y = "dgCMatrix"),
          function(x, y) {
	      storage.mode(x) <- "double"
              .Call(csc_matrix_mm, y, x, FALSE, TRUE)
          }, valueClass = "dgeMatrix")

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
                     Yis0 <- r == 0
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
                  if(!is.na(f0) && f0 == 0.) { # remain sparse
		      e1@x <- callGeneric(e1@x, e2)
		      e1
		  } else { ## non-sparse, since '0 o e2' is not 0

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


replCmat <- function (x, i, j, value)
{
    di <- dim(x)
    dn <- dimnames(x)
    i1 <- if(missing(i)) 0:(di[1] - 1:1) else .ind.prep2(i, 1, di, dn)
    i2 <- if(missing(j)) 0:(di[2] - 1:1) else .ind.prep2(j, 2, di, dn)
    dind <- c(length(i1), length(i2)) # dimension of replacement region
    lenRepl <- prod(dind)
    lenV <- length(value)
    if(lenV == 0) {
        if(lenRepl != 0)
            stop("nothing to replace with")
        else return(x)
    }
    ## else: lenV := length(value)	 is > 0
    if(lenRepl %% lenV != 0)
        stop("number of items to replace is not a multiple of replacement length")
    if(lenV > lenRepl)
        stop("too many replacement values")

    xj <- .Call(Matrix_expand_pointers, x@p)
    sel <- (!is.na(match(x@i, i1)) &
            !is.na(match( xj, i2)))

    if(sum(sel) == lenRepl) { ## all entries to be replaced are non-zero:
        value <- rep(value, length = lenRepl)
        ## Ideally we only replace them where value != 0 and drop the value==0
        ## ones; but that would have to (?) go through dgT*
        ## v0 <- 0 == value
        ## if (lenRepl == 1) and v0 is TRUE, the following is not doing anything
        ##-  --> ./dgTMatrix.R  and its  replTmat()
        ## x@x[sel[!v0]] <- value[!v0]
        x@x[sel] <- value
        return(x)
    }
    ## else go via dgT
    x <- as(x, "dgTMatrix")
    x[i,j] <- value
    as(x, "dgCMatrix")
}

### TODO (FIXME): almost the same for  "lgCMatrix" and "logical"

setReplaceMethod("[", signature(x = "dgCMatrix", i = "index", j = "missing",
                                value = "numeric"),
                 function (x, i, value) replCmat(x, i=i, value=value))

setReplaceMethod("[", signature(x = "dgCMatrix", i = "missing", j = "index",
                                value = "numeric"),
                 function (x, j, value) replCmat(x, j=j, value=value))

setReplaceMethod("[", signature(x = "dgCMatrix", i = "index", j = "index",
				value = "numeric"),
                 replCmat)




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

