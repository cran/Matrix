#### Triangular Sparse Matrices in compressed column-oriented format

setAs("dtCMatrix", "ltCMatrix",
      function(from) new("ltCMatrix", i = from@i, p = from@p,
			 uplo = from@uplo, diag = from@diag,
                         x = as.logical(from@x),
			 ## FIXME?: use from@factors smartly
			 Dim = from@Dim, Dimnames = from@Dimnames))
setAs("dtCMatrix", "ntCMatrix", # just drop 'x' slot:
      function(from) new("ntCMatrix", i = from@i, p = from@p,
			 uplo = from@uplo, diag = from@diag,
			 ## FIXME?: use from@factors smartly
			 Dim = from@Dim, Dimnames = from@Dimnames))


setAs("matrix", "dtCMatrix",
      function(from) as(as(from, "dtTMatrix"), "dtCMatrix"))

setAs("dtCMatrix", "dgCMatrix",
      function(from) {
          if (from@diag == "U")
              from <- .Call(Csparse_diagU2N, from)
          new("dgCMatrix",
              i = from@i, p = from@p, x = from@x,
              Dim = from@Dim, Dimnames = from@Dimnames)
      })

setAs("dtCMatrix", "dgTMatrix",
      function(from) {
          if (from@diag == "U") from <- .Call(Csparse_diagU2N, from)
          ## ignore triangularity in conversion to TsparseMatrix
          .Call(Csparse_to_Tsparse, from, FALSE)
      })

setAs("dgCMatrix", "dtCMatrix", # to triangular:
      function(from) as(as(as(from, "dgTMatrix"), "dtTMatrix"), "dtCMatrix"))

setAs("dtCMatrix", "dgeMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgeMatrix"))

## These are all needed because cholmod doesn't support triangular:
## (see end of ./Csparse.R )
setAs("dtCMatrix", "dtTMatrix",
      function(from) {# and this is not elegant:
          x <- as(from, "dgTMatrix")
 	  if (from@diag == "U") { ## drop diagonal entries '1':
 	      i <- x@i; j <- x@j
 	      nonD <- i != j
 	      xx <- x@x[nonD] ; i <- i[nonD] ; j <- j[nonD]
 	  } else {
 	      xx <- x@x; i <- x@i; j <- x@j
 	  }
 	  new("dtTMatrix", x = xx, i = i, j = j, Dim = x@Dim,
 	      Dimnames = x@Dimnames, uplo = from@uplo, diag = from@diag)
      })

## Now that we support triangular matrices use the inherited method.
## setAs("dtCMatrix", "TsparseMatrix", function(from) as(from, "dtTMatrix"))

setAs("dtCMatrix", "dtrMatrix",
      function(from) as(as(from, "dtTMatrix"), "dtrMatrix"))

## using  diagU2N() from ./Auxiliaries.R :
setMethod("solve", signature(a = "dtCMatrix", b = "missing"),
	  function(a, b, ...) {
	      if (a@diag == "U") {
		  if (a@uplo == "U")
		      return(.Call(dtCMatrix_upper_solve, a))
		  else
		      return(t(.Call(dtCMatrix_upper_solve, t(a))))
	      }
	      .Call(dtCMatrix_solve, a)
	  }, valueClass = "dtCMatrix")

setMethod("solve", signature(a = "dtCMatrix", b = "dgeMatrix"),
	  function(a, b, ...) {
#	      if (a@diag == "U") a <- as(diagU2N(a), "dtCMatrix")
              if (a@diag == "U") a <- .Call(Csparse_diagU2N, a)
	      .Call(dtCMatrix_matrix_solve, a, b, TRUE)
	  }, valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dtCMatrix", b = "matrix"),
	  function(a, b, ...) {
#	      if (a@diag == "U") a <- as(diagU2N(a), "dtCMatrix")
              if (a@diag == "U") a <- .Call(Csparse_diagU2N, a)
	      storage.mode(b) <- "double"
	      .Call(dtCMatrix_matrix_solve, a, b, FALSE)
	  }, valueClass = "dgeMatrix")

## Isn't this case handled by the method for (a = "Matrix', b =
## "numeric") in ./Matrix.R? Or is this method defined here for
## the as.double coercion?
setMethod("solve", signature(a = "dtCMatrix", b = "numeric"),
	  function(a, b, ...) {
	      if (a@diag == "U") a <- as(diagU2N(a), "dtCMatrix")
	      .Call(dtCMatrix_matrix_solve, a, as.matrix(as.double(b)),
		    FALSE)
	  }, valueClass = "dgeMatrix")
