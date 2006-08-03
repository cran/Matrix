setMethod("t", signature(x = "dtCMatrix"),
          function(x) {
              tg <- .Call(csc_transpose, x)
              new("dtCMatrix", Dim = tg@Dim, Dimnames = x@Dimnames[2:1],
                  p = tg@p, i = tg@i, x = tg@x, diag = x@diag,
                  uplo = ifelse(x@uplo == "U", "L", "U"))
          }, valueClass = "dtCMatrix")

setAs("dtCMatrix", "ltCMatrix", # just drop 'x' slot:
      function(from) new("ltCMatrix", i = from@i, p = from@p,
                         uplo = from@uplo, diag = from@diag,
                         ## FIXME?: use from@factors smartly
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("matrix", "dtCMatrix",
      function(from) as(as(from, "dtTMatrix"), "dtCMatrix"))

setAs("dtCMatrix", "dgCMatrix",
      function(from) {
          if(from@diag == "U") { ## add diagonal of 1's
              ##FIXME: do this smartly - directly {in C or R}
              as(as(from, "dgTMatrix"), "dgCMatrix")
          }
          else
              new("dgCMatrix",
                  i = from@i, p = from@p, x = from@x,
                  Dim = from@Dim, Dimnames = from@Dimnames)
      })

setAs("dgCMatrix", "dtCMatrix", # to triangular:
      function(from) as(as(as(from, "dgTMatrix"), "dtTMatrix"), "dtCMatrix"))

setAs("dtCMatrix", "dgTMatrix",
      function(from)
      .Call(tsc_to_dgTMatrix, from))

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

setAs("dtCMatrix", "TsparseMatrix", function(from) as(from, "dtTMatrix"))

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
	      if (a@diag == "U") a <- as(diagU2N(a), "dtCMatrix")
	      .Call(dtCMatrix_matrix_solve, a, b, TRUE)
	  }, valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dtCMatrix", b = "matrix"),
	  function(a, b, ...) {
	      if (a@diag == "U") a <- as(diagU2N(a), "dtCMatrix")
	      storage.mode(b) <- "double"
	      .Call(dtCMatrix_matrix_solve, a, b, FALSE)
	  }, valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dtCMatrix", b = "numeric"),
	  function(a, b, ...) {
	      if (a@diag == "U") a <- as(diagU2N(a), "dtCMatrix")
	      .Call(dtCMatrix_matrix_solve, a, as.matrix(as.double(b)),
		    FALSE)
	  }, valueClass = "dgeMatrix")
