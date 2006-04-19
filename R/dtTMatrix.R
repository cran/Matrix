### Coercion and Methods for Triangular Triplet Matrices

gt2tT <- function(x, uplo, diag) {
    ## coerce *gTMatrix to *tTMatrix {general -> triangular}
    i <- x@i
    j <- x@j
    sel <-
	if(uplo == "U") {
	    if(diag == "U") i < j else i <= j
	} else {
	    if(diag == "U") i > j else i >= j
	}
    i <- i[sel]
    j <- j[sel]
    if(is(x, "lMatrix"))
	new("ltTMatrix", i = i, j = j, uplo = uplo, diag = diag,
	    Dim = x@Dim, Dimnames = x@Dimnames) # no 'x' slot
    else
	new(paste(substr(class(x), 1,1), # "d", "l", "i" or "z"
		  "tTMatrix", sep=''),
	    i = i, j = j, uplo = uplo, diag = diag,
	    x = x@x[sel], Dim = x@Dim, Dimnames = x@Dimnames)
}

setAs("dtTMatrix", "dtCMatrix",
      function(from) {
          gC <- .Call("dtTMatrix_as_dgCMatrix", from, PACKAGE = "Matrix")
          new("dtCMatrix", Dim = gC@Dim, Dimnames = gC@Dimnames, p = gC@p,
              i = gC@i, x = gC@x, uplo = from@uplo, diag = from@diag)
      })

setAs("dtTMatrix", "dgTMatrix",
      function(from) {
          d <- from@Dim
          if(uDiag <- from@diag == "U") # unit diagonal, need to add '1's
              uDiag <- (n <- d[1]) > 0
          new("dgTMatrix", Dim = d, Dimnames = from@Dimnames,
              i = c(from@i, if(uDiag) 0:(n-1)),
              j = c(from@j, if(uDiag) 0:(n-1)),
              x = c(from@x, if(uDiag) rep.int(1,n)))
      })

setAs("dtTMatrix", "ltTMatrix",
      function(from) new("ltTMatrix", i = from@i, j = from@j,
                         uplo = from@uplo, diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

## Conversion to dense storage is first to a dtrMatrix
setAs("dtTMatrix", "dtrMatrix",
      function(from) .Call("dtTMatrix_as_dtrMatrix", from, PACKAGE = "Matrix"))

setAs("dtTMatrix", "matrix",
      function(from) as(as(from, "dtrMatrix"), "matrix"))

setAs("dtTMatrix", "dgeMatrix",
      function(from) as(as(from, "dtrMatrix"), "dgeMatrix"))

setAs("matrix", "dtTMatrix",
      function(from) as(as(from, "dtpMatrix"), "dtTMatrix"))


setMethod("t", signature(x = "dtTMatrix"),
          function(x)
          new("dtTMatrix", Dim = rev(x@Dim), diag = x@diag,
              i = x@j, j = x@i, x = x@x,
              uplo = if (x@uplo == "U") "L" else "U"),
          valueClass = "dtTMatrix")
