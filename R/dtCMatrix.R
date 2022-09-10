#### Triangular Sparse Matrices in compressed column-oriented format

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
setAs("matrix", "dtCMatrix",
      function(from) as(as(from, "dtTMatrix"), "dtCMatrix"))
} ## MJ

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
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

##' dtC* |-> dgC*  (provide for direct use in other coercions) :
.dtC2g <- function(from) {
    if (from@diag == "U")
        from <- .Call(Csparse_diagU2N, from)
    ## new("dgCMatrix", .....) # ---> Rather faster, no checking:
    copyClass(from, "dgCMatrix",
              sNames = c("i", "p", "x", "Dim", "Dimnames"), check = FALSE)
}
setAs("dtCMatrix", "dgCMatrix", .dtC2g)

setAs("dtCMatrix", "dsCMatrix", function(from) as(from, "symmetricMatrix"))

## FIXME: make more efficient
## -----  and  as(., "triangularMatrix") is even worse via as_Sp()
setAs("dgCMatrix", "dtCMatrix", # to triangular, needed for triu,..
      function(from) as(.Call(Csparse_to_Tsparse, from, FALSE), "dtCMatrix"))

setAs("dtCMatrix", "dgTMatrix",
      function(from) {
          if (from@diag == "U") from <- .Call(Csparse_diagU2N, from)
          ## ignore triangularity in conversion to TsparseMatrix
          .Call(Csparse_to_Tsparse, from, FALSE)
      })

setAs("dtCMatrix", "dgeMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgeMatrix"))

## These are all needed because cholmod doesn't support triangular:
## (see end of ./Csparse.R ), e.g. for triu()
setAs("dtCMatrix", "dtTMatrix",
      function(from) .Call(Csparse_to_Tsparse, from, TRUE))
##   {# and this is not elegant:
##           x <- as(from, "dgTMatrix")
##  	  if (from@diag == "U") { ## drop diagonal entries '1':
##  	      i <- x@i; j <- x@j
##  	      nonD <- i != j
##  	      xx <- x@x[nonD] ; i <- i[nonD] ; j <- j[nonD]
##  	  } else {
##  	      xx <- x@x; i <- x@i; j <- x@j
##  	  }
##  	  new("dtTMatrix", x = xx, i = i, j = j, Dim = x@Dim,
##  	      Dimnames = x@Dimnames, uplo = from@uplo, diag = from@diag)
##       })

## Now that we support triangular matrices use the inherited method.
## setAs("dtCMatrix", "TsparseMatrix", function(from) as(from, "dtTMatrix"))

setAs("dtCMatrix", "dtrMatrix",
      function(from) as(as(from, "dtTMatrix"), "dtrMatrix"))
} ## MJ
