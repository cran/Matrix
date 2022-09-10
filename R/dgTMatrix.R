## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
## Now in ./Tsparse.R
## setAs("dgTMatrix", "dgCMatrix",
##       function(from) .Call(Tsparse_to_Csparse, from, FALSE)
##       )

setAs("dgTMatrix", "dgeMatrix",
      function(from) .Call(dgTMatrix_to_dgeMatrix, from))

setAs("dgTMatrix", "matrix",
      function(from) .Call(dgTMatrix_to_matrix, from))


if(FALSE) ## special case, relatively ugly, needed ??
setAs("dgTMatrix", "dsCMatrix",
      function(from) {
          if (!isSymmetric(from))
	      stop("cannot coerce non-symmetric \"dgTMatrix\" to \"dsCMatrix\" class")
          upper <- from@i <= from@j
          uC <- as(new("dgTMatrix", Dim = from@Dim,  Dimnames = from@Dimnames,
                       i = from@i[upper],
                       j = from@j[upper], x = from@x[upper]), "dgCMatrix")
          new("dsCMatrix", Dim = uC@Dim, p = uC@p, i = uC@i, x = uC@x, uplo = "U")
      })

## This is faster:
setAs("dgTMatrix", "dtCMatrix",
      function(from) {
	  if(!(iTri <- isTriangular(from)))
	      stop("the matrix is not triangular")
	  ## else
	  stopifnot(is.character(uplo <- attr(iTri, "kind")))
	  .Call(Tsparse_to_tCsparse, from, uplo, "N")
      })

setAs("dgTMatrix", "dtTMatrix",
      function(from) check.gT2tT(from, toClass = "dtTMatrix", do.n=FALSE))
setAs("dgTMatrix", "dsTMatrix",
      function(from) check.gT2sT(from, toClass = "dsTMatrix", do.n=FALSE))
} ## MJ

## MJ: no longer needed ... methods now inherited from Matrix
if(FALSE) {
setAs("dgTMatrix", "triangularMatrix",
      function(from) check.gT2tT(from, toClass = "dtTMatrix", do.n=FALSE))
setAs("dgTMatrix", "symmetricMatrix",
      function(from) check.gT2sT(from, toClass = "dsTMatrix", do.n=FALSE))
} ## MJ

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
mat2dgT <- function(from) {
    x <- as.double(from)
    nz <- isN0(x)
    new("dgTMatrix", Dim = dim(from),
        Dimnames = .M.DN(from),
        i = row(from)[nz] - 1L,
        j = col(from)[nz] - 1L,
        x = x[nz])
}
setAs("matrix", "dgTMatrix", mat2dgT)

setAs("dgeMatrix", "dgTMatrix",
      function(from) as(as(from, "dgCMatrix"), "dgTMatrix"))
} ## MJ

## "[" methods are now in ./Tsparse.R

## "[<-" methods { setReplaceMethod()s }  too ...

## Uses the triplet convention of *adding* entries with same (i,j):
setMethod("+", signature(e1 = "dgTMatrix", e2 = "dgTMatrix"),
          function(e1, e2) {
              dimCheck(e1, e2)
              new("dgTMatrix", i = c(e1@i, e2@i), j = c(e1@j, e2@j),
                  x = c(e1@x, e2@x), Dim = e1@Dim, Dimnames = e1@Dimnames)
          })


## setMethod("writeHB", signature(obj = "dgTMatrix"),
## 	  function(obj, file, ...) callGeneric(as(obj, "CsparseMatrix"), file, ...))
