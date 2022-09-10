#### Logical Sparse Matrices in Compressed column-oriented format

### contains = "lsparseMatrix"

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
setAs("lgCMatrix", "dgCMatrix",
      function(from) new("dgCMatrix", i = from@i, p = from@p,
                         x = as.double(from@x),
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("lgCMatrix", "lgTMatrix",
      function(from) new("lgTMatrix", i = from@i, x = from@x,
                         j = .Call(Matrix_expand_pointers, from@p),
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("lgCMatrix", "lgeMatrix",
      function(from)
	  new("lgeMatrix", x = c(as(from, "matrix")), # is fast,
	      Dim = from@Dim, Dimnames = from@Dimnames))

setAs("lgCMatrix", "matrix", function(from) .Call(lgC_to_matrix, from))
## not this: .Call(Csparse_to_matrix, from)), since it goes via dense -> double precision
} ## MJ

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
setAs("matrix", "lgCMatrix", .m2lgC)
} ## MJ
