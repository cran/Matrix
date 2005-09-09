#### Logical Sparse Matrices in triplet format

### contains = "lsparseMatrix"
###             ============= ---> superclass methods in ./lsparseMatrix.R

setAs("lgTMatrix", "lgCMatrix",
      function(from) .Call("lgTMatrix_as_lgCMatrix", from))

setAs("lgTMatrix", "matrix",
      function(from) as(as(from, "lgCMatrix"), "matrix"))

setMethod("t", signature(x = "lgTMatrix"),
	  function(x) new("lgTMatrix", i = x@j, j = x@i,
			  Dim = x@Dim[2:1],
			  Dimnames= x@Dimnames[2:1]),
	  valueClass = "lgTMatrix")
