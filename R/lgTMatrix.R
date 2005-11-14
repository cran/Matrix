#### Logical Sparse Matrices in triplet format

### contains = "lsparseMatrix"
###             ============= ---> superclass methods in ./lsparseMatrix.R

setAs("lgTMatrix", "lgCMatrix",
      function(from)
      .Call("lgTMatrix_as_lgCMatrix", from, PACKAGE = "Matrix"))

setAs("lgTMatrix", "matrix",
      function(from) as(as(from, "lgCMatrix"), "matrix"))

setAs("lgTMatrix", "dgTMatrix",
      function(from)
      ## more efficient than
      ## as(as(as(sM, "lgCMatrix"), "dgCMatrix"), "dgTMatrix")
      new("dgTMatrix", i = from@i, j = from@j,
          x = rep.int(1, length(from@i)),
          ## cannot copy factors, but can we use them?
          Dim = from@Dim, Dimnames= from@Dimnames))

setMethod("t", signature(x = "lgTMatrix"),
	  function(x) new("lgTMatrix", i = x@j, j = x@i,
			  Dim = x@Dim[2:1],
			  Dimnames= x@Dimnames[2:1]),
	  valueClass = "lgTMatrix")
