#### Logical Sparse Matrices in triplet format

### contains = "lsparseMatrix"

setAs("lgTMatrix", "lgCMatrix",
      function(from) .Call("lgTMatrix_as_lgCMatrix", from))

setMethod("t", signature(x = "lgTMatrix"),
          function(x) new("lgTMatrix", i = from@j, j = from@i,
                          Dim = from@Dim[c(2,1)],
                          DimNames=from@DimNames[c(2,1)]),
          valueClass = "lgTMatrix")

setMethod("crossprod", signature(x = "lgTMatrix", y = "missing"),
	  function(x, y = NULL)
          .Call("lgCMatrix_crossprod", as(x, "lgCMatrix"), TRUE, NULL),
	  valueClass = "lsCMatrix")

setMethod("tcrossprod", signature(x = "lgTMatrix"),
	  function(x)
          .Call("lgCMatrix_crossprod", as(x, "lgCMatrix"), FALSE, NULL),
	  valueClass = "lsCMatrix")

