#### Nonzero Pattern Sparse Matrices in triplet format

### contains = "nsparseMatrix"
###             ============= ---> superclass methods in ./nsparseMatrix.R


setAs("ngTMatrix", "matrix", # go via fast C code:
      function(from) as(as(from, "ngCMatrix"), "matrix"))

setAs("matrix", "ngTMatrix",
      function(from) {
	  stopifnot(is.logical(from))
	  if(any(is.na(from)))
	      warning("'NA's coerced to 'FALSE' in coercion to logical sparse")
	  ij <- which(from, arr.ind = TRUE) - 1:1
	  if(length(ij) == 0) ij <- matrix(ij, 0, 2)
	  new("ngTMatrix",
	      i = ij[,1],
	      j = ij[,2],
	      Dim = as.integer(dim(from)),
	      Dimnames = .M.DN(from))
	  })

setAs("ngTMatrix", "dgTMatrix",
      function(from)
      ## more efficient than
      ## as(as(as(sM, "ngCMatrix"), "dgCMatrix"), "dgTMatrix")
      new("dgTMatrix", i = from@i, j = from@j,
          x = rep.int(1, length(from@i)),
          ## cannot copy factors, but can we use them?
          Dim = from@Dim, Dimnames= from@Dimnames))

setMethod("t", signature(x = "ngTMatrix"),
	  function(x) new("ngTMatrix", i = x@j, j = x@i,
			  Dim = x@Dim[2:1],
			  Dimnames= x@Dimnames[2:1]),
	  valueClass = "ngTMatrix")
