####-----------  Minimal conversion utilities  <-->  "SparseM"

### I.  The  "natural pairs"  between the two packages:

setAs("matrix.csr", "dgRMatrix",
      function(from) {
	  new("dgRMatrix",
	      x = from@ra, j = from@ja - 1:1, p = from@ia - 1:1,
	      Dim = from@dimension)
      })
setAs("dgRMatrix", "matrix.csr",
      function(from) {
	  new("matrix.csr",
	      ra = from@x, ja = from@j + 1:1, ia = from@p + 1:1,
	      dimension = from@Dim)
      })



setAs("matrix.csc", "dgCMatrix",
      function(from) {
	  new("dgCMatrix",
	      x = from@ra, i = from@ja - 1:1, p = from@ia - 1:1,
	      Dim = from@dimension)
      })
setAs( "dgCMatrix", "matrix.csc",
      function(from) {
	  new("dgCMatrix",
	      ra = from@x, ja = from@i + 1:1, ia = from@p + 1:1,
	      Dim = from@dimension)
      })

setAs("matrix.coo", "dgTMatrix",
      function(from) {
	  new("dgTMatrix",
	      x = from@ra, i = from@ia - 1:1, j = from@ja - 1:1,
	      Dim = from@dimension)
      })
setAs("dgTMatrix", "matrix.coo",
      function(from) {
	  new("dgTMatrix",
	      ra = from@x, ia = from@i + 1:1, ja = from@j + 1:1,
	      Dim = from@dimension)
      })

### II.  Enable coercion to the ``favorite'' of each package;
### ---         ----------------------------
###      i.e.,  "dgCMatrix" and  "matrix.csr"

setAs("dsparseMatrix", "matrix.csr",
      function(from) as(as(from, "dgRMatrix"), "matrix.csr"))

##
setAs("matrix.csr", "dgCMatrix",
      function(from) as(as(from, "dgRMatrix"), "dgCMatrix"))
setAs("matrix.coo", "dgCMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgCMatrix"))

## Easy coercion: just always use as( <SparseM.mat>, "Matrix") :

setAs("matrix.csr", "Matrix", function(from) as(from, "dgCMatrix")) # we favor 'dgC' so much ..
setAs("matrix.coo", "Matrix", function(from) as(from, "dgTMatrix"))
setAs("matrix.csc", "Matrix", function(from) as(from, "dgCMatrix"))


