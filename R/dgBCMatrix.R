setAs("dgBCMatrix", "dgTMatrix",
      function(from) .Call("dgBCMatrix_to_dgTMatrix", from))

setAs("dgBCMatrix", "dgCMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgCMatrix"))

setAs("dgBCMatrix", "dgeMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgeMatrix"))

setAs("dgBCMatrix", "matrix",
      function(from) as(as(from, "dgTMatrix"), "matrix"))
