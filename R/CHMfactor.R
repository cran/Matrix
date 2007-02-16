setAs("CHMfactor", "sparseMatrix",
      function(from) .Call(CHMfactor_to_sparse, from))

setMethod("image", "CHMfactor",
          function(x, ...) {
              x <- as(as(x, "sparseMatrix"), "dgTMatrix")
              callGeneric()
          })
