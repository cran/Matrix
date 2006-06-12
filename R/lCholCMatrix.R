                                        # Symbolic LDL' decompositions
                                        # of sparse matrices in
                                        # compressed, column-oriented format

setMethod("t", signature(x = "lCholCMatrix"),
          function(x) {
              as(x, "ltCMatrix") <- t(as(x, "ltCMatrix"))
              x
          },
          valueClass = "lCholCMatrix")

setMethod("solve", signature(a = "lCholCMatrix", b = "missing"),
          function(a, b)
          .Call(lCholCMatrix_solve, a),
          valueClass = "ltCMatrix")

setMethod("solve", signature(a = "lCholCMatrix", b = "lgCMatrix"),
          function(a, b)
          .Call(lCholCMatrix_lgCMatrix_solve, a, b),
          valueClass = "lgCMatrix")

