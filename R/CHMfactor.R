setAs("CHMfactor", "sparseMatrix",
      function(from) .Call(CHMfactor_to_sparse, from))

setMethod("image", "CHMfactor",
          function(x, ...) {
              x <- as(as(x, "sparseMatrix"), "dgTMatrix")
              callGeneric()
          })

setMethod("solve", signature(a = "CHMfactor", b = "ddenseMatrix"),
          function(a, b,
                   system = c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"),
                   ...)
          .Call(CHMfactor_solve, a, b,
                match(match.arg(system),
                      c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"),
                      nomatch = 0)),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "CHMfactor", b = "matrix"),
          function(a, b,
                   system = c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"),
                   ...)
          .Call(CHMfactor_solve, a, b,
                match(match.arg(system),
                      c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"),
                      nomatch = 0)),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "CHMfactor", b = "dsparseMatrix"),
          function(a, b,
                   system = c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"),
                   ...)
          .Call(CHMfactor_spsolve, a, as(b, "dgCMatrix"),
                match(match.arg(system),
                      c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"),
                      nomatch = 0)),
          valueClass = "dgCMatrix")
