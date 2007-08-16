setAs("CHMfactor", "sparseMatrix",
      function(from) .Call(CHMfactor_to_sparse, from))

setAs("CHMfactor", "Matrix", function(from) as(from, "sparseMatrix"))

setMethod("image", "CHMfactor",
          function(x, ...) {
              x <- as(as(x, "sparseMatrix"), "dgTMatrix")
              callGeneric()
          })

.CHM_solve <-
    function(a, b,
             system = c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"),
             ...)
    .Call(CHMfactor_solve, a, b,
          match(match.arg(system),
                c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"),
                nomatch = 0))

setMethod("solve", signature(a = "CHMfactor", b = "ddenseMatrix"),
          .CHM_solve, valueClass = "dgeMatrix")

setMethod("solve", signature(a = "CHMfactor", b = "matrix"),
          .CHM_solve, valueClass = "dgeMatrix")

setMethod("solve", signature(a = "CHMfactor", b = "dsparseMatrix"),
          function(a, b,
                   system = c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"),
                   ...)
          .Call(CHMfactor_spsolve, a, as(b, "dgCMatrix"),
                match(match.arg(system),
                      c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"),
                      nomatch = 0)),
          valueClass = "dgCMatrix")
