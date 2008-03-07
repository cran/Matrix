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
{
    if(length(list(...)))
	warning("arguments in", deparse(list(...)), "are disregarded")
    sysDef <- eval(formals()$system)
    .Call(CHMfactor_solve, a, b,
	  match(match.arg(system, sysDef), sysDef, nomatch = 0))
}

setMethod("solve", signature(a = "CHMfactor", b = "ddenseMatrix"),
	  .CHM_solve, valueClass = "dgeMatrix")

setMethod("solve", signature(a = "CHMfactor", b = "matrix"),
	  .CHM_solve, valueClass = "dgeMatrix")

setMethod("solve", signature(a = "CHMfactor", b = "numeric"),
	  function(a, b, ...)
	  .CHM_solve(a, matrix(if(is.double(b)) b else as.double(b),
			       length(b), 1L), ...),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "CHMfactor", b = "dsparseMatrix"),
	  function(a, b,
		   system = c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"),
		   ...) {
	      if(length(list(...)))
		  warning("arguments in", deparse(list(...)), "are disregarded")
	      sysDef <- eval(formals()$system)
	      .Call(CHMfactor_spsolve, a, as(b, "dgCMatrix"),
		    match(match.arg(system, sysDef), sysDef, nomatch = 0))
	  }, valueClass = "dgCMatrix")

## Catch-all the rest : make sure 'system' is not lost
setMethod("solve", signature(a = "CHMfactor", b = "ANY"),
	  function(a, b, system = c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"),
		   ...)
	      solve(a, as(b, "dMatrix"), system, ...))
