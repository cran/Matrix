 ### Coercion and Methods for Symmetric Matrices

setAs("dsyMatrix", "dgeMatrix",
      function(from) .Call("dsyMatrix_as_dgeMatrix", from) )

## I can't get this to work - at least inside Namespace -- FIXME
## setIs("dgeMatrix", "dsyMatrix",
##       ## R BUG:  test() doesn't see Matrix-internal functions
##       test = function(from) Matrix:::isSymmetric(from),
##       replace = function(obj, value) ## copy all slots
##       for(n in slotNames(obj)) slot(obj, n) <- slot(value, n)
##       )

setAs("dsyMatrix", "matrix",
      function(from) .Call("dsyMatrix_as_matrix", from) )

setAs("dsyMatrix", "dspMatrix",
      function(from) .Call("dsyMatrix_as_dspMatrix", from) )

setAs("dsyMatrix", "dsTMatrix",
      function(from) {
          ## This is not very efficient (FIXME)
          ij <- which(as(from,"matrix") != 0, arr.ind = TRUE)
          new("dsTMatrix", i = ij[,1], j = ij[,2],
              Dim = from@Dim, Dimnames = from@Dimnames)
      })
setAs("dsyMatrix", "dsCMatrix",
      function(from) callGeneric(as(from, "dsTMatrix")))


## Note: Just *because* we have an explicit  dtr -> dge coercion,
##       show( <ddenseMatrix> ) is not okay, and we need our own:
setMethod("show", "dsyMatrix", function(object) prMatrix(object))


setMethod("rcond", signature(x = "dsyMatrix", type = "character"),
          function(x, type, ...)
          .Call("dsyMatrix_rcond", x, type),
          valueClass = "numeric")

setMethod("rcond", signature(x = "dsyMatrix", type = "missing"),
          function(x, type, ...)
          .Call("dsyMatrix_rcond", x, "O"),
          valueClass = "numeric")

setMethod("%*%", signature(x = "dsyMatrix", y = "dgeMatrix"),
          function(x, y) .Call("dsyMatrix_dgeMatrix_mm", x, y) )

setMethod("%*%", signature(x = "dgeMatrix", y = "dsyMatrix"),
          function(x, y) .Call("dsyMatrix_dgeMatrix_mm_R", y, x) )

setMethod("solve", signature(a = "dsyMatrix", b = "missing"),
          function(a, b, ...) .Call("dsyMatrix_solve", a),
          valueClass = "dsyMatrix")

setMethod("solve", signature(a = "dsyMatrix", b = "matrix"),
          function(a, b, ...)
          .Call("dsyMatrix_matrix_solve", a, b),
          valueClass = "matrix")

setMethod("solve", signature(a = "dsyMatrix", b = "dgeMatrix"),
          function(a, b, ...)
          .Call("dsyMatrix_dgeMatrix_solve", a, b),
          valueClass = "dgeMatrix")

setMethod("norm", signature(x = "dsyMatrix", type = "character"),
          function(x, type, ...) .Call("dsyMatrix_norm", x, type),
          valueClass = "numeric")

setMethod("norm", signature(x = "dsyMatrix", type = "missing"),
          function(x, type, ...) .Call("dsyMatrix_norm", x, "O"),
          valueClass = "numeric")

## Should this create the opposite storage format - i.e. "U" -> "L"
## and vice-versa?
## MM: I think yes, since the other part can be filled arbitrarily (wrongly)
##WAS setMethod("t", signature(x = "dsyMatrix"), function(x) x)
setMethod("t", signature(x = "dsyMatrix"), t_trMatrix,
          valueClass = "dsyMatrix")

## The following has the severe effect of making
## "dsyMatrix" a subclass of "dpoMatrix" and since the reverse is
## by definition of "dpoMatrix", the class-hierarchy gets a *cycle* !
##
setIs("dsyMatrix", "dpoMatrix",
      test = function(obj)
          "try-error" != class(try(.Call("dpoMatrix_chol", obj), silent=TRUE)),
      ## MM: The following copying is necessary
      ## -- but shouldn't it be the default in such a case ??
      replace = function(obj, value) ## copy all slots
      for(n in slotNames(obj)) slot(obj, n) <- slot(value, n)
      )

## Now that we have "chol", we can define  "determinant" methods,
## exactly like in ./dsCMatrix.R
## DB - Probably figure out how to use the BunchKaufman decomposition instead
## {{FIXME: Shouldn't it be possible to have "determinant" work by
## default automatically for "Matrix"es  when there's a "chol" method available?
## -- not have to define showMethod("determinant", ...) for all classes

