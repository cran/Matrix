#### Superclass Methods for all sparse nonzero-pattern matrices

## Would be nice to do this, but it also triggers for e.g.
## ngCMatrix and hence infinite recursion
if(FALSE) ## "bug" in 'methods' -  R_FIXME ?
setAs("CsparseMatrix", "nsparseMatrix",
      function(from) .Call(Csparse_to_nz_pattern, from,
			   is(from, "triangularMatrix")))
setAs("CsparseMatrix", "nsparseMatrix",
      function(from) {
	  cld <- getClassDef(class(from))
	  if(extends(cld, "nsparseMatrix"))
	      from
	  else
	      .Call(Csparse_to_nz_pattern, from,
		    extends(cld, "triangularMatrix"))
      })


###------- Work via  as(*, lgC) : ------------

## For multiplication operations, sparseMatrix overrides other method
## selections.  Coerce a ddensematrix argument to a lsparseMatrix.

setMethod("%*%", signature(x = "nsparseMatrix", y = "ndenseMatrix"),
          function(x, y) callGeneric(x, as(y, "nsparseMatrix")))

setMethod("%*%", signature(x = "ndenseMatrix", y = "nsparseMatrix"),
          function(x, y) callGeneric(as(x, "nsparseMatrix"), y))

setMethod("crossprod", signature(x = "nsparseMatrix", y = "ndenseMatrix"),
          function(x, y = NULL) callGeneric(x, as(y, "nsparseMatrix")))

setMethod("crossprod", signature(x = "ndenseMatrix", y = "nsparseMatrix"),
          function(x, y = NULL) callGeneric(as(x, "nsparseMatrix"), y))

## and coerce lsparse* to lgC*
setMethod("%*%", signature(x = "nsparseMatrix", y = "nsparseMatrix"),
          function(x, y) callGeneric(as(x, "ngCMatrix"), as(y, "ngCMatrix")))

setMethod("crossprod", signature(x = "nsparseMatrix", y = "nsparseMatrix"),
          function(x, y = NULL)
          callGeneric(as(x, "ngCMatrix"), as(y, "ngCMatrix")))

## Use "Matrix" method !as(. , "lMatrix")
## setMethod("!", "nsparseMatrix",
##           ## turns FALSE to TRUE --> dense matrix
##           function(e1) !as(e1, "ngeMatrix"))
