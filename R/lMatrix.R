## MJ: no longer ... prefer more efficient methods defined for denseMatrix,
##     .sparseMatrix, and diagonalMatrix separately, going via C utilities
##     R_(dense|sparse)_as_kind() and R_diagonal_as_sparse()
if(FALSE) {
## NOTE: This is *VERY* parallel to  ("dMatrix" -> "nMatrix") in ./dMatrix.R :
setAs("lMatrix", "nMatrix",
      function(from) {
	  if(anyNA(from@x) && ((.w <- isTRUE(getOption("Matrix.warn"))) ||
				   isTRUE(getOption("Matrix.verbose")))) {
	      (if(.w) warning else message)(
		  "\"lMatrix\" object with NAs coerced to \"nMatrix\":  NA |-> TRUE")
	      from@x[is.na(from@x)] <- TRUE
	  }
	  ## ==> from@x are in {TRUE, FALSE}
	  cld <- getClassDef(cl <- MatrixClass(class(from)))
	  if(extends(cld, "diagonalMatrix")) # no "ndi*" class
	      ## should not happen, setAs(diagonalMatrix -> nMatrix) in ./diagMatrix.R:
	      return(di2nMat(from))
	  ## else
	  isSp <- extends(cld, "sparseMatrix")
	  if(isSp && !all(from@x)) {
	      from <- drop0(from) # was drop0(from, cld)
	      if(cl != (c. <- class(from)))
		  cld <- getClassDef(cl <- c.)
	  }
	  sNams <- slotNames(cld)
	  copyClass(from, sub("^l", "n", cl),
		    if(isSp) sNams[sNams != "x"] else sNams)
      })

## and the reverse as well :

setAs("nMatrix", "lMatrix",
      function(from) {
	  cld <- getClassDef(cl <- MatrixClass(class(from)))
	  r <- copyClass(from, sub("^n", "l", cl), slotNames(cld))
	  if(extends(cld, "sparseMatrix"))
	      r@x <- rep.int(TRUE, length(if(!extends(cld, "RsparseMatrix"))
					  from@i else from@j))
	  r
      })

setAs("dMatrix", "lMatrix",
      function(from) {
	  cld <- getClassDef(newCl <- class2(class(from), "l"))
	  sNams <- slotNames(cld)
	  r <- copyClass(from, newCl, sNames = sNams[sNams != "x"])
	  r@x <- as.logical(from@x)
	  r
      })

setAs("lMatrix", "dMatrix",
      function(from) {
	  cld <- getClassDef(cl <- MatrixClass(class(from)))
	  sNams <- slotNames(cld)
	  r <- copyClass(from, newCl = sub("^l", "d", cl),
			 sNames = sNams[sNams != "x"])
	  r@x <- as.double(from@x)
	  r
      })

## needed at least for lsparse* :
setAs("lMatrix", "dgCMatrix",
      function(from) as(as(from, "lgCMatrix"), "dgCMatrix"))
} ## MJ
