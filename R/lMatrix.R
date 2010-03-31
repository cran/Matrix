setAs("matrix", "lMatrix",
      function(from) { storage.mode(from) <- "logical" ; Matrix(from) })

## NOTE: This is *VERY* parallel to  ("dMatrix" -> "nMatrix") in ./dMatrix.R :
setAs("lMatrix", "nMatrix",
      function(from) {
	  if(any(is.na(from@x)))
	      stop("\"lMatrix\" object with NAs cannot be coerced to \"nMatrix\"")
	  ## i.e. from@x are only TRUE or FALSE
	  cld <- getClassDef(cl <- MatrixClass(class(from)))
	  if(extends(cld, "diagonalMatrix")) { # have no "ndi*" etc class
	      cl <- class(from <- as(from, "sparseMatrix"))
	      isSp <- TRUE
	  } else {
	      isSp <- extends(cld, "sparseMatrix")
	      if(isSp && !all(from@x)) {
		  from <- drop0(from) # was drop0(from, cld)
		  if(cl != (c. <- class(from)))
		      cld <- getClassDef(cl <- c.)
	      }
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
	  cld <- getClassDef(newCl <- class2(cl <- class(from), "l"))
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


## "ldi: is both "sparseMatrix" and "lMatrix" but not "lsparseMatrix"
setMethod("which", "ldiMatrix",
	  function(x, arr.ind) {
	      n <- x@Dim[1L]
	      i <- if(x@diag == "U") seq_len(n) else which(x@x)
	      if(arr.ind) cbind(i,i, deparse.level = 0) else i + n*(i - 1L) })

## "FIXME": use .Internal(which(x@x)) and arrayInd() for *dense*
setMethod("which", "nMatrix",
	  function(x, arr.ind) which(as(x, "sparseMatrix"), arr.ind=arr.ind))
setMethod("which", "lMatrix",
	  function(x, arr.ind) which(as(x, "sparseMatrix"), arr.ind=arr.ind))

.which.via.vec <- function(m) sort.int(as(m, "sparseVector")@i, method="quick")

setMethod("which", "nsparseMatrix",
	  function(x, arr.ind) {
	      if(arr.ind) which(as(x, "TsparseMatrix"), arr.ind=TRUE)
	      else .which.via.vec(x)
	  })
setMethod("which", "lsparseMatrix",
	  function(x, arr.ind) {
	      if(arr.ind) which(as(x, "TsparseMatrix"), arr.ind=TRUE)
	      else which(as(x, "sparseVector"))
	  })

which.ngT <- function(x, arr.ind)
    if(arr.ind) cbind(x@i, x@j) + 1L else .which.via.vec(x)
setMethod("which", "ngTMatrix", which.ngT)
setMethod("which", "ntTMatrix", function(x, arr.ind)
	  which.ngT(.Call(Tsparse_diagU2N, x), arr.ind))
setMethod("which", "nsTMatrix", function(x, arr.ind)
	  which.ngT(as(x, "generalMatrix"), arr.ind))

which.lgT <- function(x, arr.ind) {
    iT <- is1(x@x)
    if(arr.ind) cbind(x@i[iT], x@j[iT]) + 1L
    else sort.int(as(x, "sparseVector")@i[iT], method="quick")
}
setMethod("which", "lgTMatrix", which.lgT)
setMethod("which", "ltTMatrix", function(x, arr.ind)
	  which.lgT(.Call(Tsparse_diagU2N, x), arr.ind))
setMethod("which", "lsTMatrix", function(x, arr.ind)
	  which.lgT(as(x, "generalMatrix"), arr.ind))




## all() methods ---> ldenseMatrix.R and lsparseMatrix.R

setMethod("any", signature(x = "lMatrix"),
	  function(x, ..., na.rm = FALSE)
	  ## logical unit-triangular has TRUE diagonal:
	  (prod(dim(x)) >= 1 && is(x, "triangularMatrix") && x@diag == "U") ||
	  any(x@x, ..., na.rm = na.rm))
