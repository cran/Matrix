#### Sparse Matrices in Compressed column-oriented format

### contains = "dsparseMatrix", "CsparseMatrix"

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
setAs("dgCMatrix", "ngCMatrix", function(from) .C2nC(from, FALSE))

## rather use Csparse* to lsparse* in ./lsparseMatrix.R ,
## but this is for "back-compatibility" (have had tests for it..):
setAs("dgCMatrix", "lgCMatrix",
      function(from) { ## FIXME use .Call() too!
	  r <- new("lgCMatrix")
	  r@x <- as.logical(from@x)
	  ## and copy the other slots
	  for(nm in c("i", "p", "Dim", "Dimnames"))
	      slot(r, nm) <- slot(from, nm)
	  r
      })
} ## MJ

## MJ: no longer needed ... now inherited from Matrix
if(FALSE) {
setMethod("image", "dgCMatrix", function(x, ...) image(as(x, "dgTMatrix"), ...))
} ## MJ

## Group Methods, see ?Arith (e.g.)
## -----
##
## "Arith" is now in ./Ops.R
##
## "Math" and "Math2"  in ./Math.R



## "[<-" methods { setReplaceMethod()s }  are now in ./Csparse.R

## setMethod("writeHB", signature(obj = "dgCMatrix"),
## 	  function(obj, file, ...) {
## 	      .Deprecated("writeMM")
## 	      .Call(Matrix_writeHarwellBoeing, obj,
## 		    as.character(file), "DGC")
## 	  })

##-> ./colSums.R  for colSums,... rowMeans

## MJ: no longer needed ... now inherited from CsparseMatrix
if(FALSE) {
setMethod("t", signature(x = "dgCMatrix"),
	  function(x) .Call(Csparse_transpose, x, FALSE),
	  valueClass = "dgCMatrix")
} ## MJ
