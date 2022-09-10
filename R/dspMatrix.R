## METHODS FOR CLASS: dspMatrix
## dense (packed) symmetric matrices with 'x' slot of type "double"
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
dsp2dsy <- function(from) .Call(dspMatrix_as_dsyMatrix, from)
dsp2C <- function(from) dsy2C(.Call(dspMatrix_as_dsyMatrix, from))
setAs("dspMatrix", "dsyMatrix", dsp2dsy)

## dge <--> dsp   via  dsy
.dense2sp <- function(from) .dsy2dsp(.dense2sy(from))
setAs("dgeMatrix", "dspMatrix", .dense2sp)
setAs("matrix", "dspMatrix",
      function(from) .dense2sp(..2dge(from)))

## S3-matrix <--> dsp   via  dsy
setAs("dspMatrix", "matrix", function(from) .dsy2mat(dsp2dsy(from)))

dsp2C <- function(from) dsy2C(unpack(from))
setAs("dspMatrix", "CsparseMatrix", dsp2C)
setAs("dspMatrix", "sparseMatrix", dsp2C)
} ## MJ

## MJ: no longer needed ... replacement in ./packedMatrix.R
if (FALSE) {
setMethod("t", signature(x = "dspMatrix"),
          function(x) .dsy2dsp(t(dsp2dsy(x))), # FIXME inefficient
          valueClass = "dspMatrix")

setMethod("diag", signature(x = "dspMatrix"),
	  function(x, nrow, ncol) .Call(dspMatrix_getDiag, x))
setMethod("diag<-", signature(x = "dspMatrix"),
	  function(x, value) .Call(dspMatrix_setDiag, x, value))
} ## MJ
