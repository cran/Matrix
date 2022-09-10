## METHODS FOR CLASS: dgeMatrix
## dense general matrices with 'x' slot of type "double"
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
## ..2dge() -> ./Auxiliaries.R
setAs("matrix",  "dgeMatrix", function(from) ..2dge(from))
setAs("numLike", "dgeMatrix", function(from) ..2dge(from))

ge2mat <- function(from) array(from@x, dim = from@Dim, dimnames = from@Dimnames)
setAs("dgeMatrix", "matrix", ge2mat)

setMethod("as.vector", "dgeMatrix",
          function(x, mode) as.vector(x@x, mode))
} ## MJ

## MJ: no longer needed ... replacement in ./unpackedMatrix.R
if (FALSE) {
..get.diag <- function(x, nrow, ncol, names=TRUE) {
    ##         vvvvvvvvvvvvvvvvv here just a place holder, replaced in .mkSpec.diag()
    y <- .Call(dgeMatrix_getDiag, x) # double or logical
    if(names) {
        nms <- dimnames(x)
        if(is.list(nms) && !any(vapply(nms, is.null, NA)) &&
           identical((nm <- nms[[1L]][im <- seq_len(min(dim(x)))]), nms[[2L]][im]))
            names(y) <- nm
    }
    y
}
.mkSpec.diag <- function(symb) {
    rr <- ..get.diag
    body(rr)[[2]][[3]][[2]] <- symb
    rr
}
.dge.diag <- .mkSpec.diag(quote(dgeMatrix_getDiag))

setMethod("t", signature(x = "dgeMatrix"), t_geMatrix)
setMethod("diag", signature(x = "dgeMatrix"), .dge.diag)
setMethod("diag<-", signature(x = "dgeMatrix"),
	  function(x, value) .Call(dgeMatrix_setDiag, x, value))
} ## MJ
