## METHODS FOR CLASS: dsyMatrix
## dense (unpacked) symmetric matrices with 'x' slot of type "double"
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: no longer needed ... replacement in ./symmetricMatrix.R, ./denseMatrix.R
##     { where .dense2sy(), .dsy2mat() and .dsy2dsp(), which
##     have all been exported, are defined as simple aliases ... }
if(FALSE) {
.dense2sy <- function(from, ...) {
    if(isSymmetric(from, ...)) # < with tolerance!
	.Call(dense_to_symmetric, from, "U", FALSE)
    else
	stop("not a symmetric matrix; consider forceSymmetric() or symmpart()")
}
..dense2sy <- function(from) {
    ## NB: The alternative, 'zero tolerance' { <=> isSymmetric(*, tol=0) }
    ##     breaks too much previous code -- though it would be much faster --
    if(isSymmetric(from)) # < with tolerance!
	.Call(dense_to_symmetric, from, "U", FALSE)
    else
	stop("not a symmetric matrix; consider forceSymmetric() or symmpart()")
}
.dsy2mat <- function(from, keep.dimnames = TRUE) {
    .Call(dsyMatrix_as_matrix, from, keep.dimnames)
}
..dsy2mat <- function(from) {
    .Call(dsyMatrix_as_matrix, from, TRUE)
}
.dsy2dsp <- function(from) {
    .Call(dsyMatrix_as_dspMatrix, from)
}

setAs("dgeMatrix", "dsyMatrix", ..dense2sy)
setAs(   "matrix", "dsyMatrix", function(from) .dense2sy(..2dge(from)))
setAs("dsyMatrix",    "matrix", ..dsy2mat)
setAs("dsyMatrix", "dspMatrix", .dsy2dsp)

dsy2T <- function(from) { # 'dsT': only store upper *or* lower
    uplo <- from@uplo
    if(any0(dim(from))) {
	ij <- matrix(0L, 0,2) ; m <- from@x
    } else {
	## FIXME!	 working via "matrix" is *not* efficient:
	## the "other triangle" is filled, compared with 0, and then trashed:
	m <- .dense2m(from)
	ij <- which(m != 0, arr.ind = TRUE, useNames = FALSE)
	ij <- ij[if(uplo == "U") ij[,1] <= ij[,2] else ij[,1] >= ij[,2], , drop = FALSE]
    }
    new("dsTMatrix", i = ij[,1] - 1L, j = ij[,2] - 1L,
	x = as.vector(m[ij]), uplo = uplo,
	Dim = from@Dim, Dimnames = from@Dimnames)
}
dsy2C <- function(from) .T2Cmat(dsy2T(from), isTri=FALSE)

setAs("dsyMatrix", "dsTMatrix", dsy2T)
setAs("dsyMatrix", "dsCMatrix", dsy2C)
} ## MJ

## MJ: no longer needed ... replacement in ./unpackedMatrix.R
if(FALSE) {
.dsy.diag <- function(x, nrow, ncol, names=TRUE) {
    if(min(dim(x)) == 0L) return(numeric(0L))
    y <- .Call(dgeMatrix_getDiag, x)
    if(names) {
        nms <- dimnames(x) # rely on method for "symmetricMatrix" to symmetrize
        if(is.list(nms) && length(nms) == 2L)
            names(y) <- nms[[1L]]
    }
    y
}
## *Should* create the opposite storage format:  "U" -> "L"  and vice-versa:
setMethod("t", signature(x = "dsyMatrix"), t_trMatrix,
          valueClass = "dsyMatrix")
setMethod("diag", signature(x = "dsyMatrix"), .dsy.diag)
setMethod("diag<-", signature(x = "dsyMatrix"),
	  function(x, value) .Call(dgeMatrix_setDiag, x, value))
} ## MJ
