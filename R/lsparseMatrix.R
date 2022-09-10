## METHODS FOR CLASS: lsparseMatrix (virtual)
## sparse matrices with 'x' slot of type "logical" (TRUE, FALSE, or NA)
## ... but _excluding_ ldiMatrix (FIXME?)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: no longer needed ... now inherited from Matrix
if(FALSE) {
setMethod("image", "lsparseMatrix",
          function(x, ...) image(as(x, "dMatrix"), ...))
} ## MJ

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
C2l <- function(from) {
    if(extends(cld <- getClassDef(class(from)), "lsparseMatrix"))
	return(from)
    ## else
    if(!(is.n <- extends(cld, "nsparseMatrix"))) {
        ## len.x <- length(from@x)
        from <- .Call(Csparse_drop, from, 0)
        ## did.drop <- length(from@x) != len.x
    }
    r <- as(.C2nC(from, extends(cld, "triangularMatrix")), "lsparseMatrix")
    if(!is.n && any(ina <- is.na(from@x))) { ## NAs must remain NA
        ## since we dropped, we "know"  that the 'x' slots match:
        stopifnot(length(from@x) == length(r@x))
        is.na(r@x) <- ina
    }
    r
}

setAs("CsparseMatrix", "lMatrix", C2l)
setAs("CsparseMatrix", "lsparseMatrix", C2l)

setAs("lsparseMatrix", "dsparseMatrix", function(from) as(from, "dMatrix"))

setAs("lsparseMatrix", "matrix",
      function(from) as(as(from, "ldenseMatrix"), "matrix"))
} ## MJ
