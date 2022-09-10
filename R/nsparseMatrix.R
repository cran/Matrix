## METHODS FOR CLASS: lsparseMatrix (virtual)
## sparse matrices with no 'x' slot (nonzero pattern)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: no longer needed ... now inherited from Matrix
if(FALSE) {
setMethod("image", "nsparseMatrix",
          function(x, ...) image(as(x, "dMatrix"), ...))
} ## MJ

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
.C2nC <- function(from, isTri = is(from, "triangularMatrix"))
    .Call(Csparse_to_nz_pattern, from, isTri)

setAs("CsparseMatrix", "nsparseMatrix", function(from) .C2nC(from))
setAs("CsparseMatrix", "nMatrix",       function(from) .C2nC(from))

setAs("nsparseMatrix", "dsparseMatrix", function(from) as(from, "dMatrix"))
} ## MJ
