#### "Namespace private" Auxiliaries  such as method functions
#### (called from more than one place --> need to be defined early)

## For %*% (M = Matrix; v = vector (double or integer {complex maybe?}):
.M.v <- function(x, y) callGeneric(x, as.matrix(y))
.v.M <- function(x, y) callGeneric(rbind(x), y)

.has.DN <- ## has non-trivial Dimnames slot?
    function(x) !identical(list(NULL,NULL), x@Dimnames)

## chol() via "dpoMatrix"
cholMat <- function(x, pivot, LINPACK) {
    px <- as(x, "dpoMatrix")
    if(identical(TRUE, validObject(px, test=TRUE)))
        chol(px)
    else stop("'x' is not positive definite -- chol() undefined.")
}
