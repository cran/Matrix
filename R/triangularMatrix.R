#### Methods for the virtual class 'triangularMatrix' of triangular matrices
#### Note that specific methods are in (8 different) ./?t?Matrix.R

.tril.tr <- function(x, k = 0, ...) {  # are always square
    k <- as.integer(k[1])
    dd <- dim(x)
    stopifnot(-dd[1] <= k, k <= dd[1])  # had k <= 0
    if(k <= 0 && x@uplo == "L")
        x
    else ## more to do
        callNextMethod()
}

.triu.tr <- function(x, k = 0, ...) {  # are always square
    k <- as.integer(k[1])
    dd <- dim(x)
    stopifnot(-dd[1] <= k, k <= dd[1])  # had k >= 0
    if(k >= 0 && x@uplo == "U")
        x
    else ## more to do
        callNextMethod()
}

## In order to evade method dispatch ambiguity (with [CTR]sparse* and ddense*),
## but still remain "general"
## we use this hack instead of signature  x = "triagonalMatrix" :

for(cls in names(getClass("triangularMatrix")@subclasses))
    if(length(grep(".t.Matrix", cls)) == 1) { # not for "Cholesky"
        setMethod("tril", cls, .tril.tr)
        setMethod("triu", cls, .triu.tr)
    }


