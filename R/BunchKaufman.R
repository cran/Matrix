## METHODS FOR GENERIC: BunchKaufman
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("BunchKaufman", signature(x = "dsyMatrix"),
	  function(x, ...) .Call(dsyMatrix_trf, x, 2L))

setMethod("BunchKaufman", signature(x = "dspMatrix"),
	  function(x, ...) .Call(dspMatrix_trf, x, 2L))

setMethod("BunchKaufman", signature(x = "matrix"),
	  function(x, uplo = NULL, ...) .Call(matrix_trf, x, 2L, uplo))


## METHODS FOR CLASS: p?BunchKaufman
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(FALSE) {
## returning:
##
## list(D, P[n], U[n], ..., P[1], U[n])
##     where A = U' D U and U = P[n] U[n] ... P[1] U[1]
##
## OR
##
## list(P[1], L[1], ..., P[n], L[n], D)
##     where A = L D L' and L = P[1] L[1] ... P[n] L[n]
##
## as described in the documentation for LAPACK 'ds[yp]trf'
setMethod("expand", signature(x = "BunchKaufman"),
          function(x, ...) .Call(BunchKaufman_expand, x))

setMethod("expand", signature(x = "pBunchKaufman"),
          function(x, ...) .Call(BunchKaufman_expand, x))
}

if(FALSE) {
library(Matrix)
set.seed(1)

X <- new("dsyMatrix", Dim = c(6L, 6L), x = rnorm(36L))
Y <- t(X)

as(bkX <- BunchKaufman(X), "dtrMatrix")
as(bkY <- BunchKaufman(Y), "dtrMatrix")

DU <- .Call("BunchKaufman_expand", bkX)
D <- DU[[1L]]
U <- Reduce(`%*%`, DU[-1L])
## FIXME: 'DU' looks correct ... but is actually wrong {second test fails}??
stopifnot(identical(DU, .Call("BunchKaufman_expand", pack(bkX))),
          all.equal(as(t(U) %*% D %*% U, "matrix"), as(X, "matrix")))

LD <- .Call("BunchKaufman_expand", bkY)
D <- LD[[length(LD)]]
L <- Reduce(`%*%`, LD[-length(LD)])
stopifnot(identical(LD, .Call("BunchKaufman_expand", pack(bkY))),
          all.equal(as(L %*% D %*% t(L), "matrix"), as(Y, "matrix")))
}
