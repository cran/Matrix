# Efficient Khatri-Rao product for large sparse matrices
# Assumes two matrices in CsparseMatrix format
# Written by Michael Cysouw <cysouw@mac.com>

## MM: there's  a "public" Matlab version, at
## http://www.mathworks.com/matlabcentral/fileexchange/28872-khatri-rao-product/content/kr.m
## with documentation
##
## %  Khatri-Rao product.
##
## %   kr(A,B) returns the Khatri-Rao product of two matrices A and B, of
## %   dimensions I-by-K and J-by-K respectively. The result is an I*J-by-K
## %   matrix formed by the matching columnwise Kronecker products, i.e.
## %   the k-th column of the Khatri-Rao product is defined as
## %   kron(A(:,k),B(:,k)).


KhatriRao <- function(X, Y = X, FUN = "*", make.dimnames = FALSE)
{
    stopifnot((p <- ncol(X)) == ncol(Y))
    X <- as(X,"CsparseMatrix")
    Y <- as(Y,"CsparseMatrix")
    is.n <- (nX <- is(X, "nMatrix")) & (nY <- is(Y, "nMatrix"))

    xn <- diff(      X@p)
    yn <- diff(yp <- Y@p) ## both of length p
    newp <- as.integer(diffinv(xn*yn))


    rep.yn <- rep.int(yn,xn)
    xn.yp <- xn[ as.logical(yn) ] # xn "where" Y is present
    non0 <- length(xn.yp) > 0L && any(xn.yp != 0L)
    i1 <- rep.int(X@i, rep.yn)
    i2 <-
	if(non0) {
	    yj <- .Call(Matrix_expand_pointers, yp)## as(Y,"TsparseMatrix")@j
	    yj <- factor(yj) # for 2x split() below
	    unlist(rep(split.default(Y@i,yj), xn.yp))
	}
	else integer()
    n1 <- nrow(X); n2 <- nrow(Y)
    newi <- i1*n2 + i2
    dim <- as.integer(c(n1*n2, p))

    dns <- if (make.dimnames) { ## this is not good enough:  dnx, dny may be NULL
	list(as.vector(outer(rownames(Y),rownames(X), FUN = "paste", sep = ":")),
	     colnames(X))
    } else list(NULL,NULL)

    if(is.n)
	new("ngCMatrix", Dim=dim, Dimnames=dns, i = newi, p = newp)
    else { ## at least one of 'X' and 'Y' has an "x" slot:
	if(nX) X <- as(X, "lgCMatrix")
	x1 <- rep.int(X@x, rep.yn)
	x2 <- if(non0) {
		  if(nY) Y <- as(Y, "lgCMatrix")
		  unlist(rep(split.default(Y@x,yj), xn.yp))
	      } else if(nY) logical() else Y@x[0]
	new("dgCMatrix", Dim=dim, Dimnames=dns, i = newi, p = newp,
	    x = match.fun(FUN) (x1,x2))
    }
}
