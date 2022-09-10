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

KhatriRao <- function(X, Y = X, FUN = "*",
                      sparseY = TRUE, make.dimnames = FALSE) {
    stopifnot((p <- ncol(X)) == ncol(Y))
    ## TODO? speedup when X = Diagonal(<n>)
    X <- asCspN(X) # not-U-diag CsparseMatrix -->  can use @x, @p, @i
    ## xn, yn; number of non-zero values per column
    xn <- diff(X@p)
    n1 <- nrow(X); n2 <- nrow(Y)
    if(sparseY) {
        Y <- asCspN(Y)
        yn <- diff(yp <- Y@p) ## both of length p
    } else {
        if(!is.matrix(Y)) Y <- as.matrix(Y)# needed e.g. in c(Y) below)
        yn <- rep(n2, p)
    }
    newp <- as.integer(diffinv(xn*yn))
    ## indices for new values
    ## dense: = newp <- as.integer(diffinv(xn*nrow(Y)))
    xn.yp <- xn[ as.logical(yn) ] # xn "where" Y is present
    yj <- if(sparseY)
              .Call(Matrix_expand_pointers, yp)## as(Y,"TsparseMatrix")@j
          else rep(seq(p)-1L, each = n2)
    yj <- factor(yj) # for split() below
    non0 <- length(xn.yp) > 0L && any(xn.yp != 0L)
    rep.yn <- rep.int(yn, xn)
    i1 <- rep.int(X@i, rep.yn)
    i2 <-
	if(non0) {
            if(sparseY)
                unlist(rep(split.default(Y@i,yj), xn.yp), use.names=FALSE)
            else rep.int(seq.int(n2)-1L, p)
	}
	else integer()
    newi <- i1*n2 + i2
    dim <- as.integer(c(n1*n2, p))

    dns <- if (make.dimnames) { ## this is not good enough:  dnx, dny may be NULL
	list(as.vector(outer(rownames(Y),rownames(X), FUN = "paste", sep = ":")),
	     colnames(X))
    } else list(NULL,NULL)

    if((nX <- is(X, "nMatrix")) & (nY <- is(Y, "nMatrix")))
	new("ngCMatrix", Dim=dim, Dimnames=dns, i = newi, p = newp)
    else { ## at least one of 'X' and 'Y' has an "x" slot:
	if(nX) X <- .sparse2g(..sparse2l(X))
	x1 <- rep.int(X@x, rep.yn)
	x2 <- if(non0) {
		  if(nY && sparseY) Y <- .sparse2g(..sparse2l(X))
                  yx <- if(sparseY) Y@x else c(Y)
		  unlist(rep(split.default(yx, yj), xn.yp), use.names=FALSE)
	      } else if(nY) logical() else (if(sparseY) Y@x else c(Y))[0]
        if(!is.double(x1)) x1 <- as.double(x1) ## or if we had "igCMatrix" ...
	new("dgCMatrix", Dim=dim, Dimnames=dns, i = newi, p = newp,
	    x = match.fun(FUN) (x1,x2))
    }
}
