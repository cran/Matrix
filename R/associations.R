## From: Michael Cysouw <cysouw@uni-marburg.de>
## Subject: Re: Sparse version of "cor"
## Date: Tue, 18 Feb 2014 10:34:15 +0100
## To: Martin Maechler <maechler@stat.math.ethz.ch>

## ... in addition to  cor.sparse(), the following purely sparse
## association measures seem interesting :

##--------------------------------


# cosine similarity matrix between columns of X, Y
# different weightings of normalisations can be used

	# allow for different weighting functions: can also be defined externally!
	# IDF is taken as default
	idf <- function(x) { log( (1+dim(x)[2]) / (1+rowSums(x)) ) }
	# inverse square root
	isqrt <- function(x) { rowSums(x) ^ -0.5 }

	# allow for different normalisation functions: can also be defined externally!
	# Euclidean 2-norm is taken as default
	norm2 <- function(x) { colSums(x^2) ^ -0.5 }
	# Alternatively take 1-norm
	norm1 <- function(x) { colSums(x) ^ -1 }

cos.sparse <- function(X, Y = NULL, norm = norm2 , weight = idf) {

	X <- as(X,"dgCMatrix")
	if (!is.null(Y)) {
		Y <- as(Y,"dgCMatrix")
	}

	if (!is.null(weight)) {
		Wx <- Diagonal( x = match.fun(weight)(X) )
		X <- Wx %*% X
		if (!is.null(Y)) {
			Wy <- Diagonal( x = match.fun(weight)(Y) )
			Y <- Wy %*% Y
		}
	}

	if (!is.null(norm)) {
		Nx <- Diagonal( x = match.fun(norm)(X) )
		X <- X %*% Nx
		if (!is.null(Y)) {
			Ny <- Diagonal( x = match.fun(norm)(Y) )
			Y <- Y %*% Ny
		}
	}

	if (is.null(Y)) {
		return(crossprod(X))
		} else {
		return(crossprod(X,Y))
	}
}

# association matrix between columns of X, Y. The argument "n" is to specify the N for the calculation of the expectation. Default, this is the number of columns of the matrix, but I have had situations (in combination with KhatriRao!) in which I needed different values, so I added the option.

	# allow for different functions to be specified: can be defined externally!
	# poisson-based association measure
	poisson <- function(o,e) { sign(o-e) * (o * log(o/e) - (o-e)) }
	# pointwise mutual information, aka "log-odds" in bioinformatics
	pmi <- function(o,e) { log(o/e) }
	# good old pearson residuals
	residuals <- function(o,e) { (o-e) / sqrt(e) }

assoc.sparse <- function(X, Y = NULL, method = "pmi", n = dim(X)[1]) {

	# observed coocurrences O
	X <- as(X,"dgCMatrix")
	if (is.null(Y)) {
		O <- crossprod(X)
            } else {
		Y <- as(Y,"dgCMatrix")
		O <- crossprod(X,Y)
	}

	# column sums as diagonal matrices
	Fx <- Diagonal( x = colSums(X) )
	if (is.null(Y)) {
		Fy <- Fx
		} else {
		Fy <- Diagonal( x = colSums(Y) )
	}

	# expected cooccurrences E only for cells with Observed!=0
	P <- as(O,"nMatrix")/n
	E <- Fx %*% P %*% Fy
	if (is.null(Y)) {
		E <- as(E,"symmetricMatrix")
	}

	P@x <- match.fun(method)(O@x,E@x)
	return(P)
}
