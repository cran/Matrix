
## nearcor.R :
## Copyright (2007) Jens Oehlschlägel
## GPL licence, no warranty, use at your own risk

nearPD <-
    ## Computes the nearest correlation matrix to an approximate
    ## correlation matrix, i.e. not positive semidefinite.

    function(x               # n-by-n approx covariance/correlation matrix
             , corr = FALSE
             , eig.tol   = 1e-6 # defines relative positiveness of eigenvalues compared to largest
             , conv.tol  = 1e-7 # convergence tolerance for algorithm
             , posd.tol  = 1e-8 # tolerance for enforcing positive definiteness
             , do2eigen = TRUE  # if TRUE do a sfsmisc::posdefify() eigen step
             , maxit    = 100 # maximum number of iterations allowed
             , verbose   = FALSE # set to TRUE to verbose convergence
             )
{
    stopifnot(isSymmetric(x))

    ## MM: these are there in 'Matrix' ! :
    ## Inf norm   == norm(*, "I")
    ## inorm <- function(x) max(rowSums(abs(x)))
    ## Froebenius norm == norm(*, "F")
    ## fnorm <- function(x) sqrt(sum(diag(t(x) %*% x)))

    n <- ncol(x)
    ## U should be like x, but filled with '0' -- following also works for 'Matrix':
    U <- x; U[] <- 0
    X <- x
    iter <- 0 ; converged <- FALSE; conv <- Inf

    while (iter < maxit && !converged) {
        Y <- X
        T <- Y - U

        ## project onto PSD matrices
        e <- eigen(Y, symmetric = TRUE)
        Q <- e$vectors
        d <- e$values
        D <- diag(d)

        ## create mask from relative positive eigenvalues
        p <-  d > eig.tol*d[1]

        ## MM{FIXME}: Make the following more efficient {do *NOT* multiply with diag!}
        ## use p mask to only compute 'positive' part
        X <- Q[,p,drop = FALSE] %*% D[p,p,drop = FALSE] %*% t(Q[,p,drop = FALSE])

        ## update Dykstra's correction
        U <- X - T

        ## project onto unit diag matrices
        X <- (X + t(X))/2
        if(corr) diag(X) <- 1

        conv <- norm(Y-X, "I") / norm(Y, "I")
        iter <- iter + 1
        if (verbose)
            cat(sprintf("iter %3d : conv.crit= %11g\n", iter, conv))

        converged <- (conv <= conv.tol)
    }

    if(!converged) {
        warning("nearPD() did not converge in ", iter, " iterations")
    }

    ## force symmetry
    X <- (X + t(X))/2
    if(do2eigen) {
        ## begin from posdefify(sfsmisc)
        e <- eigen(X, symmetric = TRUE)
        d <- e$values
        Eps <- posd.tol * abs(d[1])
        if (d[n] < Eps) {
            d[d < Eps] <- Eps
            Q <- e$vectors
            o.diag <- diag(X)
            X <- Q %*% (d * t(Q))
            D <- sqrt(pmax(Eps, o.diag)/diag(X))
            X[] <- D * X * rep(D, each = n)
        }
    } ## end from posdefify(sfsmisc)

    ## unneeded(?!): X <- (X + t(X))/2
    if(corr) diag(X) <- 1

    structure(list(mat =
		   new("dpoMatrix", x = as.vector(X),
		       Dim = c(n,n), Dimnames = .M.DN(x)),
                   corr = corr, normF = norm(x-X, "F"), iterations = iter,
		   rel.tol = conv, converged = converged),
	      class = "nearPD")
}
