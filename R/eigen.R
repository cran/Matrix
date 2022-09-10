#### eigen() , Schur() etc
#### =====     =====

## eigen() is not even generic, and we haven't any C code,
## NOTE  base::eigen()  "magically"  can work via as.matrix()
if(.Matrix.avoiding.as.matrix) {
    ## ---- IFF  as.matrix(.)  <==>  as(., "matrix")  [which we consider _deprecating_]
    ## FIXME: Code for *sparse* !! [RcppEigen ~??~]
    setMethod("eigen", signature(x = "Matrix", only.values = "missing"),
	      function(x, symmetric, only.values, EISPACK) # << must match generic
		  base::eigen(as(x,"matrix"), symmetric, FALSE))
    setMethod("eigen", signature(x = "Matrix", only.values = "logical"),
	      function(x, symmetric, only.values, EISPACK)
		  base::eigen(as(x,"matrix"), symmetric, only.values))

    ## base::svd()  using  as.matrix() :=  asRbasematrix()
    setMethod("svd", "Matrix",
	      function(x, ...) base::svd(as(x, "matrix"), ...))
}


## FIXME: C-level code should scan first for non-finite values
##        since R-level test needs an allocation

setMethod("Schur", signature(x = "Matrix", vectors = "missing"),
          function(x, vectors, ...) Schur(x, TRUE, ...))

setMethod("Schur", signature(x = "matrix", vectors = "missing"),
          function(x, vectors, ...) Schur(x, TRUE, ...))

setMethod("Schur", signature(x = "dgeMatrix", vectors = "logical"),
	  function(x, vectors, ...) {
              if(!all(is.finite(x@x)))
                  stop("'x' has non-finite values")
              cl <- .Call(dgeMatrix_Schur, x, vectors, TRUE)
              if(all(cl$WI == 0)) {
                  vals <- cl$WR
                  T <- triu(cl$T)
              } else {
                  vals <- complex(real = cl$WR, imaginary = cl$WI)
                  T <- .m2ge(cl$T)
              }
              if(vectors)
                  new("Schur", Dim = x@Dim, Q = .m2ge(cl$Z), T = T,
                      EValues = vals)
              else list(T = T, EValues = vals)
          })

setMethod("Schur", signature(x = "dsyMatrix", vectors = "logical"),
	  function(x, vectors, ...) {
              e <- eigen(x, symmetric = TRUE, only.values = !vectors)
              vals <- as.double(e$values)
              T <- new("ddiMatrix", Dim = x@Dim, x = vals)
              if(vectors)
                  new("Schur", Dim = x@Dim, Q = .m2ge(e$vectors), T = T,
                      EValues = vals)
              else list(T = T, EValues = vals)
          })

setMethod("Schur", signature(x = "matrix", vectors = "logical"),
	  function(x, vectors, ...) {
              storage.mode(x) <- "double"
              if(!all(is.finite(x)))
                  stop("'x' has non-finite values")
              cl <- .Call(dgeMatrix_Schur, x, vectors, FALSE)
              vals <-
                  if(all(cl$WI == 0))
                      cl$WR
                  else complex(real = cl$WR, imaginary = cl$WI)
              if(vectors)
                  list(Q = cl$Z, T = cl$T, EValues = vals)
              else list(T = cl$T, EValues = vals)
          })

## FIXME: don't coerce from sparse to dense
setMethod("Schur", signature(x = "generalMatrix", vectors = "logical"),
	  function(x, vectors, ...)
              Schur(as(as(x, "dMatrix"), "unpackedMatrix"), vectors, ...))

## FIXME: don't coerce from sparse to dense
setMethod("Schur", signature(x = "symmetricMatrix", vectors = "logical"),
	  function(x, vectors, ...)
              Schur(as(as(x, "dMatrix"), "unpackedMatrix"), vectors, ...))

## Giving the _unsorted_ Schur factorization
setMethod("Schur", signature(x = "diagonalMatrix", vectors = "logical"),
	  function(x, vectors, ...) {
              d <- x@Dim
              if(x@diag != "N") {
                  vals <- rep.int(1, d[1L])
                  T <- new("ddiMatrix", Dim = d, diag = "U")
              } else {
                  vals <- x@x
                  if(!all(is.finite(vals)))
                      stop("'x' has non-finite values")
                  T <- new("ddiMatrix", Dim = d, x = vals)
              }
              if(vectors) {
                  Q <- new("ddiMatrix", Dim = d, diag = "U")
                  new("Schur", Dim = d, Q = Q, T = T, EValues = vals)
              } else list(T = T, EValues = vals)
          })

setMethod("Schur", signature(x = "triangularMatrix", vectors = "logical"),
	  function(x, vectors, ...) {
              cld <- getClassDef(class(x))
              if(!extends(cld, "nMatrix") &&
                 (anyNA(x) || (extends(cld, "dMatrix") && any(is.infinite(x)))))
                  ## any(is.finite(<sparseMatrix>)) would allocate too much
                  stop("'x' has non-finite values")
              n <- (d <- x@Dim)[1L]
              vals <- diag(x, names = FALSE)
              if(x@uplo == "U" || n == 0L) {
                  if(vectors) {
                      Q <- new("ddiMatrix", Dim = d, diag = "U")
                      new("Schur", Dim = d, Q = Q, T = x, EValues = vals)
                  } else list(T = x, EValues = vals)
              } else {
                  perm <- n:1L
                  vals <- vals[perm]
                  T <- triu(x[perm, perm, drop = FALSE])
                  if(vectors) {
                      Q <- new("pMatrix", Dim = d, perm = perm)
                      new("Schur", Dim = d, Q = Q, T = T, EValues = vals)
                  } else list(T = x, EValues = vals)
              }
          })
