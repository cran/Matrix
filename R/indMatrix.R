#### Index Matrices -- Coercion and Methods (--> ../man/indMatrix-class.Rd )

## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("integer", "indMatrix",
      function(from) {
          if ((m <- length(from)) == 0L)
              return(new("indMatrix"))
          if(anyNA(from))
              stop("'perm' slot cannot contain NA")
          r <- range(from)
          if(r[1L] < 1L)
              stop("elements of 'perm' slot must be positive integers")
          new("indMatrix",
              Dim = c(m, r[2L]),
              Dimnames = list(names(from), NULL),
              perm = from)
      })

setAs("numeric", "indMatrix",
      function(from) {
          if ((m <- length(from)) == 0L)
              return(new("indMatrix"))
          if(anyNA(from))
              stop("'perm' slot cannot contain NA")
          r <- range(from)
          if((r2 <- r[2L]) > .Machine$integer.max)
              stop("elements of 'perm' slot cannot exceed 2^31-1")
          if(r[1L] < 1 || any(from != (from.i <- as.integer(from))))
              stop("elements of 'perm' slot must be positive integers")
          new("indMatrix",
              Dim = c(m, as.integer(r2)),
              Dimnames = list(names(from), NULL),
              perm = from.i)
      })

setAs("list", "indMatrix",
      ## Here, 'from' must be a list of the form 'list(perm, ncol)' ...
      ## needed for the 'max(perm) < ncol' and 'length(perm) == 0L' cases
      function(from) {
          if(length(from) != 2L)
              stop("only lists of length 2 can be coerced to indMatrix")
          n <- n.i <- from[[2L]]
          if(!is.numeric(n) || length(n) != 1L ||
             is.na(n) || n < 0 || n > .Machine$integer.max ||
             (!is.integer(n) && n != (n.i <- as.integer(n))))
              stop("<list>[[2]] must be a nonnegative integer less than or equal to 2^31-1")
          perm <- perm.i <- from[[1L]]
          if(!is.numeric(perm))
              stop("<list>[[1]] must be numeric")
          if(anyNA(perm))
              stop("<list>[[1]] cannot contain NA")
          r <- range(perm)
          if(r[2L] > n.i)
              stop("elements of <list>[[1]] cannot exceed <list>[[2]]")
          if(r[1L] < 1 || (!is.integer(perm) &&
                           any(perm != (perm.i <- as.integer(perm)))))
              stop("elements of <list>[[1]] must be positive integers")
          new("indMatrix",
              Dim = c(length(perm.i), n.i),
              Dimnames = list(names(perm.i), NULL),
              perm = perm.i)
      })

setAs("nsparseMatrix", "indMatrix",
      function(from) {
	  from <- as(as(from, "RsparseMatrix"), "generalMatrix")
          p <- from@p
          m <- length(p) - 1L
          if(m > 0L && any(p != 0:m))
              stop("matrix must have exactly one nonzero element in each row")
          new("indMatrix", Dim = from@Dim, Dimnames = from@Dimnames,
              perm = from@j + 1L)
      })

setAs("Matrix", "indMatrix",
      function(from) as(as(from, "nsparseMatrix"), "indMatrix"))

setAs("matrix", "indMatrix",
      function(from) as(as(from, "nsparseMatrix"), "indMatrix"))


## ~~~~ COERCIONS FROM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.ind2dge <- function(from) {
    x <- double(prod(d <- from@Dim))
    if((n <- d[1L]) > 0L)
        x[seq_len(n) + (from@perm - 1L) * as.double(n)] <- 1
    new("dgeMatrix", Dim = d, Dimnames = from@Dimnames, x = x)
}
.ind2lge <- function(from) {
    x <- logical(prod(d <- from@Dim))
    if((n <- d[1L]) > 0L)
        x[seq_len(n) + (from@perm - 1L) * as.double(n)] <- TRUE
    new("lgeMatrix", Dim = d, Dimnames = from@Dimnames, x = x)
}
.ind2nge <- function(from) {
    x <- logical(prod(d <- from@Dim))
    if((n <- d[1L]) > 0L)
        x[seq_len(n) + (from@perm - 1L) * as.double(n)] <- TRUE
    new("ngeMatrix", Dim = d, Dimnames = from@Dimnames, x = x)
}
.ind2n.p <- function(from) {
    from <-
        if(isSymmetric(from))
            forceSymmetric(from)
        else if(!(it <- isTriangular(from)))
            stop("matrix is not symmetric or triangular")
        else if (attr(it, "kind") == "U")
            triu(from)
        else tril(from)
    .Call(R_sparse_as_dense, from, TRUE)
}
.ind2dgC <- function(from) {
    perm <- from@perm
    d <- from@Dim
    new("dgCMatrix", Dim = d, Dimnames = from@Dimnames,
        p = c(0L, cumsum(tabulate(perm, d[2L]))), i = sort.list(perm) - 1L,
        x = rep.int(1, length(perm)))
}
.ind2lgC <- function(from) {
    perm <- from@perm
    d <- from@Dim
    new("lgCMatrix", Dim = d, Dimnames = from@Dimnames,
        p = c(0L, cumsum(tabulate(perm, d[2L]))), i = sort.list(perm) - 1L,
        x = rep.int(TRUE, length(perm)))
}
.ind2ngC <- function(from) {
    perm <- from@perm
    d <- from@Dim
    new("ngCMatrix", Dim = d, Dimnames = from@Dimnames,
        p = c(0L, cumsum(tabulate(perm, d[2L]))), i = sort.list(perm) - 1L)
}
.ind2ngR <- function(from) {
    perm <- from@perm
    new("ngRMatrix", Dim = from@Dim, Dimnames = from@Dimnames,
        p = 0:length(perm), j = perm - 1L)
}
.ind2ngT <- function(from) {
    perm <- from@perm
    new("ngTMatrix", Dim = from@Dim, Dimnames = from@Dimnames,
        i = if((m <- length(perm)) > 0L) 0:(m-1L) else integer(0L),
        j = perm - 1L)
}
.ind2diag <- function(from) {
    if (!isDiagonal(from))
        stop("matrix is not diagonal; consider Diagonal(x=diag(.))")
    new("ldiMatrix", Dim = from@Dim, Dimnames = from@Dimnames, diag = "U")
}
.ind2p <- function(from) new("pMatrix", from)

setAs("indMatrix",    "denseMatrix", .ind2nge)
setAs("indMatrix", "unpackedMatrix", .ind2nge)
setAs("indMatrix",   "packedMatrix", .ind2n.p)
setAs("indMatrix",         "matrix", .ind2m)
setAs("indMatrix",         "vector", .ind2v)

setAs("indMatrix",        "dMatrix", .ind2dgC)
setAs("indMatrix",  "dsparseMatrix", .ind2dgC)
setAs("indMatrix",   "ddenseMatrix", .ind2dge)
setAs("indMatrix",        "lMatrix", .ind2lgC)
setAs("indMatrix",  "lsparseMatrix", .ind2lgC)
setAs("indMatrix",   "ldenseMatrix", .ind2lge)
setAs("indMatrix",        "nMatrix", .ind2ngC)
setAs("indMatrix",  "nsparseMatrix", .ind2ngC)
setAs("indMatrix",   "ndenseMatrix", .ind2nge)

setAs("indMatrix",  "generalMatrix", .ind2ngC)
## setAs("indMatrix", "triangularMatrix", .) # inherited from Matrix
## setAs("indMatrix",  "symmetricMatrix", .) # inherited from Matrix

setAs("indMatrix",  "CsparseMatrix", .ind2ngC)
setAs("indMatrix",  "RsparseMatrix", .ind2ngR)
setAs("indMatrix",  "TsparseMatrix", .ind2ngT)
setAs("indMatrix", "diagonalMatrix", .ind2diag)
setAs("indMatrix",        "pMatrix", .ind2p)

setMethod("as.vector", signature(x = "indMatrix"),
          function(x, mode) as.vector(.ind2v(x), mode))
setMethod("as.numeric", signature(x = "indMatrix"),
          function(x, ...) as.double(.ind2v(x), mode))
setMethod("as.logical", signature(x = "indMatrix"),
          function(x, ...) .ind2v(x))

## DEPRECATED IN 1.4-2; see ./zzz.R
if(FALSE) {
setAs("indMatrix", "ngTMatrix", .ind2ngT)
setAs("indMatrix", "ngeMatrix", .ind2nge)
} ## DEPRECATED IN 1.4-2; see ./zzz.R

rm(.ind2dge, .ind2lge, .ind2nge, .ind2n.p,
   .ind2dgC, .ind2lgC, .ind2ngC, .ind2ngR, # .ind2ngT,
   .ind2diag, .ind2p)

if(!.Matrix.supporting.cached.methods) {
rm(.ind2ngT)
}


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("isSymmetric", signature(object = "indMatrix"),
	  function(object, checkDN = TRUE, ...) {
	      d <- object@Dim
	      if((n <- d[1L]) != d[2L])
		  return(FALSE)
              if(checkDN) {
                  ca <- function(check.attributes = TRUE, ...) check.attributes
                  if(ca(...) && !isSymmetricDN(object@Dimnames))
                      return(FALSE)
              }
              i <- seq_len(n)
              j <- object@perm
              all(j[j] == i)
	  })

setMethod("isTriangular", signature(object = "indMatrix"),
          function(object, upper = NA, ...) {
              d <- object@Dim
	      if((n <- d[1L]) != d[2L])
		  return(FALSE)
              i <- seq_len(n)
              j <- object@perm
              if(is.na(upper)) {
                  if(all(j >= i))
                      return(`attr<-`(TRUE, "kind", "U"))
                  if(all(i <= j))
                      return(`attr<-`(TRUE, "kind", "L"))
                  FALSE
              } else if(upper) {
                  all(j >= i)
              } else {
                  all(i <= j)
              }
          })

setMethod("isDiagonal", signature(object = "indMatrix"),
          function(object) {
              d <- object@Dim
	      if((n <- d[1L]) != d[2L])
		  return(FALSE)
              all(object@perm == seq_len(n))
          })

setMethod("t", signature(x = "indMatrix"),
          function(x) {
              j <- x@perm
              new("ngTMatrix",
                  Dim = x@Dim[2:1], Dimnames = x@Dimnames[2:1],
                  i = j - 1L,
                  j = if((m <- length(j)) > 0L) 0:(m-1L) else integer(0L))
          })

setMethod("diag", signature(x = "indMatrix"),
          function(x, nrow, ncol, names) {
              if((m <- min(x@Dim)) == 0L)
                  return(logical(0L))
              i <- seq_len(m)
              y <- x@perm[i] == i
              if(names &&
                 !any(vapply(dn <- x@Dimnames, is.null, NA)) &&
                 identical(nms <- dn[[1L]][i], dn[[2L]][i]))
                  names(y) <- nms
              y
          })

setMethod("diag<-", signature(x = "indMatrix"),
          function(x, value) {
              x <- as(x, "nsparseMatrix")
              callGeneric()
          })

setMethod("band", signature(x = "indMatrix"),
          function(x, k1, k2, ...) {
              x <- as(x, "nsparseMatrix")
              callGeneric()
          })

setMethod("triu", signature(x = "indMatrix"),
          function(x, k, ...) {
              x <- as(x, "nsparseMatrix")
              callGeneric()
          })

setMethod("tril", signature(x = "indMatrix"),
          function(x, k, ...) {
              x <- as(x, "nsparseMatrix")
              callGeneric()
          })

setMethod("forceSymmetric", signature(x = "indMatrix", uplo = "missing"),
          function(x, uplo) {
              x <- as(x, "nsparseMatrix")
              callGeneric()
          })

setMethod("forceSymmetric", signature(x = "indMatrix", uplo = "character"),
          function(x, uplo) {
              x <- as(x, "nsparseMatrix")
              callGeneric()
          })

setMethod("%*%", signature(x = "matrix", y = "indMatrix"),
	  function(x, y) x %*% as(y, "lMatrix"))
setMethod("%*%", signature(x = "Matrix", y = "indMatrix"),
	  function(x, y) x %*% as(y, "lMatrix"))

setMethod("%*%", signature(x = "indMatrix", y = "matrix"),
	  function(x, y) { mmultCheck(x,y); y[x@perm ,] })
setMethod("%*%", signature(x = "indMatrix", y = "Matrix"),
	  function(x, y) { mmultCheck(x,y); y[x@perm ,] })


setMethod("crossprod", signature(x = "indMatrix", y = "matrix"),
	  function(x, y) as(t(x), "lMatrix") %*% y)
setMethod("crossprod", signature(x = "indMatrix", y = "Matrix"),
	  function(x, y) as(t(x), "lMatrix") %*% y)
setMethod("crossprod", signature(x = "indMatrix", y = "indMatrix"),
	  function(x, y) {
	      mmultCheck(x,y, 2L)
              ## xy <- interaction(x@perm, y@perm)
              ## this is wrong if any of the columns in X or Y are empty because interaction()
              ## drops non-occuring levels from a non-factor. Explicitly defining a factor with
              ## levels 1:ncol(<indMatrix>) avoids that.
              nx <- x@Dim[2L]
              ny <- y@Dim[2L]
	      ## xy <- interaction(factor(x@perm, levels=seq_len(nx)),
	      ##   		   factor(y@perm, levels=seq_len(ny)))
	      ## much faster (notably for large x,y):
	      xy <- x@perm + nx*as.double(y@perm-1L)
	      Matrix(tabulate(xy, nbins = nx*ny), nrow = nx, ncol = ny,
		     dimnames = list(x@Dimnames[[2L]], y@Dimnames[[2L]]))
	  })

setMethod("tcrossprod", signature(x = "matrix", y = "indMatrix"),
	  function(x, y) { mmultCheck(x,y, 3L); x[, y@perm] })
setMethod("tcrossprod", signature(x = "Matrix", y = "indMatrix"),
	  function(x, y) { mmultCheck(x,y, 3L); x[, y@perm] })
setMethod("tcrossprod", signature(x = "indMatrix", y = "indMatrix"),
	  function(x, y) { mmultCheck(x,y, 3L); x[, y@perm] })

setMethod("crossprod", signature(x = "indMatrix", y = "missing"),
	  function(x, y=NULL) Diagonal(x = tabulate(x@perm, nbins=x@Dim[2L])))

setMethod("tcrossprod", signature(x = "indMatrix", y = "missing"),
	  function(x, y=NULL) x[,x@perm])


setMethod("kronecker", signature(X = "indMatrix", Y = "indMatrix"),
	  function (X, Y, FUN = "*", make.dimnames = FALSE, ...) {
	      if (FUN != "*") stop("kronecker method must use default 'FUN'")
	      if(any(as.double(X@Dim)*Y@Dim >= .Machine$integer.max))
		  stop("resulting matrix dimension would be too large")
	      ## Explicitly defining a factor with levels 1:ncol(.) avoids that
	      ## interaction() drops non-occuring levels when any of the
	      ## columns in X or Y are empty:
	      ## perm <-  as.integer(interaction(factor(rep(X@perm, each =Y@Dim[1]),
	      ##                                        levels=seq_len(X@Dim[2])),
	      ##                                 factor(rep.int(Y@perm, times=X@Dim[1]),
	      ##                                        levels=seq_len(Y@Dim[2])),
	      ##                                 lex.order=TRUE))
	      ## much faster (notably for large X, Y):
	      fX <- rep    (X@perm-1L, each  = Y@Dim[1])
	      fY <- rep.int(Y@perm-1L, times = X@Dim[1])
	      new("indMatrix", perm = 1L + fY + Y@Dim[2] * fX,
		  Dim = X@Dim*Y@Dim)
	  })


setMethod("[", signature(x = "indMatrix", i = "index", j = "missing",
			 drop = "logical"),
	  function (x, i, j, ..., drop)
      {
	  n <- length(newperm <- x@perm[i])
	  if(drop && n == 1) { ## -> logical unit vector
	      newperm == seq_len(x@Dim[2])
	  } else { ## stay matrix
	      if(!is.null((DN <- x@Dimnames)[[1]])) DN[[1]] <- DN[[1]][i]
	      new("indMatrix", perm = newperm,
		  Dim = c(n, x@Dim[2]), Dimnames = DN)
	  }
      })


.indMatrix.sub <- function(x, i, j, ..., value) {
    x <- as(x, "TsparseMatrix")
    callGeneric()
}
for (.i in c("missing", "index"))
    for (.j in c("missing", "index"))
        setReplaceMethod("[", signature(x = "indMatrix", i = .i, j = .j),
                         .indMatrix.sub)
rm(.indMatrix.sub, .i, .j)
