## METHODS FOR GENERIC: solve
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.solve.checkDim1 <-
function(nrow.a, ncol.a) {
    if(nrow.a != ncol.a)
        stop("'a' is not square")
}

.solve.checkDim2 <-
function(nrow.a, nrow.b) {
    if(nrow.a != nrow.b)
        stop("dimensions of 'a' and 'b' are inconsistent")
}

.solve.checkCond <-
function(a, tol, rcond.a = rcond(a)) {
    if(tol > 0 && a@Dim[1L] > 0L && rcond.a < tol)
        stop(gettextf("'a' is computationally singular, rcond(a)=%g", rcond.a),
             domain = NA)
}

## Factorize  A  as  P1 A P2 = L U .  Then :
##
##     kappa(A) <= kappa(P1') * kappa(L) * kappa(U) * kappa(P2')
##              ==          1 * kappa(L) * kappa(U) * 1
##
## If  U  is diagonally dominant, i.e., if there exists  a > 1  such that :
##
##     abs(U[i,i]) >= a * sum(abs(U[i,-i]))        i = 1, ... , n
##
## then :
##
##     kappa(U) <= ((a + 1) / (a - 1)) * max(abs(diag(U))) / min(abs(diag(U)))
##
## The bound contracts as  a --> Inf  and in the limit we have:
##
##     kappa(A) / kappa(L) <= kappa(U) <= max(abs(diag(U))) / min(abs(diag(U)))
##
.solve.checkCondBound <-
function(u, tol, rad.u = range(abs(diag(u, names = FALSE)))) {
    if(tol > 0 && u@Dim[1L] > 0L) {
        r <- rad.u[1L] / rad.u[2L]
        if(r < tol)
            stop(gettextf("'a' is computationally singular, min(d)/max(d)=%g, d=abs(diag(U))",
                          r),
                 domain = NA)
    }
}

.solve.checkDiagonal <-
function(diag.a) {
    zero <- diag.a == 0
    if(any(zero, na.rm = TRUE))
        stop(gettextf("matrix is exactly singular, D[i,i]=0, i=%d",
                      which.max(zero)),
             domain = NA)
}

.solve.checkInd <-
function(perm, margin, n) {
    if(anyDuplicated.default(perm)) {
        i <- seq_len(n)[-perm][1L]
        if(margin == 1L)
            stop(gettextf("matrix is exactly singular, J[,j]=0, j=%d", i),
                 domain = NA)
        else
            stop(gettextf("matrix exactly singular, J[i,]=0, i=%d", i),
                 domain = NA)
    }
}


########################################################################
##  1. MatrixFactorization incl. triangularMatrix
########################################################################

for(.cl in c("MatrixFactorization", "triangularMatrix")) {

setMethod("solve", signature(a = .cl, b = "numLike"),
          function(a, b, ...)
              drop(solve(a, .m2ge(b, "d"), ...)))

setMethod("solve", signature(a = .cl, b = "matrix"),
          function(a, b, ...)
              solve(a, .m2ge(b, "d"), ...))

setMethod("solve", signature(a = .cl, b = "denseMatrix"),
          function(a, b, ...)
              solve(a, .dense2g(b, "d"), ...))

setMethod("solve", signature(a = .cl, b = "CsparseMatrix"),
          function(a, b, ...)
              solve(a, .sparse2g(b, "d"), ...))

setMethod("solve", signature(a = .cl, b = "RsparseMatrix"),
          function(a, b, ...)
              solve(a, .sparse2g(.CR2RC(b), "d"), ...))

setMethod("solve", signature(a = .cl, b = "TsparseMatrix"),
          function(a, b, ...)
              solve(a, .sparse2g(.T2C(b), "d"), ...))

setMethod("solve", signature(a = .cl, b = "diagonalMatrix"),
          function(a, b, ...)
              solve(a, .diag2sparse(b, "dgC"), ...))

setMethod("solve", signature(a = .cl, b = "indMatrix"),
          function(a, b, ...)
              solve(a, as(b, "dMatrix"), ...))

setMethod("solve", signature(a = .cl, b = "dgeMatrix"),
          function(a, b, ...)
              solve(a, .dense2sparse(b, "C"), ...))

setMethod("solve", signature(a = .cl, b = "dgCMatrix"),
          function(a, b, ...)
              solve(a, .sparse2dense(b, FALSE), ...))

}

setMethod("solve", signature(a = "denseLU", b = "missing"),
          function(a, b, ...)
              .Call(denseLU_solve, a, NULL))

setMethod("solve", signature(a = "denseLU", b = "dgeMatrix"),
          function(a, b, ...)
              .Call(denseLU_solve, a, b))

setMethod("solve", signature(a = "BunchKaufman", b = "missing"),
          function(a, b, ...)
              .Call(BunchKaufman_solve, a, NULL, FALSE))

setMethod("solve", signature(a = "BunchKaufman", b = "dgeMatrix"),
          function(a, b, ...)
              .Call(BunchKaufman_solve, a, b, FALSE))

setMethod("solve", signature(a = "pBunchKaufman", b = "missing"),
          function(a, b, ...)
              .Call(BunchKaufman_solve, a, NULL, TRUE))

setMethod("solve", signature(a = "pBunchKaufman", b = "dgeMatrix"),
          function(a, b, ...)
              .Call(BunchKaufman_solve, a, b, TRUE))

setMethod("solve", signature(a = "Cholesky", b = "missing"),
          function(a, b, ...)
              .Call(Cholesky_solve, a, NULL, FALSE))

setMethod("solve", signature(a = "Cholesky", b = "dgeMatrix"),
          function(a, b, ...)
              .Call(Cholesky_solve, a, b, FALSE))

setMethod("solve", signature(a = "pCholesky", b = "missing"),
          function(a, b, ...)
              .Call(Cholesky_solve, a, NULL, TRUE))

setMethod("solve", signature(a = "pCholesky", b = "dgeMatrix"),
          function(a, b, ...)
              .Call(Cholesky_solve, a, b, TRUE))

setMethod("solve", signature(a = "pCholesky", b = "missing"),
          function(a, b, ...)
              .Call(Cholesky_solve, a, NULL, TRUE))

setMethod("solve", signature(a = "pCholesky", b = "dgeMatrix"),
          function(a, b, ...)
              .Call(Cholesky_solve, a, b, TRUE))

setMethod("solve", signature(a = "sparseLU", b = "missing"),
          function(a, b, tol = .Machine$double.eps, sparse = TRUE, ...) {
              .solve.checkCondBound(a@U, tol)
              .Call(sparseLU_solve, a, NULL, sparse)
          })

setMethod("solve", signature(a = "sparseLU", b = "dgeMatrix"),
          function(a, b, tol = .Machine$double.eps, ...) {
              .solve.checkCondBound(a@U, tol)
              .Call(sparseLU_solve, a, b, FALSE)
          })

setMethod("solve", signature(a = "sparseLU", b = "dgCMatrix"),
          function(a, b, tol = .Machine$double.eps, ...) {
              .solve.checkCondBound(a@U, tol)
              .Call(sparseLU_solve, a, b, TRUE)
          })

setMethod("solve", signature(a = "sparseQR", b = "missing"),
          function(a, b, sparse = TRUE, ...) {
              r <- qr.coef(a, diag(a@Dim[2L]))
              if(is.na(sparse) || sparse) .dense2sparse(r, "C") else r
          })

setMethod("solve", signature(a = "sparseQR", b = "dgeMatrix"),
          function(a, b, ...)
              qr.coef(a, b))

setMethod("solve", signature(a = "sparseQR", b = "dgCMatrix"),
          function(a, b, ...)
              .dense2sparse(qr.coef(a, .sparse2dense(b, FALSE)), "C"))

setMethod("solve", signature(a = "CHMfactor", b = "missing"),
          function(a, b, sparse = TRUE,
                   system = c("A","LDLt","LD","DLt","L","Lt","D","P","Pt"), ...) {
              if((is.na(sparse) || sparse) && !missing(system)) {
                  ## Do the "right" thing :
                  if(identical(system, "D")) {
                      r <- new("ddiMatrix")
                      r@Dim <- a@Dim
                      r@Dimnames <- a@Dimnames[2:1]
                      if(.CHM.is.LDL(a))
                          r@x <- a@x[a@p + 1L]
                      else r@diag <- "U"
                      return(r)
                  } else if(identical(system, "P") || identical(system, "Pt")) {
                      r <- new("pMatrix")
                      r@Dim <- a@Dim
                      r@Dimnames <- a@Dimnames[2:1]
                      if(system == "Pt")
                          r@margin <- 2L
                      r@perm <- a@perm + 1L
                      return(r)
                  }
              }
              .Call(CHMfactor_solve, a, NULL, sparse, system)
          })

setMethod("solve", signature(a = "CHMfactor", b = "dgeMatrix"),
          function(a, b,
                   system = c("A","LDLt","LD","DLt","L","Lt","D","P","Pt"), ...)
              .Call(CHMfactor_solve, a, b, FALSE, system))

setMethod("solve", signature(a = "CHMfactor", b = "dgCMatrix"),
          function(a, b,
                   system = c("A","LDLt","LD","DLt","L","Lt","D","P","Pt"), ...)
              .Call(CHMfactor_solve, a, b, TRUE, system))

setMethod("solve", signature(a = "dtrMatrix", b = "missing"),
          function(a, b, tol = .Machine$double.eps, ...) {
              .solve.checkCond(a, tol)
              .Call(dtrMatrix_solve, a, NULL, FALSE)
          })

setMethod("solve", signature(a = "dtrMatrix", b = "dgeMatrix"),
          function(a, b, tol = .Machine$double.eps, ...) {
              .solve.checkCond(a, tol)
              .Call(dtrMatrix_solve, a, b, FALSE)
          })

setMethod("solve", signature(a = "dtpMatrix", b = "missing"),
          function(a, b, tol = .Machine$double.eps, ...) {
              .solve.checkCond(a, tol)
              .Call(dtrMatrix_solve, a, NULL, TRUE)
          })

setMethod("solve", signature(a = "dtpMatrix", b = "dgeMatrix"),
          function(a, b, tol = .Machine$double.eps, ...) {
              .solve.checkCond(a, tol)
              .Call(dtrMatrix_solve, a, b, TRUE)
          })

setMethod("solve", signature(a = "dtCMatrix", b = "missing"),
          function(a, b, sparse = TRUE, ...) {
              if(a@diag != "N")
                  a <- ..diagU2N(a)
              .Call(dtCMatrix_solve, a, NULL, sparse)
          })

setMethod("solve", signature(a = "dtCMatrix", b = "dgeMatrix"),
          function(a, b, sparse = FALSE, ...) {
              if(a@diag != "N")
                  a <- ..diagU2N(a)
              if(is.na(sparse) || sparse)
                  b <- .dense2sparse(b, "C")
              .Call(dtCMatrix_solve, a, b, sparse)
          })

setMethod("solve", signature(a = "dtCMatrix", b = "dgCMatrix"),
          function(a, b, sparse = TRUE, ...) {
              if(a@diag != "N")
                  a <- ..diagU2N(a)
              if(!(is.na(sparse) || sparse))
                  b <- .sparse2dense(b, FALSE)
              .Call(dtCMatrix_solve, a, b, sparse)
          })

for(.cl in c("dtrMatrix", "dtpMatrix", "dtCMatrix"))
setMethod("solve", signature(a = .cl, b = "triangularMatrix"),
          function(a, b, ...) {
              r <- solve(a, as(b, "generalMatrix"), ...)
              if(a@uplo == b@uplo) {
                  r <- if(a@uplo == "U") triu(r) else tril(r)
                  if(a@diag != "N" && b@diag != "N")
                      r <- ..diagN2U(r, sparse = .isCRT(r))
              }
              r
          })
rm(.cl)

## MJ: truly an exceptional case ...
setMethod("solve", signature(a = "Schur", b = "ANY"),
          function(a, b, ...) {
              Q <- a@Q
              T <- a@T
              if(missing(b)) {
                  r <- Q %*% solve(T, t(Q))
                  r@Dimnames <- a@Dimnames[2:1]
                  r
              } else {
                  db <- dim(b)
                  dnb <- dimnames(b)
                  r <- Q %*% solve(T, crossprod(Q, b))
                  r@Dimnames <- c(a@Dimnames[2L],
                                  if(is.null(dnb)) list(NULL) else dnb[2L])
                  if(is.null(db)) drop(r) else r
              }
          })


########################################################################
##  2. denseMatrix excl. triangularMatrix
########################################################################

setMethod("solve", signature(a = "denseMatrix", b = "ANY"),
          function(a, b, ...) {
              a <- ..dense2d(a)
              if(missing(b)) solve(a, ...) else solve(a, b, ...)
          })

setMethod("solve", signature(a = "dgeMatrix", b = "ANY"),
          function(a, b, tol = .Machine$double.eps, ...) {
              d <- a@Dim
              .solve.checkDim1(d[1L], d[2L])
              .solve.checkCond(a, tol)
              trf <- lu(a, warnSing = FALSE)
              if(missing(b)) solve(trf, ...) else solve(trf, b, ...)
          })

for(.cl in c("dsyMatrix", "dspMatrix"))
setMethod("solve", signature(a = .cl, b = "ANY"),
          function(a, b, tol = .Machine$double.eps, ...) {
              .solve.checkCond(a, tol)
              trf <- BunchKaufman(a, warnSing = FALSE)
              if(missing(b)) solve(trf, ...) else solve(trf, b, ...)
          })
rm(.cl)

for(.cl in c("dpoMatrix", "dppMatrix"))
setMethod("solve", signature(a = .cl, b = "ANY"),
          function(a, b, tol = .Machine$double.eps, ...) {
              .solve.checkCond(a, tol)
              trf <- Cholesky(a, perm = FALSE)
              if(missing(b)) solve(trf, ...) else solve(trf, b, ...)
          })
rm(.cl)


########################################################################
##  3. CsparseMatrix excl. triangularMatrix
########################################################################

setMethod("solve", signature(a = "CsparseMatrix", b = "ANY"),
          function(a, b, ...) {
              a <- ..sparse2d(a)
              if(missing(b)) solve(a, ...) else solve(a, b, ...)
          })

setMethod("solve", signature(a = "dgCMatrix", b = "missing"),
          function(a, b, sparse = TRUE, ...) {
              trf <- lu(a, errSing = TRUE)
              solve(trf, sparse = sparse, ...)
          })

setMethod("solve", signature(a = "dgCMatrix", b = "numLike"),
          function(a, b, ...) {
              trf <- lu(a, errSing = TRUE)
              solve(trf, b, ...)
          })

setMethod("solve", signature(a = "dgCMatrix", b = "matrix"),
          function(a, b, sparse = FALSE, ...) {
              trf <- lu(a, errSing = TRUE)
              if(is.na(sparse) || sparse)
                  b <- .m2sparse(b, "dgC")
              solve(trf, b, ...)
          })

setMethod("solve", signature(a = "dgCMatrix", b = "denseMatrix"),
          function(a, b, sparse = FALSE, ...) {
              trf <- lu(a, errSing = TRUE)
              if(is.na(sparse) || sparse)
                  b <- as(b, "CsparseMatrix")
              solve(trf, b, ...)
          })

setMethod("solve", signature(a = "dgCMatrix", b = "sparseMatrix"),
          function(a, b, sparse = TRUE, ...) {
              trf <- lu(a, errSing = TRUE)
              if(!(is.na(sparse) || sparse))
                  b <- as(b, "unpackedMatrix")
              solve(trf, b, ...)
          })

setMethod("solve", signature(a = "dsCMatrix", b = "missing"),
          function(a, b, sparse = TRUE, ...) {
              trf <- tryCatch(
                  Cholesky(a, perm = TRUE, LDL = TRUE, super = FALSE),
                  error = function(e) lu(a, errSing = TRUE))
              solve(trf, sparse = sparse, ...)
          })

setMethod("solve", signature(a = "dsCMatrix", b = "numLike"),
          function(a, b, ...) {
              trf <- tryCatch(
                  Cholesky(a, perm = TRUE, LDL = TRUE, super = FALSE),
                  error = function(e) lu(a, errSing = TRUE))
              solve(trf, b, ...)
          })

setMethod("solve", signature(a = "dsCMatrix", b = "matrix"),
          function(a, b, sparse = FALSE, ...) {
              trf <- tryCatch(
                  Cholesky(a, perm = TRUE, LDL = TRUE, super = FALSE),
                  error = function(e) lu(a, errSing = TRUE))
              if(is.na(sparse) || sparse)
                  b <- .m2sparse(b, "dgC")
              solve(trf, b, ...)
          })

setMethod("solve", signature(a = "dsCMatrix", b = "denseMatrix"),
          function(a, b, sparse = FALSE, ...) {
              trf <- tryCatch(
                  Cholesky(a, perm = TRUE, LDL = TRUE, super = FALSE),
                  error = function(e) lu(a, errSing = TRUE))
              if(is.na(sparse) || sparse)
                  b <- as(b, "CsparseMatrix")
              solve(trf, b, ...)
          })

setMethod("solve", signature(a = "dsCMatrix", b = "sparseMatrix"),
          function(a, b, sparse = TRUE, ...) {
              trf <- tryCatch(
                  Cholesky(a, perm = TRUE, LDL = TRUE, super = FALSE),
                  error = function(e) lu(a, errSing = TRUE))
              if(!(is.na(sparse) || sparse))
                  b <- as(b, "unpackedMatrix")
              solve(trf, b, ...)
          })


########################################################################
##  4. RsparseMatrix excl. triangularMatrix
########################################################################

## TODO: implement triangular solver for dtRMatrix, so that we can handle
##       A = <dgRMatrix>  and  A' = .tCR2RC(A)  like so:
##
##                   P1 A' P2 = L U
##       A x = b  <==================>  x = P1' inv(L') inv(U') P2' b
##

setMethod("solve", signature(a = "RsparseMatrix", b = "ANY"),
          function(a, b, ...) {
              a <- ..sparse2d(.CR2RC(a))
              if(missing(b)) solve(a, ...) else solve(a, b, ...)
          })


########################################################################
##  5. TsparseMatrix excl. triangularMatrix
########################################################################

setMethod("solve", signature(a = "TsparseMatrix", b = "ANY"),
          function(a, b, ...) {
              a <- ..sparse2d(.T2C(a))
              if(missing(b)) solve(a, ...) else solve(a, b, ...)
          })


########################################################################
##  6. diagonalMatrix
########################################################################

setMethod("solve", signature(a = "diagonalMatrix", b = "ANY"),
          function(a, b, ...) {
              a <- ..diag2d(a)
              if(missing(b)) solve(a, ...) else solve(a, b, ...)
          })

setMethod("solve", signature(a = "ddiMatrix", b = "missing"),
          function(a, b, ...) {
              if(a@diag == "N") {
                  x <- a@x
                  .solve.checkDiagonal(x)
                  a@x <- 1 / x
              }
              a@Dimnames <- a@Dimnames[2:1]
              a
          })

setMethod("solve", signature(a = "ddiMatrix", b = "numLike"),
          function(a, b, ...) {
              m <- length(b)
              .solve.checkDim2(a@Dim[1L], m)
              r <-
                  if(a@diag == "N") {
                      x <- a@x
                      .solve.checkDiagonal(x)
                      as.double(b) / x
                  } else as.double(b)
              names(r) <- a@Dimnames[[2L]]
              r
          })

setMethod("solve", signature(a = "ddiMatrix", b = "matrix"),
          function(a, b, ...) {
              d <- dim(b)
              .solve.checkDim2(a@Dim[1L], d[1L])
              dn <- dimnames(b)
              r <- new("dgeMatrix")
              r@Dim <- d
              r@Dimnames <- c(a@Dimnames[2L],
                              if(is.null(dn)) list(NULL) else dn[2L])
              r@x <-
              if(a@diag == "N") {
                  x <- a@x
                  .solve.checkDiagonal(x)
                  as.double(b) / x
              } else as.double(b)
              r
          })

setMethod("solve", signature(a = "ddiMatrix", b = "Matrix"),
          function(a, b, ...) {
              .solve.checkDim2(a@Dim[1L], b@Dim[1L])
              if(a@diag == "N") {
                  x <- a@x
                  .solve.checkDiagonal(x)
                  a@x <- 1 / x
              }
              a@Dimnames <- a@Dimnames[2:1]
              a %*% b
          })


########################################################################
##  7. indMatrix
########################################################################

setMethod("solve", signature(a = "indMatrix", b = "ANY"),
          function(a, b, ...) {
              d <- a@Dim
              .solve.checkDim1(d[1L], d[2L])
              perm <- a@perm
              margin <- a@margin
              .solve.checkInd(perm, margin, d[1L])
              p <- new("pMatrix")
              p@Dim <- d
              p@Dimnames <- a@Dimnames
              p@perm <- perm
              p@margin <- margin
              if(missing(b)) solve(p, ...) else solve(p, b, ...)
          })

setMethod("solve", signature(a = "pMatrix", b = "missing"),
          function(a, b, ...) {
              a@Dimnames <- a@Dimnames[2:1]
              a@margin <- if(a@margin == 1L) 2L else 1L
              a
          })

setMethod("solve", signature(a = "pMatrix", b = "numLike"),
          function(a, b, ...) {
              m <- length(b)
              .solve.checkDim2(a@Dim[1L], m)
              perm <- if(a@margin == 1L) invertPerm(a@perm) else a@perm
              r <- as.double(b)[perm]
              names(r) <- a@Dimnames[[2L]]
              r
          })

setMethod("solve", signature(a = "pMatrix", b = "matrix"),
          function(a, b, ...) {
              d <- dim(b)
              .solve.checkDim2(a@Dim[1L], d[1L])
              dn <- dimnames(b)
              perm <- if(a@margin == 1L) invertPerm(a@perm) else a@perm
              r <- new("dgeMatrix")
              r@Dim <- d
              r@Dimnames <- c(a@Dimnames[2L],
                              if(is.null(dn)) list(NULL) else dn[2L])
              r@x <- as.double(b[perm, , drop = FALSE])
              r
          })

setMethod("solve", signature(a = "pMatrix", b = "Matrix"),
          function(a, b, ...) {
              .solve.checkDim2(a@Dim[1L], b@Dim[1L])
              perm <- if(a@margin == 1L) invertPerm(a@perm) else a@perm
              r <- b[perm, , drop = FALSE]
              r@Dimnames <- c(a@Dimnames[2L], b@Dimnames[2L])
              r
          })


########################################################################
##  8. Other ... a=matrix  OR  b=sparseVector
########################################################################

## for now ... fast for this special case ...
.spV2dgC <- function(x) {
    if(is.double(m <- x@length)) {
        if(m >= .Machine$integer.max + 1)
            stop("dimensions cannot exceed 2^31-1")
        m <- as.integer(m)
    }
    i <- as.integer(x@i) - 1L
    nnz <- length(i)
    r <- new("dgCMatrix")
    r@Dim <- c(m, 1L)
    r@p <- c(0L, nnz)
    r@i <- i
    r@x <-
        if(!.hasSlot(x, "x"))
            rep.int(1, nnz)
        else if(is.complex(y <- x@x))
            stop("cannot coerce zsparseVector to dgCMatrix")
        else y
    r
}

## for now ... fast for this special case ...
.spV2dge <- function(x) {
    m <- x@length
    if(is.double(m)) {
        if(m >= .Machine$integer.max + 1)
            stop("dimensions cannot exceed 2^31-1")
        m <- as.integer(m)
    }
    r <- new("dgeMatrix")
    r@Dim <- c(m, 1L)
    r@x <- replace(double(m), x@i,
        if(!.hasSlot(x, "x"))
            1
        else if(is.complex(y <- x@x))
            stop("cannot coerce zsparseVector to dgeMatrix")
        else y)
    r
}

setMethod("solve", signature(a = "Matrix", b = "sparseVector"),
          function(a, b, ...)
              solve(a, .spV2dgC(b), ...)) # FIXME? drop(.)?

setMethod("solve", signature(a = "MatrixFactorization", b = "sparseVector"),
          function(a, b, ...)
              solve(a, .spV2dgC(b), ...)) # FIXME? drop(.)?

setMethod("solve", signature(a = "matrix", b = "Matrix"),
          function(a, b, ...)
              solve(.m2ge(a, "d"), b, ...))

setMethod("solve", signature(a = "matrix", b = "sparseVector"),
          function(a, b, ...)
              solve(.m2ge(a, "d"), .spV2dge(b), ...)) # FIXME? drop(.)?


########################################################################
##  9. Exported solvers
########################################################################

## a=dgCMatrix
## b=vector, matrix, or Matrix
## x=dg[Ce]Matrix
.solve.dgC.lu <- function(a, b, tol = .Machine$double.eps, check = TRUE) {
    if(check && !is(a, "dgCMatrix"))
        a <- as(as(as(a, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    trf <- lu(a, errSing = TRUE)
    solve(trf, b, tol = tol)
}

## a=dgCMatrix
## b=vector or 1-column matrix
## x=list(L, coef, Xty, resid)
.solve.dgC.chol <- function(a, b, check = TRUE) { # -> MatrixModels
    if(check && !is(a, "dgCMatrix"))
        a <- as(as(as(a, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    .Call(dgCMatrix_cholsol, a, b)
}

## *The* interface to cs_qrsol()
## a=dgCMatrix
## b=vector or 1-column matrix  {FIXME in ../src/dgCMatrix.c}
## x=double vector
.solve.dgC.qr <- function(a, b, order = 3L, check = TRUE) { # -> MatrixModels
    if(check && !is(a, "dgCMatrix"))
        a <- as(as(as(a, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    .Call(dgCMatrix_qrsol, a, b, order)
}
