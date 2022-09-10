## METHODS FOR GENERIC: solve
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## FIXME:
## ? solve(<Matrix>, <vector>) should be a dimensionless vector
## ? solve(<denseMatrix>, <sparseVector>) should be a dimensionless vector
## ? solve(<sparseMatrix>, <sparseVector>) should be a _sparseVector_
## * methods should behave consistently in the 0-by-0 case
## * dimension checking should happen before any potentially large allocations
## * many of these do not handle '[dD]imnames' correctly;
##   we need something like:
##   - dimnames(solve(a, b)) = c(dimnames(a)[2L], dimnames(b)[2L])
##   - dimnames(solve(a   )) = dimnames(a)[2:1]
##   -    names(solve(a, b)) = dimnames(a)[[2L]]    {{ for vector 'b' }}
##   names(dimnames(base::solve.default(a, b))) is NULL always; bug??


## ~~~~ denseMatrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## TODO: we could avoid the proliferation of methods seen here
##       by dispatching at C-level with R_check_class_etc() and
##       having at R-level one method for all ddenseMatrix 'a' ...

## Dispatch to methods for ddenseMatrix (below)
setMethod("solve", signature(a = "denseMatrix", b = "ANY"),
	  function(a, b, ...) {a <- ..dense2d(a); callGeneric()})

setMethod("solve", signature(a = "dgeMatrix", b = "missing"),
	  function(a, b, ...) .Call(dgeMatrix_solve, a))
setMethod("solve", signature(a = "dtrMatrix", b = "missing"),
	  function(a, b, ...) .Call(dtrMatrix_solve, a))
setMethod("solve", signature(a = "dtpMatrix", b = "missing"),
	  function(a, b, ...) .Call(dtpMatrix_solve, a))
setMethod("solve", signature(a = "dsyMatrix", b = "missing"),
          function(a, b, ...) .Call(dsyMatrix_solve, a))
setMethod("solve", signature(a = "dspMatrix", b = "missing"),
	  function(a, b, ...) .Call(dspMatrix_solve, a))
setMethod("solve", signature(a = "dpoMatrix", b = "missing"),
          function(a, b, ...) .Call(dpoMatrix_solve, a))
setMethod("solve", signature(a = "dppMatrix", b = "missing"),
          function(a, b, ...) .Call(dppMatrix_solve, a))

for(.cl in c("numLike", "matrix")) {
setMethod("solve", signature(a = "dgeMatrix", b = .cl),
	  function(a, b, ...) .Call(dgeMatrix_matrix_solve, a, b))
setMethod("solve", signature(a = "dtrMatrix", b = .cl),
	  function(a, b, ...) .Call(dtrMatrix_matrix_solve, a, b))
setMethod("solve", signature(a = "dtpMatrix", b = .cl),
	  function(a, b, ...) .Call(dtpMatrix_matrix_solve, a, b))
setMethod("solve", signature(a = "dsyMatrix", b = .cl),
          function(a, b, ...) .Call(dsyMatrix_matrix_solve, a, b))
setMethod("solve", signature(a = "dspMatrix", b = .cl),
	  function(a, b, ...) .Call(dspMatrix_matrix_solve, a, b))
setMethod("solve", signature(a = "dpoMatrix", b = .cl),
          function(a, b, ...) .Call(dpoMatrix_matrix_solve, a, b))
setMethod("solve", signature(a = "dppMatrix", b = .cl),
          function(a, b, ...) .Call(dppMatrix_matrix_solve, a, b))
}

setMethod("solve", signature(a = "dgeMatrix", b = "Matrix"),
	  function(a, b, ...)
              .Call(dgeMatrix_matrix_solve, a, as(b, "denseMatrix")))
setMethod("solve", signature(a = "dtrMatrix", b = "Matrix"),
	  function(a, b, ...)
              .Call(dtrMatrix_matrix_solve, a, as(b, "denseMatrix")))
setMethod("solve", signature(a = "dtpMatrix", b = "Matrix"),
	  function(a, b, ...)
              .Call(dtpMatrix_matrix_solve, a, as(b, "denseMatrix")))
setMethod("solve", signature(a = "dsyMatrix", b = "Matrix"),
	  function(a, b, ...)
              .Call(dsyMatrix_matrix_solve, a, as(b, "denseMatrix")))
setMethod("solve", signature(a = "dspMatrix", b = "Matrix"),
	  function(a, b, ...)
              .Call(dspMatrix_matrix_solve, a, as(b, "denseMatrix")))
setMethod("solve", signature(a = "dpoMatrix", b = "Matrix"),
          function(a, b, ...)
              .Call(dpoMatrix_matrix_solve, a, as(b, "denseMatrix")))
setMethod("solve", signature(a = "dppMatrix", b = "Matrix"),
          function(a, b, ...)
              .Call(dppMatrix_matrix_solve, a, as(b, "denseMatrix")))


## ~~~~ sparseMatrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## .... diagonalMatrix .................................................

## Dispatch to methods for ddiMatrix (below)
setMethod("solve", signature(a = "diagonalMatrix", b = "ANY"),
          function(a, b, ...) {a <- ..diag2d(a); callGeneric()})

setMethod("solve", signature(a = "ddiMatrix", b = "missing"),
	  function(a, b, ...) {
              if(a@diag == "N") {
                  ## FIXME?
                  ## if(!all(isN0(a@x)))
                  ##     stop("'a' is exactly singular")
                  a@x <- 1 / a@x
              }
              a@Dimnames <- a@Dimnames[2:1]
	      a
	  })

setMethod("solve", signature(a = "ddiMatrix", b = "numLike"),
	  function(a, b, ...) {
              if(length(b) != a@Dim[2L])
                  stop("dimensions of 'a' and 'b' are incompatible")
              r <-
                  if(a@diag == "N") {
                      ## FIXME?
                      ## if(!all(isN0(a@x)))
                      ##     stop("'a' is exactly singular")
                      as.double(b) / a@x
                  } else as.double(b)
              names(r) <- a@Dimnames[[2L]]
	      .m2ge(r)
	  })

setMethod("solve", signature(a = "ddiMatrix", b = "matrix"),
	  function(a, b, ...) {
              if(dim(b)[1L] != a@Dim[2L])
                  stop("dimensions of 'a' and 'b' are incompatible")
              if(a@diag == "N") {
                  ## FIXME?
                  ## if(!all(isN0(a@x)))
                  ##     stop("'a' is exactly singular")
                  b <- b / a@x
              } else storage.mode(b) <- "double"
              r <- .m2ge(b)
              r@Dimnames <- c(a@Dimnames[2L], r@Dimnames[2L])
              r
	  })

setMethod("solve", signature(a = "ddiMatrix", b = "Matrix"),
	  function(a, b, ...) {
              if(b@Dim[1L] != a@Dim[2L])
                  stop("dimensions of 'a' and 'b' are incompatible")
              if(a@diag == "N") {
                  ## FIXME?
                  ## if(!all(isN0(a@x)))
                  ##     stop("'a' is exactly singular")
                  a@x <- 1 / a@x
              }
              a@Dimnames <- a@Dimnames[2:1]
              a %*% b
	  })


## .... indMatrix ......................................................

## Dispatch to methods for pMatrix (below)  {{ nonsingular <=> permutation }}
setMethod("solve", signature(a = "indMatrix", b = "ANY"),
          function(a, b, ...) {a <- as(a, "pMatrix"); callGeneric()})

setMethod("solve", signature(a = "pMatrix", b = "missing"),
	  function(a, b, ...) t(a))

setMethod("solve", signature(a = "pMatrix", b = "numLike"),
	  function(a, b, ...) {
              r <- as.double(b)[invPerm(a@perm)]
              names(r) <- a@Dimnames[[2L]]
              .m2ge(r)
          })

setMethod("solve", signature(a = "pMatrix", b = "matrix"),
	  function(a, b, ...) {
              storage.mode(b) <- "double"
              r <- b[invPerm(a@perm), , drop = FALSE]
              dimnames(r) <-
                  c(a@Dimnames[2L],
                    if(is.null(dnb <- dimnames(b))) list(NULL) else dnb[2L])
              .m2ge(r)
          })

setMethod("solve", signature(a = "pMatrix", b = "Matrix"),
	  function(a, b, ...) {
              b <- as(b, "dMatrix")
              r <- b[invPerm(a@perm), , drop = FALSE]
              r@Dimnames <- c(a@Dimnames[2L], b@Dimnames[2L])
              r
          })


## .... [CRT]sparseMatrix ..............................................

## Dispatch to methods for d.CMatrix (below)
setMethod("solve", signature(a = "CsparseMatrix", b = "ANY"),
          function(a, b, ...) {
              a <- ..sparse2d(a)
              callGeneric()
          })
setMethod("solve", signature(a = "RsparseMatrix", b = "ANY"),
          function(a, b, ...) {
              a <- ..sparse2d(as(a, "CsparseMatrix"))
              callGeneric()
          })
setMethod("solve", signature(a = "TsparseMatrix", b = "ANY"),
          function(a, b, ...) {
              a <- ..sparse2d(.T2C(a))
              callGeneric()
          })

## TODO: implement triangular solver for RsparseMatrix, so that
##       solve(<d[gt]RMatrix>) can be fast also ... note that for
##       dgRMatrix 'a' we can utilize the LU factorization of t(a)


## .... dgCMatrix ......................................................

## a=dgCMatrix
## b=vector, matrix, or Matrix
## x=dg[Ce]Matrix or double vector ... x=dgCMatrix iff b=sparseMatrix
.solve.dgC.sparse.lu <- function(a, b, tol = .Machine$double.eps) {
    ## MM: see also solveSparse() in ~/R/MM/Pkg-ex/Matrix/Doran-A.R
    n <- (da <- a@Dim)[2L]
    if(da[1L] != n)
        stop("'a' is not square")
    i.v <- is.null(db <- dim(b))
    if((if(i.v) length(b) else db[1L]) != n)
        stop("dimensions of 'a' and 'b' are incompatible")

    lu.a <- lu(a) # error if near-singular

    ## MJ: Is the ratio below really a good proxy for 1/kappa(A) ??

    if(tol > 0) {
        rU <- range(abs(diag(lu.a@U)))
        if(rU[1L] / rU[2L] < tol)
            stop(gettextf("A = LU is computationally singular: min(d)/max(d) = %9.4g, d = abs(diag(U))",
        		  rU[1L] / rU[2L]),
        	 domain = NA)
    }

    ## Solve for  X  in  A X = P' L U Q X = B  ...
    ## 1. compute P B
    pb <- if(i.v) b[lu.a@p + 1L] else b[lu.a@p + 1L, , drop = FALSE]
    ## 2. compute Q X = U^{-1} L^{-1} (P B)
    qx <- solve(lu.a@U, solve(lu.a@L, pb))
    ## 3. compute X = Q' (Q X)
    if(i.v) {
        x <- qx[invPerm(lu.a@q, zero.p = TRUE)]
        names(x) <- a@Dimnames[[2L]]
        .m2ge(x)
    } else if(db[2L] == 1L) {
        ## FIXME: inexplicably, splines:::interpSpline.default() expects
        ##        solve(<dgCMatrix>, <1-column dgeMatrix>, sparse = TRUE)
        ##        to return a dense vector ??  [ as of r82776 ]
        x <- qx[invPerm(lu.a@q, zero.p = TRUE), , drop = TRUE]
        names(x) <- a@Dimnames[[2L]]
        x
    } else {
        x <- qx[invPerm(lu.a@q, zero.p = TRUE), , drop = FALSE]
        x@Dimnames <- c(a@Dimnames[2L],
                        if(is.null(dnb <- dimnames(b))) list(NULL) else dnb[2L])
        x
    }
}

## a=dgCMatrix
## b=vector, matrix, or denseMatrix
## x=dgeMatrix or double vector
.solve.dgC.dense.lu <- function(a, b) {
    n <- (da <- a@Dim)[2L]
    if(da[1L] != n)
        stop("'a' is not square")
    if((if(is.null(db <- dim(b))) length(b) else db[1L]) != n)
        stop("dimensions of 'a' and 'b' are incompatible")

    x <- .Call(dgCMatrix_matrix_solve, a, b, FALSE)
    x@Dimnames <- c(a@Dimnames[2L],
                    if(is.null(dnb <- dimnames(b))) list(NULL) else dnb[2L])
    x
}

setMethod("solve", signature(a = "dgCMatrix", b = "missing"),
	  function(a, b, sparse = NA, ...) {
              x <-
                  if(is.na(sparse) && isSymmetric(a, checkDN = FALSE))
                      solve(forceSymmetric(a), ...)
                  else if(!is.na(sparse) && sparse)
                      solve(a, .sparseDiagonal(a@Dim[2L], shape = "g"),
                            sparse = TRUE, ...)
                  else .solve.dgC.dense.lu(a, diag(a@Dim[2L]))
              x@Dimnames <- a@Dimnames[2:1]
              x
          })

for(.cl in c("numLike", "matrix"))
setMethod("solve", signature(a = "dgCMatrix", b = .cl),
	  function(a, b, sparse = FALSE, ...) {
              if(!is.na(sparse) && sparse)
                  solve(a, .m2sparse(b, "dgC"), sparse = TRUE, ...)
              else .solve.dgC.dense.lu(a, b)
          })
rm(.cl)

setMethod("solve", signature(a = "dgCMatrix", b = "denseMatrix"),
	  function(a, b, sparse = FALSE, ...) {
              if(!is.na(sparse) && sparse)
                  solve(a, ..sparse2d(.sparse2g(as(b, "CsparseMatrix"))),
                        sparse = TRUE, ...)
              else .solve.dgC.dense.lu(a, b)
          })

setMethod("solve", signature(a = "dgCMatrix", b = "sparseMatrix"),
	  function(a, b, sparse = NA, tol = .Machine$double.eps, ...) {
              if(is.na(sparse) && isSymmetric(a, checkDN = FALSE)) {
                  a@Dimnames <- a@Dimnames[c(2L, 2L)]
                  solve(forceSymmetric(a), b, ...)
              } else if(!is.na(sparse) && sparse)
                  .solve.dgC.sparse.lu(a, b, tol = tol)
              else solve(a, as(b, "denseMatrix"), sparse = FALSE, ...)
          })


## .... dsCMatrix ......................................................

## a=dsCMatrix
## b=d[gt]CMatrix
## x=dgCMatrix
.solve.dsC.sparse.chol <- function(a, b, LDL = NA) {
    if(b@Dim[1L] != a@Dim[2L])
        stop("dimensions of 'a' and 'b' are incompatible")

    ## Do _not_ catch warnings directly, as CHOLMOD must free()
    w <- list()
    x <- withCallingHandlers(tryCatch(.Call(dsCMatrix_Csparse_solve, a, b, LDL),
                                      error = identity),
                             warning = function(cond) {
                                 w[[length(w) + 1L]] <<- cond
                                 tryInvokeRestart("muffleWarning")
                             })
    if(.solve.dsC.status(sys.call(), x, w))
        .solve.dgC.sparse.lu(.sparse2g(a), b)
    else {
        x@Dimnames <- c(dimnames(a)[2L], b@Dimnames[2L])
        x
    }
}

## a=dsCMatrix
## b=vector, matrix, or denseMatrix
## x=dgeMatrix or double vector
.solve.dsC.dense.chol <- function(a, b, LDL = NA) {
    if((if(is.null(db <- dim(b))) length(b) else db[1L]) != a@Dim[2L])
        stop("dimensions of 'a' and 'b' are incompatible")

    ## Do _not_ catch warnings directly, as CHOLMOD must free()
    w <- list()
    x <- withCallingHandlers(tryCatch(.Call(dsCMatrix_matrix_solve, a, b, LDL),
                                      error = identity),
                             warning = function(cond) {
                                 w[[length(w) + 1L]] <<- cond
                                 tryInvokeRestart("muffleWarning")
                             })
    if(.solve.dsC.status(sys.call(), x, w))
        .solve.dgC.dense.lu(.sparse2g(a), b)
    else {
        x@Dimnames <- c(dimnames(a)[2L],
                        if(is.null(dnb <- dimnames(b))) list(NULL) else dnb[2L])
        x
    }
}

.solve.dsC.status <- function(call, result, conditionList) {
    e <- inherits(result, "error")
    w <- length(conditionList) > 0L
    if((status <- e || w) && (v <- Matrix.verbose()) >= 1) {
        fmt <- "%s(): Cholesky factorization failed %s... trying LU factorization ..."
        name <- as.character(call[[1L]])
        extras <-
            if(v >= 2)
                paste0(c("with warnings:",
                         if(e) conditionMessage(result),
                         if(w) unlist(lapply(conditionList, conditionMessage),
                                      FALSE, FALSE),
                         ""),
                       collapse = "\n")
            else ""
        message(gettextf(fmt, name, extras), domain = NA)
    }
    status
}

setMethod("solve", signature(a = "dsCMatrix", b = "missing"),
	  function(a, b, LDL = NA, ...) {
              b <- .sparseDiagonal(a@Dim[2L], shape = "g")
              x <- solve(a, b, LDL = LDL, ...)
              x@Dimnames <- dimnames(a)
              x
          })

for(.cl in c("numLike", "matrix", "denseMatrix"))
setMethod("solve", signature(a = "dsCMatrix", b = .cl),
	  function(a, b, LDL = NA, ...)
              .solve.dsC.dense.chol(a, b, LDL = LDL))
rm(.cl)

setMethod("solve", signature(a = "dsCMatrix", b = "sparseMatrix"),
	  function(a, b, LDL = NA, ...)
              .solve.dsC.sparse.chol(
                  a, ..sparse2d(.sparse2g(as(b, "CsparseMatrix"))), LDL = LDL))


## .... dtCMatrix ......................................................

## FIXME: dtCMatrix_sparse_solve() can return an invalid dtCMatrix:
##
## a <- new("dtCMatrix", Dim = c(5L, 5L), diag = "U",
##          p = c(0L, 0L, 0:2, 5L), i = c(1L, 0:3), x = rep(1, 5))
## b <- .trDiagonal(n, unitri = FALSE)
## validObject(.Call(dtCMatrix_sparse_solve, a, b))
##
## This must be fixed at C-level so that we do not rely on .sortCsparse()
## to produce a valid object.

## a=dtCMatrix
## b=d[gt]CMatrix
## x=dgCMatrix
.solve.dtC.sparse <- function(a, b) {
    if(b@Dim[1L] != a@Dim[2L])
        stop("dimensions of 'a' and 'b' are incompatible")
    x <- .sortCsparse(.Call(dtCMatrix_sparse_solve, a, b))
    x@Dimnames <- c(a@Dimnames[2L], b@Dimnames[2L])
    x
}

## a=dtCMatrix
## b=double matrix or dgeMatrix
## x=dgeMatrix
.solve.dtC.dense <- function(a, b) {
    if(dim(b)[1L] != a@Dim[2L])
        stop("dimensions of 'a' and 'b' are incompatible")
    x <- .Call(dtCMatrix_matrix_solve, a, b, isS4(b))
    x@Dimnames <- c(a@Dimnames[2L],
                    if(is.null(dnb <- dimnames(b))) list(NULL) else dnb[2L])
    x
}

setMethod("solve", signature(a = "dtCMatrix", b = "missing"),
	  function(a, b, ...) {
              b <- .trDiagonal(a@Dim[2L], uplo = a@uplo, unitri = FALSE)
              x <- solve(a, b, ...)
              x@Dimnames <- a@Dimnames[2:1]
              x
          })

setMethod("solve", signature(a = "dtCMatrix", b = "numLike"),
	  function(a, b, ...) {
              b <- cbind(as.double(b), deparse.level = 0L)
              .solve.dtC.dense(a, b)
          })

setMethod("solve", signature(a = "dtCMatrix", b = "matrix"),
	  function(a, b, ...) {
              storage.mode(b) <- "double"
              .solve.dtC.dense(a, b)
	  })

setMethod("solve", signature(a = "dtCMatrix", b = "sparseMatrix"),
	  function(a, b, ...)
              solve(a, ..sparse2d(as(b, "CsparseMatrix")), ...))

setMethod("solve", signature(a = "dtCMatrix", b = "dgCMatrix"),
	  function(a, b, ...)
              .solve.dtC.sparse(a, b))

setMethod("solve", signature(a = "dtCMatrix", b = "dtCMatrix"),
	  function(a, b, ...) {
              x <- .solve.dtC.sparse(a, b)
              if(a@uplo != b@uplo)
                  x
              else if(a@uplo == "U")
                  triu(x)
              else tril(x)
          })

setMethod("solve", signature(a = "dtCMatrix", b = "dsCMatrix"),
	  function(a, b, ...)
              .solve.dtC.sparse(a, .sparse2g(b)))

setMethod("solve", signature(a = "dtCMatrix", b = "denseMatrix"),
	  function(a, b, ...)
              solve(a, ..dense2d(b), ...))

setMethod("solve", signature(a = "dtCMatrix", b = "dgeMatrix"),
	  function(a, b, ...)
              .solve.dtC.dense(a, b))

setMethod("solve", signature(a = "dtCMatrix", b = "dtrMatrix"),
	  function(a, b, ...) {
              x <- .solve.dtC.dense(a, .dense2g(b))
              if(a@uplo != b@uplo)
                  x
              else if(a@uplo == "U")
                  triu(x)
              else tril(x)
          })

setMethod("solve", signature(a = "dtCMatrix", b = "dtpMatrix"),
	  function(a, b, ...) {
              x <- .solve.dtC.dense(a, .dense2g(b))
              if(a@uplo != b@uplo)
                  x
              else .Call(unpackedMatrix_pack, x, TRUE, TRUE, a@uplo == "U")
          })

setMethod("solve", signature(a = "dtCMatrix", b = "dsyMatrix"),
	  function(a, b, ...)
              .solve.dtC.dense(a, .dense2g(b)))

setMethod("solve", signature(a = "dtCMatrix", b = "dspMatrix"),
	  function(a, b, ...)
              .solve.dtC.dense(a, .dense2g(b)))


## ~~~~ MatrixFactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("solve", signature(a = "MatrixFactorization", b = "ANY"),
	  function(a, b, ...)
              .bail.out.2("solve", class(a), class(b)))

setMethod("solve", signature(a = "MatrixFactorization", b = "missing"),
	  function(a, b, ...)
              solve(a, .sparseDiagonal(a@Dim[2L], shape = "g"), ...))


## .... CHMfactor ......................................................

setMethod("solve", signature(a = "CHMfactor", b = "missing"),
	  function(a, b, system = c("A", "LDLt",
                                    "LD", "DLt", "L", "Lt", "D",
                                    "P", "Pt"),
		   ...) {
	      system.def <- eval(formals()$system)
              system <- match.arg(system, system.def)
              b <- .sparseDiagonal(a@Dim[2L], shape = "g")
              x <- .Call(CHMfactor_spsolve,
                         a, b, match(system, system.def, 0L))
              switch(system,
                     A =, LDLt = .M2symm(x),
                     LD =, DLt =, L =, Lt =, D = .M2tri(x),
                     P =, Pt = as(x, "pMatrix"))
	  })

for(.cl in c("numLike", "matrix", "denseMatrix"))
setMethod("solve", signature(a = "CHMfactor", b = .cl),
	  function(a, b, system = c("A", "LDLt",
                                    "LD", "DLt", "L", "Lt", "D",
                                    "P", "Pt"),
                   ...) {
              system.def <- eval(formals()$system)
              system <- match.arg(system, system.def)
              .Call(CHMfactor_solve,
                    a, b, match(system, system.def, 0L))
          })

setMethod("solve", signature(a = "CHMfactor", b = "sparseMatrix"),
	  function(a, b, system = c("A", "LDLt",
                                    "LD", "DLt", "L", "Lt", "D",
                                    "P", "Pt"),
		   ...) {
               system.def <- eval(formals()$system)
               system <- match.arg(system, system.def)
               b <- ..sparse2d(.sparse2g(as(b, "CsparseMatrix")))
               ## cholmod_spsolve() in ../src/CHOLMOD/Cholesky/cholmod_spsolve.c
               .Call(CHMfactor_spsolve,
                     a, b, match(system, system.def, 0L))
	  })


## .... denseLU ........................................................

setMethod("solve", signature(a = "denseLU", b = "missing"),
	  function(a, b, ...) {
	      ll <- expand(a)
	      solve(ll$U, solve(ll$L, ll$P))
	  })


## .... sparseLU .......................................................

## TODO?


## .... sparseQR .......................................................

setMethod("solve", signature(a = "sparseQR", b = "ANY"),
	  function(a, b, ...)
              qr.coef(a, b))

setMethod("solve", signature(a = "sparseQR", b = "missing"),
	  function(a, b, ...)
              qr.coef(a, .sparseDiagonal(a@Dim[2L], shape = "g")))


## .... Schur ..........................................................

## TODO?


## ~~~~ 'Matrix' class on RHS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("solve", signature(a = "matrix", b = "Matrix"),
	  function(a, b, ...) solve(as(a, "dMatrix"), b, ...))

setMethod("solve", signature(a = "matrix", b = "sparseVector"),
	  function(a, b, ...) solve(as(a, "dMatrix"), as(b, "sparseMatrix"), ...))

setMethod("solve", signature(a = "Matrix", b = "sparseVector"),
	  function(a, b, ...) solve(a, as(b, "sparseMatrix"), ...))

setMethod("solve", signature(a = "MatrixFactorization", b = "sparseVector"),
	  function(a, b, ...) solve(a, as(b, "sparseMatrix"), ...))


## ~~~~ Exported solvers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## a=dgCMatrix
## b=vector, matrix, or Matrix
## x=d[gt][^RT]Matrix or double vector
.solve.dgC.lu <- function(a, b, tol = .Machine$double.eps, check = TRUE) {
    if(check && !is(a, "dgCMatrix"))
        a <- as(as(as(a, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    .solve.dgC.sparse.lu(a, b, tol = tol)
}

## a=dgCMatrix
## b=vector or 1-column matrix
## x=dgeMatrix
.solve.dgC.chol <- function(x, y, check = TRUE) { # -> MatrixModels
    if(check && !is(x, "dgCMatrix"))
        x <- as(as(as(x, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    .Call(dgCMatrix_cholsol, x, y)
}

## a=dgCMatrix
## b=vector or 1-column matrix
## x=dgeMatrix
.solve.dgC.qr <- function(x, y, order = 1L, check = TRUE) { # -> MatrixModels
    if(check && !is(x, "dgCMatrix"))
        x <- as(as(as(x, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    .Call(dgCMatrix_qrsol, x, y, order)
}
