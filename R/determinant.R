## METHODS FOR GENERIC: determinant
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Constructor for "det" objects, used liberally below
.mkDet <- function(modulus = sum(log(abs(x))),
                   logarithm = TRUE,
                   sign = if(prod(base::sign(x)) < 0) -1L else 1L,
                   x) {
    if(!logarithm)
        modulus <- exp(modulus)
    attr(modulus, "logarithm") <- logarithm
    val <- list(modulus = modulus, sign = sign)
    class(val) <- "det"
    val
}

## 'base::det' calls 'base::determinant', which is not S4 generic,
## so we define own our 'det' calling 'Matrix::determinant' ...
det <- base::det
environment(det) <- environment() # the Matrix namespace


########################################################################
##  1. MatrixFactorization
########################################################################

setMethod("determinant", signature(x = "MatrixFactorization", logarithm = "missing"),
          function(x, logarithm = TRUE, ...)
              determinant(x, TRUE, ...))

## FIXME: if we knew the specific class of 'T', then we could optimize
## knowing that 'T' is block upper triangular with 1-by-1 and 2-by-2
## diagonal blocks
setMethod("determinant", signature(x = "Schur", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              determinant(x@T, logarithm, ...))

setMethod("determinant", signature(x = "denseLU", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(denseLU_determinant, x, logarithm))

setMethod("determinant", signature(x = "sparseLU", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(sparseLU_determinant, x, logarithm))

setMethod("determinant", signature(x = "sparseQR", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(sparseQR_determinant, x, logarithm))

setMethod("determinant", signature(x = "BunchKaufman", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(BunchKaufman_determinant, x, logarithm, FALSE))

setMethod("determinant", signature(x = "pBunchKaufman", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(BunchKaufman_determinant, x, logarithm, TRUE))

setMethod("determinant", signature(x = "Cholesky", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(Cholesky_determinant, x, logarithm, FALSE))

setMethod("determinant", signature(x = "pCholesky", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(Cholesky_determinant, x, logarithm, TRUE))

setMethod("determinant", signature(x = "CHMfactor", logarithm = "logical"),
          function(x, logarithm = TRUE, sqrt = TRUE, ...) {
              if(missing(sqrt)) {
                  w <- getOption("Matrix.warnSqrtDefault",
                                 .MatrixEnv[["warnSqrtDefault"]])
                  if(is.atomic(w) && length(w) == 1L &&
                     ((w.na <- is.na(w <- as.integer(w))) || w > 0L)) {
                      if(w.na)
                          on.exit(options(Matrix.warnSqrtDefault = 0L))
                      else if(w > 1L) {
                          oop <- options(warn = 2L)
                          on.exit(options(oop))
                      }
                      warning("the default value of argument 'sqrt' of method 'determinant(<CHMfactor>, <logical>)' may change from TRUE to FALSE as soon as the next release of Matrix; set 'sqrt' when programming")
                  }
              }
              .Call(CHMfactor_determinant, x, logarithm, sqrt)
          })


########################################################################
##  2. Matrix
########################################################################

setMethod("determinant", signature(x = "Matrix", logarithm = "missing"),
          function(x, logarithm = TRUE, ...)
              determinant(x, TRUE, ...))

setMethod("determinant", signature(x = "Matrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              determinant(as(x, "dMatrix"), logarithm, ...))

## .... GENERAL ........................................................

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              d <- x@Dim
              if(d[1L] != d[2L])
                  stop("determinant of non-square matrix is undefined")
              trf <- lu(x, warnSing = FALSE)
              determinant(trf, logarithm, ...)
          })

setMethod("determinant", signature(x = "dgCMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              d <- x@Dim
              if((n <- d[1L]) != d[2L])
                  stop("determinant of non-square matrix is undefined")
              trf <- lu(x, errSing = FALSE)
              if(isS4(trf))
                  determinant(trf, logarithm, ...)
              else .mkDet(if(anyNA(x@x)) NaN else -Inf, logarithm, 1L)
          })

setMethod("determinant", signature(x = "dgRMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              determinant(.tCRT(x), logarithm, ...))

setMethod("determinant", signature(x = "dgTMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              determinant(.M2C(x), logarithm, ...))

setMethod("determinant", signature(x = "indMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              d <- x@Dim
              if((n <- d[1L]) != d[2L])
                  stop("determinant of non-square matrix is undefined")
              if(anyDuplicated.default(perm <- x@perm))
                  .mkDet(-Inf, logarithm, 1L)
              else .mkDet(0, logarithm, signPerm(perm))
          })

setMethod("determinant", signature(x = "pMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .mkDet(0, logarithm, signPerm(x@perm)))


## .... TRIANGULAR .....................................................

setMethod("determinant", signature(x = "triangularMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              if(x@diag == "N")
                  .mkDet(x = diag(x, names = FALSE), logarithm = logarithm)
              else .mkDet(0, logarithm, 1L)
          })

setMethod("determinant", signature(x = "diagonalMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              if(x@diag == "N")
                  .mkDet(x = x@x, logarithm = logarithm)
              else .mkDet(0, logarithm, 1L)
          })


## .... SYMMETRIC ......................................................

for(.cl in c("dsyMatrix", "dspMatrix"))
setMethod("determinant", signature(x = .cl, logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              trf <- BunchKaufman(x, warnSing = FALSE)
              determinant(trf, logarithm, ...)
          })
rm(.cl)

for(.cl in c("dpoMatrix", "dppMatrix"))
setMethod("determinant", signature(x = .cl, logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              trf <- tryCatch(
                  Cholesky(x, perm = FALSE),
                  error = function(e) BunchKaufman(x, warnSing = FALSE))
              determinant(trf, logarithm, ...)
          })
rm(.cl)

setMethod("determinant", signature(x = "dsCMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              trf <- tryCatch(
                  Cholesky(x, perm = TRUE, LDL = TRUE, super = FALSE),
                  error = function(e) lu(x, errSing = FALSE))
              if(isS4(trf))
                  determinant(trf, logarithm, sqrt = FALSE, ...)
              else .mkDet(if(anyNA(x@x)) NaN else -Inf, logarithm, 1L)
          })

setMethod("determinant", signature(x = "dsRMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              determinant(.tCRT(x), logarithm, ...))

setMethod("determinant", signature(x = "dsTMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              determinant(.M2C(x), logarithm, ...))

## MJ: unused
if(FALSE) {
ldet1.dsC <- function(x, ...)
    .Call(CHMfactor_ldetL2, Cholesky(x, ...))

## ~3% faster than ldet1:
ldet2.dsC <- function(x, ...) {
    Ch <- Cholesky(x, super = FALSE, ...)
    .Call(diag_tC, Ch, "sumLog")
}

## <1% faster than ldet2:
ldet3.dsC <- function(x, perm = TRUE)
    .Call(dsCMatrix_LDL_D, x, perm = perm, "sumLog")
} ## MJ
