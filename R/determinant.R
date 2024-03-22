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

setMethod("determinant", c(x = "MatrixFactorization", logarithm = "missing"),
          function(x, logarithm = TRUE, ...)
              determinant(x, TRUE, ...))

## FIXME: if we knew the specific class of 'T', then we could optimize
## knowing that 'T' is block upper triangular with 1-by-1 and 2-by-2
## diagonal blocks
setMethod("determinant", c(x = "Schur", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              determinant(x@T, logarithm, ...))

setMethod("determinant", c(x = "denseLU", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(denseLU_determinant, x, logarithm))

setMethod("determinant", c(x = "sparseLU", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(sparseLU_determinant, x, logarithm))

setMethod("determinant", c(x = "sparseQR", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(sparseQR_determinant, x, logarithm))

for(.cl in c("BunchKaufman", "pBunchKaufman"))
setMethod("determinant", c(x = .cl, logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(BunchKaufman_determinant, x, logarithm))
rm(.cl)

for(.cl in c("Cholesky", "pCholesky"))
setMethod("determinant", c(x = .cl, logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(Cholesky_determinant, x, logarithm))
rm(.cl)

setMethod("determinant", c(x = "CHMfactor", logarithm = "logical"),
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
                      warning(gettextf("the default value of argument '%s' of method '%s(<%s>, <%s>)' may change from %s to %s as soon as the next release of Matrix; set '%s' when programming",
                                       "sqrt", "determinant", "CHMfactor", "logical", "TRUE", "FALSE", "sqrt"),
                              domain = NA)
                  }
              }
              .Call(CHMfactor_determinant, x, logarithm, sqrt)
          })


########################################################################
##  2. Matrix
########################################################################

setMethod("determinant", c(x = "Matrix", logarithm = "missing"),
          function(x, logarithm = TRUE, ...)
              determinant(x, TRUE, ...))

setMethod("determinant", c(x = "Matrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              determinant(.M2kind(x, ","), logarithm, ...))

## .... GENERAL ........................................................

setMethod("determinant", c(x = "dgeMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              d <- x@Dim
              if(d[1L] != d[2L])
                  stop("determinant of non-square matrix is undefined")
              trf <- lu(x, warnSing = FALSE)
              determinant(trf, logarithm, ...)
          })

setMethod("determinant", c(x = "dgCMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              d <- x@Dim
              if(d[1L] != d[2L])
                  stop("determinant of non-square matrix is undefined")
              trf <- lu(x, errSing = FALSE)
              if(isS4(trf))
                  determinant(trf, logarithm, ...)
              else .mkDet(if(anyNA(x@x)) NaN else -Inf, logarithm, 1L)
          })

setMethod("determinant", c(x = "dgRMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              determinant(.tCRT(x), logarithm, ...))

setMethod("determinant", c(x = "dgTMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              determinant(.M2C(x), logarithm, ...))

setMethod("determinant", c(x = "indMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              d <- x@Dim
              if(d[1L] != d[2L])
                  stop("determinant of non-square matrix is undefined")
              if(anyDuplicated.default(perm <- x@perm))
                  .mkDet(-Inf, logarithm, 1L)
              else .mkDet(0, logarithm, signPerm(perm))
          })

setMethod("determinant", c(x = "pMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .mkDet(0, logarithm, signPerm(x@perm)))


## .... SYMMETRIC ......................................................

for(.cl in c("dsyMatrix", "dspMatrix"))
setMethod("determinant", c(x = .cl, logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              trf <- BunchKaufman(x, warnSing = FALSE)
              determinant(trf, logarithm, ...)
          })
rm(.cl)

for(.cl in c("dpoMatrix", "dppMatrix"))
setMethod("determinant", c(x = .cl, logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              trf <- tryCatch(
                  Cholesky(x, perm = FALSE),
                  error = function(e) BunchKaufman(x, warnSing = FALSE))
              determinant(trf, logarithm, ...)
          })
rm(.cl)

setMethod("determinant", c(x = "dsCMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              trf <- tryCatch(
                  Cholesky(x, perm = TRUE, LDL = TRUE, super = FALSE),
                  error = function(e) lu(x, errSing = FALSE))
              if(isS4(trf))
                  determinant(trf, logarithm, sqrt = FALSE, ...)
              else .mkDet(if(anyNA(x@x)) NaN else -Inf, logarithm, 1L)
          })

setMethod("determinant", c(x = "dsRMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              determinant(.tCRT(x), logarithm, ...))

setMethod("determinant", c(x = "dsTMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              determinant(.M2C(x), logarithm, ...))


## .... TRIANGULAR .....................................................

for(.cl in c("triangularMatrix", "diagonalMatrix"))
setMethod("determinant", c(x = .cl, logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              if(x@diag == "N")
                  .mkDet(x = diag(x, names = FALSE), logarithm = logarithm)
              else .mkDet(0, logarithm, 1L))
rm(.cl)
