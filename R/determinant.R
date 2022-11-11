## METHODS FOR GENERIC: determinant
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Constructor for "det" objects, used liberally below
.mkDet <- function(d, logarithm = TRUE,
                   ldet = sum(log(abs(d))),
                   sig = if(prod(sign(d)) < 0) -1L else 1L) {
    ##             ^^^ -1 or 1, never 0
    modulus <- if(logarithm) ldet else exp(ldet)
    attr(modulus, "logarithm") <- logarithm
    val <- list(modulus = modulus, sign = sig)
    class(val) <- "det"
    val
}


## ~~~~ GENERAL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Compute product of determinants in LU factorization

.det.dge <- function(x, logarithm, ...)
    .Call(dgeMatrix_determinant, x, logarithm)

.det.dgC <- function(x, logarithm, ...) {
    d <- x@Dim
    if((n <- d[1L]) != d[2L])
        stop("determinant of non-square matrix is undefined")
    if(n <= 1L)
        return(.mkDet(diag(x, names = FALSE), logarithm))
    ll <- lu(x, errSing = FALSE)
    if(identical(ll, NA))
        ## LU factorization failed due to singularity
	return(.mkDet(ldet = if(anyNA(x)) NaN else -Inf,
                      logarithm = logarithm, sig = 1L))
    r <- .mkDet(diag(ll@U), logarithm)
    ## det(x)
    ## = det(P L U Q)
    ## = det(P) * 1 * det(U) * det(Q)
    ## where det(P), det(Q) are in {-1,1}
    r$sign <- r$sign * signPerm(ll@p + 1L) * signPerm(ll@q + 1L)
    r
}


## ~~~~ TRIANGULAR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Compute product of diagonal elements

.det.tri <- function(x, logarithm, ...) {
    if(x@diag == "N")
        .mkDet(diag(x, names = FALSE), logarithm)
    else .mkDet(, logarithm, ldet = 0, sig = 1L)
}

.det.diag <- function(x, logarithm, ...) {
    if(x@diag == "N")
        .mkDet(x@x, logarithm)
    else .mkDet(, logarithm, ldet = 0, sig = 1L)
}


## ~~~~ SYMMETRIC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Compute product of determinants in Bunch-Kaufman factorization

.det.dsy <- function(x, logarithm, ...)
    .Call(dsyMatrix_determinant, x, logarithm)

.det.dsp <- function(x, logarithm, ...)
    .Call(dspMatrix_determinant, x, logarithm)

## Compute product of determinants in Cholesky factorization

.det.dpo <- function(x, logarithm, ...) {
    cholx <- .Call(dpoMatrix_trf, x, 2L)
    .mkDet(, logarithm, ldet = 2 * sum(log(abs(diag(cholx)))), sig = 1L)
}

.det.dpp <- function(x, logarithm, ...) {
    cholx <- .Call(dppMatrix_trf, x, 2L)
    .mkDet(, logarithm, ldet = 2 * sum(log(abs(diag(cholx)))), sig = 1L)
}

.det.dpC <- function(x, logarithm, ...) {
    cholx <- .Call(dsCMatrix_Cholesky, x, TRUE, TRUE, FALSE, 0)
    .mkDet(.Call(diag_tC, cholx, res.kind = "diag"), logarithm)
}

## Compute product of determinants in Cholesky factorization;
## if that fails due to non-positive definiteness, then go via "general"

.det.dsC <- function(x, logarithm, ...) {
    if(x@Dim[1L] <= 1L)
        .mkDet(diag(x, names = FALSE), logarithm)
    else
        tryCatch(suppressWarnings(.det.dpC(x, logarithm)),
                 error = function(e) .det.dgC(.sparse2g(x), logarithm))
}


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("determinant",
          signature(x = "Matrix", logarithm = "missing"),
          function(x, logarithm, ...) determinant(x, TRUE, ...))

setMethod("determinant",
          signature(x = "MatrixFactorization", logarithm = "missing"),
          function(x, logarithm, ...) determinant(x, TRUE, ...))

setMethod("determinant",
          signature(x = "Matrix", logarithm = "logical"),
          function(x, logarithm, ...) determinant(as(x, "dMatrix"), logarithm, ...))

setMethod("determinant",
          signature(x = "triangularMatrix", logarithm = "logical"),
          .det.tri)

setMethod("determinant",
          signature(x = "diagonalMatrix", logarithm = "logical"),
	  .det.diag)

setMethod("determinant",
          signature(x = "pMatrix", logarithm = "logical"),
	  function(x, logarithm, ...) {
              if(x@Dim[1L] <= 1L)
                  .mkDet(diag(x, names = FALSE), logarithm)
              else .mkDet(, logarithm, ldet = 0, sig = signPerm(x@perm))
	  })

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "logical"),
          .det.dge)

setMethod("determinant", signature(x = "dsyMatrix", logarithm = "logical"),
          .det.dsy)

setMethod("determinant", signature(x = "dspMatrix", logarithm = "logical"),
          .det.dsp)

setMethod("determinant", signature(x = "dpoMatrix", logarithm = "logical"),
          .det.dpo)

setMethod("determinant", signature(x = "dppMatrix", logarithm = "logical"),
          .det.dpp)

setMethod("determinant", signature(x = "dgCMatrix", logarithm = "logical"),
          .det.dgC)

setMethod("determinant", signature(x = "dgRMatrix", logarithm = "logical"),
          function(x, logarithm, ...) .det.dgC(.tCR2RC(x), logarithm))

setMethod("determinant", signature(x = "dgTMatrix", logarithm = "logical"),
          function(x, logarithm, ...) .det.dgC(.T2C(x), logarithm))

setMethod("determinant", signature(x = "dsCMatrix", logarithm = "logical"),
          .det.dsC)

setMethod("determinant", signature(x = "dsRMatrix", logarithm = "logical"),
          function(x, logarithm, ...) .det.dsC(.tCR2RC(x), logarithm))

setMethod("determinant", signature(x = "dsTMatrix", logarithm = "logical"),
          function(x, logarithm, ...) .det.dsC(.T2C(x), logarithm))

setMethod("determinant",
          signature(x = "CHMfactor", logarithm = "logical"),
          function(x, logarithm, ...)
              .mkDet(, logarithm,
                     ldet = 0.5 * .Call(CHMfactor_ldetL2, x), sig = 1L))

rm(.det.tri, .det.diag, .det.dsy, .det.dsp)

if(.Matrix.supporting.cached.methods) {
mkDet <- .mkDet
}

## MJ: used only in tests
if(TRUE) {
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

## MJ: no longer needed ... replacement above
if(FALSE) {
setMethod("determinant", signature(x = "dsCMatrix", logarithm = "logical"),
	  function(x, logarithm, ...) {
	  if(x@Dim[1] <= 1L)
	      return(mkDet(diag(x), logarithm))
	  Chx <- tryCatch(suppressWarnings(Cholesky(x, LDL=TRUE)),
                          error = function(e) NULL)
	  ## or
	  ## ldet <- .Call("CHMfactor_ldetL2", Chx) # which would also work
	  ##				     when Chx <- Cholesky(x, super=TRUE)
          ## ldet <- tryCatch(.Call(dsCMatrix_LDL_D, x, perm=TRUE, "sumLog"),
	  ## if(is.null(ldet))

          if(is.null(Chx))  ## we do *not* have a positive definite matrix
	      detSparseLU(x, logarithm)
	  else {
              d <- .Call(diag_tC, Chx, res.kind = "diag")
	      mkDet(d, logarithm=logarithm)
          }
      })
} ## MJ
