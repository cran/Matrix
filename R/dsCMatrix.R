#### Symmetric Sparse Matrices in compressed column-oriented format

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
setAs("matrix", "dsCMatrix",
      function(from) as(.m2dgC(from), "symmetricMatrix"))
} ## MJ

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
setAs("dgCMatrix", "dsCMatrix",
      function(from) { ## r2130 ... | 2008-03-14 | added deprecation warning
	  warning("as(.,\"dsCMatrix\") is deprecated (since 2008); do use as(., \"symmetricMatrix\")")
	  as(from, "symmetricMatrix")
      })

setAs("dsCMatrix", "dgTMatrix", # needed for show(), image()
      function(from)
      ## pre-Cholmod -- FIXME: get rid of
      .Call(dsCMatrix_to_dgTMatrix, from))

setAs("dsCMatrix", "dgeMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgeMatrix"))

setAs("dsCMatrix", "matrix",
      function(from) as(as(from, "generalMatrix"), "matrix"))

setAs("dsCMatrix", "lsCMatrix",
      function(from) new("lsCMatrix", i = from@i, p = from@p, uplo = from@uplo,
                         x = as.logical(from@x),
                         Dim = from@Dim, Dimnames = from@Dimnames))
setAs("dsCMatrix", "nsCMatrix",
      function(from) new("nsCMatrix", i = from@i, p = from@p, uplo = from@uplo,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("dsCMatrix", "dgCMatrix",
      function(from) .Call(Csparse_symmetric_to_general, from))

setAs("dsCMatrix", "dsyMatrix",
      function(from) as(from, "denseMatrix"))
} ## MJ

## MJ: unused
if(FALSE) {
##' Check if \code{name} (== "[sS][pP][dD]Cholesky") fits the values of the
##' logicals (perm, LDL, super).
##' @param name a string such as "sPdCholesky"
##' @param perm also known as \code{pivot}
##' @param LDL
##' @param super
##' @return logical: TRUE if the name matches
.chkName.CHM <- function(name, perm, LDL, super)
    .Call(R_chkName_Cholesky, name, perm, LDL, super)
## ../src/dsCMatrix.c

.CHM.factor.name <- function(perm, LDL, super)
    .Call(R_chm_factor_name, perm, LDL, super)
} ## MJ

## MJ: no longer needed ... methods now inherited from CsparseMatrix
if(FALSE) {
## have rather tril() and triu() methods than
## setAs("dsCMatrix", "dtCMatrix", ....)
setMethod("tril", "dsCMatrix",
	  function(x, k = 0, ...) {
	      if(x@uplo == "L" && k == 0)
		  ## same internal structure (speedup potential !?)
		  new("dtCMatrix", uplo = x@uplo, i = x@i, p = x@p,
		      x = x@x, Dim = x@Dim, Dimnames = x@Dimnames)
	      else tril(.sparse2g(x), k = k, ...)
	  })

setMethod("triu", "dsCMatrix",
	  function(x, k = 0, ...) {
	      if(x@uplo == "U" && k == 0)
		  ## same internal structure (speedup potential !?)
		  new("dtCMatrix", uplo = x@uplo, i = x@i, p = x@p,
		      x = x@x, Dim = x@Dim, Dimnames = x@Dimnames)
	      else triu(.sparse2g(x), k = k, ...)
	  })
} ## MJ

## MJ: no longer needed ... method now inherited from CsparseMatrix
if(FALSE) {
setMethod("t", signature(x = "dsCMatrix"),
          function(x) .Call(Csparse_transpose, x, FALSE),
          valueClass = "dsCMatrix")
} ## MJ

### These two are very similar, the first one has the advantage
### to be applicable to 'Chx' directly:

## "used" currently only in ../tests/factorizing.R
.diag.dsC <- function(x, Chx = Cholesky(x, LDL=TRUE), res.kind = "diag") {
    force(Chx)
    if(!missing(Chx)) stopifnot(.isLDL(Chx), is.integer(Chx@p), is.double(Chx@x))
    .Call(diag_tC, Chx, res.kind)
    ##    ^^^^^^^ from ../src/Csparse.c
    ## => res.kind in ("trace", "sumLog", "prod", "min", "max", "range", "diag", "diagBack")
}

## MJ: unused
if(FALSE) {
## here, we  *could* allow a 'mult = 0' factor :
.CHM.LDL.D <- function(x, perm = TRUE, res.kind = "diag") {
    .Call(dsCMatrix_LDL_D, x, perm, res.kind)
    ##    ^^^^^^^^^^^^^^^^ from ../src/dsCMatrix.c
}
} ## MJ

## FIXME:  kind = "diagBack" is not yet implemented
##	would be much more efficient, but there's no CHOLMOD UI (?)
##
## Note: for det(), permutation is unimportant;
##       for diag(), apply *inverse* permutation
##    	q <- p ; q[q] <- seq_along(q); q

## setMethod("writeHB", signature(obj = "dsCMatrix"),
##           function(obj, file, ...) {
##               .Deprecated("writeMM")
##               .Call(Matrix_writeHarwellBoeing,
##                     if (obj@uplo == "U") t(obj) else obj,
##                     as.character(file), "DSC")
##           })
