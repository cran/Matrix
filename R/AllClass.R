## README:
##
## Validity methods should assume that methods for superclasses have passed,
## following validObject(). We should _not_ be testing, e.g., length(Dim),
## typeof(Dimnames), etc. repeatedly ...
##
## When checking whether a class is validated correctly, take care to follow
## the 'contains' recursively!!


## To be used in initialize() method for Matrix, or other constructors
## NB: This must be defined _here_ and _not_ be migrated to ./Auxiliaries.R
fixupDN <- function(dn) .Call(R_DimNames_fixup, dn)

## MJ: no longer
if(FALSE) {
.fixupDimnames <- function(dnms) {
    N.N <- list(NULL, NULL)
    if(is.null(dnms) || identical(dnms, N.N)) return(N.N)
    ## else
    if(any(i0 <- lengths(dnms) == 0) && !all(vapply(dnms[i0], is.null, NA)))
	## replace character(0) etc, by NULL :
	dnms[i0] <- list(NULL)
    ## coerce, e.g. integer dimnames to character: -- as  R's matrix(..):
    if(any(i0 <- vapply(dnms, function(d) !is.null(d) && !is.character(d), NA)))
	dnms[i0] <- lapply(dnms[i0], as.character)
    dnms
}
} ## MJ


########################################################################
##  1. Matrix
########################################################################

## ====== Virtual Subclasses ===========================================

## ------ The Mother Class 'Matrix' ------------------------------------

## Virtual class of all Matrix objects
setClass("Matrix", contains = "VIRTUAL",
	 slots = c(Dim = "integer", Dimnames = "list"),
	 prototype = prototype(Dim = integer(2L), Dimnames = list(NULL, NULL)),
	 validity = function(object) .Call(Matrix_validate, object))

## Matrix_validate() allows Dimnames[[i]] to be a vector of type
## other than "character" and, moreover, to be a vector of length
## zero rather than NULL.  fixupDN() takes care of the coercions.
setMethod("initialize", "Matrix",
          function(.Object, ...) {
              .Object <- callNextMethod()
              ## Suboptimal if ...names() is NULL but that will "never"
              ## happen if ...length() is nonzero:
              if(...length() && any(...names() == "Dimnames"))
                  .Object@Dimnames <- fixupDN(.Object@Dimnames)
              .Object
          })

if(FALSE) {
## This method would allow 'Dimnames' to (possibly) define 'Dim'.
## However, DimNames_validate() in ../src/validity.c and the way
## it is used in Matrix_validate() would need to change.
setMethod("initialize", "Matrix",
          function(.Object, ...) {
              .Object <- callNextMethod()
              if(...length() && any((nms <- ...names()) == "Dimnames")) {
                  .Object@Dimnames <- DN <- fixupDN(.Object@Dimnames)
                  if(!(any(nms == "Dim") ||
                       is.null(DN[[1L]]) || is.null(DN[[2L]])))
                      .Object@Dim <- lengths(DN, use.names = FALSE)
              }
              .Object
          })
}


## ------ Virtual by structure -----------------------------------------

## Virtual class of composite matrices,
## i.e., those for which it makes sense to define a factorization
setClass("compMatrix", contains = c("Matrix", "VIRTUAL"),
	 slots = c(factors = "list"),
         validity = function(object) .Call(compMatrix_validate, object))

## Virtual class of matrices that are not symmetric, triangular, _or diagonal_
setClass("generalMatrix", contains = c("compMatrix", "VIRTUAL"))

## Virtual class of triangular matrices
setClass("triangularMatrix", contains = c("Matrix", "VIRTUAL"),
	 slots = c(uplo = "character", diag = "character"),
	 prototype = prototype(uplo = "U", diag = "N"),
	 validity = function(object) .Call(triangularMatrix_validate, object))

## Virtual class of symmetric matrices
setClass("symmetricMatrix", contains = c("compMatrix", "VIRTUAL"),
	 slots = c(uplo = "character"),
	 prototype = prototype(uplo = "U"),
	 validity = function(object) .Call(symmetricMatrix_validate, object))


## ------ Virtual by kind ----------------------------------------------

## Virtual class of _n_onzero pattern matrices
## NB: only subclass ndenseMatrix requires an 'x' slot
setClass("nMatrix", contains = c("Matrix", "VIRTUAL"))

## Virtual class of logical matrices,
## * typically the result of comparisons, e.g., <dMatrix> <relop> <dMatrix>,
##   hence NA are allowed and distinct from TRUE, in contrast with nMatrix
setClass("lMatrix", contains = c("Matrix", "VIRTUAL"),
         slots = c(x = "logical"),
         validity = function(object) .Call(lMatrix_validate, object))

## Virtual class of integer matrices
setClass("iMatrix", contains = c("Matrix", "VIRTUAL"),
         slots = c(x = "integer"),
         validity = function(object) .Call(iMatrix_validate, object))

## Virtual class of double matrices
setClass("dMatrix", contains = c("Matrix", "VIRTUAL"),
         slots = c(x = "numeric"),
         validity = function(object) .Call(dMatrix_validate, object))

## Virtual class of complex matrices
## * initial 'z' is derived from the names of LAPACK routines
setClass("zMatrix", contains = c("Matrix", "VIRTUAL"),
         slots = c(x = "complex"),
         validity = function(object) .Call(zMatrix_validate, object))


## ------ Virtual Dense ------------------------------------------------

## Virtual class of dense matrices
## * includes "unpacked" _and_ "packed" matrices
## * included diagonal matrices until 0.999375-11 (2008-07)
setClass("denseMatrix", contains = c("Matrix", "VIRTUAL"))


## ...... Virtual Dense ... by storage .................................

## Virtual class of dense, "unpacked" matrices, s.t. length(.@x) == m*n
setClass("unpackedMatrix", contains = c("denseMatrix", "VIRTUAL"),
         validity = function(object) .Call(unpackedMatrix_validate, object))

## Virtual class of dense, "packed" matrices, s.t. length(.@x) == n*(n+1)/2
setClass("packedMatrix", contains = c("denseMatrix", "VIRTUAL"),
         slots = c(uplo = "character"),
         prototype = prototype(uplo = "U"),
	 validity = function(object) .Call(packedMatrix_validate, object))


## ...... Virtual Dense ... by kind ....................................

## Virtual class of dense, _n_onzero pattern matrices
setClass("ndenseMatrix", contains = c("nMatrix", "denseMatrix", "VIRTUAL"),
	 slots = c(x = "logical"),
         validity = function(object) .Call(ndenseMatrix_validate, object))

## Virtual class of dense, logical matrices
setClass("ldenseMatrix", contains = c("lMatrix", "denseMatrix", "VIRTUAL"))

if(FALSE) { # --NOT YET--
## Virtual class of dense, integer matrices
setClass("idenseMatrix", contains = c("iMatrix", "denseMatrix", "VIRTUAL"))
} # --NOT YET--

## Virtual class of dense, double matrices
setClass("ddenseMatrix", contains = c("dMatrix", "denseMatrix", "VIRTUAL"))

if(FALSE) { # --NOT YET--
## Virtual class of dense, complex matrices
setClass("zdenseMatrix", contains = c("zMatrix", "denseMatrix", "VIRTUAL"))
} # --NOT YET--


## ....... Virtual Dense ... class intersections .......................
##                               {for method dispatch}

if(FALSE) {
## This is "natural" but gives WARNINGs when other packages use "it"
setClass("geMatrix", contains = c("denseMatrix", "generalMatrix", "VIRTUAL"))
} else {
## This may work better for other packages
## --> setClassUnion() ... below
}


## ------ Virtual Sparse -----------------------------------------------

## Virtual class of sparse matrices
## * includes diagonal matrices since 0.999375-11 (2008-07)
setClass("sparseMatrix", contains = c("Matrix", "VIRTUAL"))


## ...... Virtual Sparse ... by storage ................................

## Virtual class of sparse matrices in compressed sparse column (CSC) format
setClass("CsparseMatrix", contains = c("sparseMatrix", "VIRTUAL"),
	 slots = c(i = "integer", p = "integer"),
	 prototype = prototype(p = 0L), # to be valid
         validity = function(object) .Call(CsparseMatrix_validate, object))

## Virtual class of sparse matrices in compressed sparse row (CSR) format
setClass("RsparseMatrix", contains = c("sparseMatrix", "VIRTUAL"),
	 slots = c(p = "integer", j = "integer"),
	 prototype = prototype(p = 0L), # to be valid
	 validity = function(object) .Call(RsparseMatrix_validate, object))

## Virtual class of sparse matrices in triplet format
setClass("TsparseMatrix", contains = c("sparseMatrix", "VIRTUAL"),
	 slots = c(i = "integer", j = "integer"),
	 validity = function(object) .Call(TsparseMatrix_validate, object))

## Virtual class of diagonal matrices
setClass("diagonalMatrix", contains = c("sparseMatrix", "VIRTUAL"),
         slots = c(diag = "character"),
	 prototype = prototype(diag = "N"),
         validity = function(object) .Call(diagonalMatrix_validate, object))

if(FALSE) { # --NOT YET--
## These methods would allow initialization of zero matrices _without_ 'p',
## as in the call new("dgCMatrix", Dim = c(6L, 6L)).  However, they would
## also incur a small performance penalty on all other new("..[CR]Matrix")
## calls.
setMethod("initialize", "CsparseMatrix",
          function(.Object, ...) {
              ## Suboptimal if ...names() is NULL or if 'Dim' is missing
              ## but that will "never" happen if ...length() is nonzero:
              if(...length() &&
                 all((nms <- ...names()) != "p") &&
                 length(w <- which(nms == "Dim")) &&
                 !is.character(validDim(d <- ...elt(w[1L]))))
                  callNextMethod(.Object, ..., p = integer(d[2L] + 1))
              else callNextMethod()
          })

setMethod("initialize", "RsparseMatrix",
          function(.Object, ...) {
              ## Suboptimal if ...names() is NULL or if 'Dim' is missing
              ## but that will "never" happen if ...length() is nonzero:
              if(...length() &&
                 all((nms <- ...names()) != "p") &&
                 length(w <- which(nms == "Dim")) &&
                 !is.character(validDim(d <- ...elt(w[1L]))))
                  callNextMethod(.Object, ..., p = integer(d[1L] + 1))
              else callNextMethod()
          })
} # --NOT YET--


## ...... Virtual Sparse ... by kind ...................................

## Virtual class of sparse, _n_onzero pattern matrices
## * these are the "pattern" matrices from "symbolic analysis" of sparse OPs
setClass("nsparseMatrix", contains = c("nMatrix", "sparseMatrix", "VIRTUAL"))

## Virtual class of sparse, logical matrices
setClass("lsparseMatrix", contains = c("lMatrix", "sparseMatrix", "VIRTUAL"))

if(FALSE) { # --NOT YET--
## Virtual class of sparse, integer matrices
setClass("isparseMatrix", contains = c("iMatrix", "sparseMatrix", "VIRTUAL"))
} # --NOT YET--

## Virtual class of sparse, double matrices
setClass("dsparseMatrix", contains = c("dMatrix", "sparseMatrix", "VIRTUAL"))

if(FALSE) { # --NOT YET--
## Virtual class of sparse, complex matrices
setClass("zsparseMatrix", contains = c("zMatrix", "sparseMatrix", "VIRTUAL"))
} # --NOT YET--


## ...... Virtual Sparse ... class intersections .......................
##                               {for method dispatch}

if(FALSE) {
## This is "natural" but gives WARNINGs when other packages use "it"
setClass("nCsparseMatrix",
         contains = c("nsparseMatrix", "CsparseMatrix", "VIRTUAL"))
setClass("lCsparseMatrix",
         contains = c("lsparseMatrix", "CsparseMatrix", "VIRTUAL"))
setClass("iCsparseMatrix",
         contains = c("isparseMatrix", "CsparseMatrix", "VIRTUAL"))
setClass("dCsparseMatrix",
         contains = c("dsparseMatrix", "CsparseMatrix", "VIRTUAL"))
setClass("zCsparseMatrix",
         contains = c("zsparseMatrix", "CsparseMatrix", "VIRTUAL"))
} else {
## These may work better for other packages
## --> setClassUnion() ... below
}


## ====== Non-Virtual Subclasses =======================================

## ------ Non-Virtual Dense --------------------------------------------

## ...... Dense, _n_onzero pattern .....................................

## Unpacked, general
setClass("ngeMatrix",
         contains = c("unpackedMatrix", "ndenseMatrix", "generalMatrix"))

## Unpacked, triangular
setClass("ntrMatrix",
	 contains = c("unpackedMatrix", "ndenseMatrix", "triangularMatrix"))

## Unpacked, symmetric
setClass("nsyMatrix",
         contains = c("unpackedMatrix", "ndenseMatrix", "symmetricMatrix"))

## Packed, triangular
setClass("ntpMatrix",
	 contains = c("packedMatrix", "ndenseMatrix", "triangularMatrix"))

## Packed, symmetric
setClass("nspMatrix",
	 contains = c("packedMatrix", "ndenseMatrix", "symmetricMatrix"))


## ...... Dense, logical ...............................................

## Unpacked, general
setClass("lgeMatrix",
         contains = c("unpackedMatrix", "ldenseMatrix", "generalMatrix"))

## Unpacked, triangular
setClass("ltrMatrix",
	 contains = c("unpackedMatrix", "ldenseMatrix", "triangularMatrix"))

## Unpacked, symmetric
setClass("lsyMatrix",
         contains = c("unpackedMatrix", "ldenseMatrix", "symmetricMatrix"))

## Packed, triangular
setClass("ltpMatrix",
	 contains = c("packedMatrix", "ldenseMatrix", "triangularMatrix"))

## Packed, symmetric
setClass("lspMatrix",
	 contains = c("packedMatrix", "ldenseMatrix", "symmetricMatrix"))


## ...... Dense, double ................................................

## Unpacked, general
setClass("dgeMatrix",
         contains = c("unpackedMatrix", "ddenseMatrix", "generalMatrix"))

## Unpacked, triangular
setClass("dtrMatrix",
	 contains = c("unpackedMatrix", "ddenseMatrix", "triangularMatrix"))

## Unpacked, symmetric
setClass("dsyMatrix",
         contains = c("unpackedMatrix", "ddenseMatrix", "symmetricMatrix"))

## Unpacked, symmetric, positive semidefinite
setClass("dpoMatrix", contains = "dsyMatrix",
	 validity = function(object) .Call(dpoMatrix_validate, object))

## Unpacked, symmetric, positive semidefinite, correlation
setClass("corMatrix", contains = "dpoMatrix",
         slots = c(sd = "numeric"),
	 validity = function(object) .Call(corMatrix_validate, object))

## Packed, triangular
setClass("dtpMatrix",
	 contains = c("packedMatrix", "ddenseMatrix", "triangularMatrix"))

## Packed, symmetric
setClass("dspMatrix",
	 contains = c("packedMatrix", "ddenseMatrix", "symmetricMatrix"))

## Packed, symmetric, positive semidefinite
setClass("dppMatrix", contains = "dspMatrix",
         validity = function(object) .Call(dppMatrix_validate, object))

if(FALSE) { # --NOT YET--
## Packed, symmetric, positive semidefinite, correlation
setClass("copMatrix", contains = "dppMatrix",
         slots = c(sd = "numeric"),
	 validity = function(object) .Call(copMatrix_validate, object))
} # --NOT YET--


## ------ Non-Virtual Sparse -------------------------------------------

## NB: We should _not_ have .t[CRT]Matrix inherit from .g[CRT]Matrix,
##     because a .t[CRT]Matrix could be less than fully stored if diag = "U".
##     Methods for .g[CRT]Matrix applied to such .t[CRT]Matrix" could produce
##     incorrect results, even though all slots are present.

## ...... Sparse, nonzero pattern ......................................

## NB: Unlike [^n]sparseMatrix (below), there is no 'x' slot to validate here.

## CSC, general
setClass("ngCMatrix",
	 contains = c("CsparseMatrix", "nsparseMatrix", "generalMatrix"))

## CSC, triangular
setClass("ntCMatrix",
	 contains = c("CsparseMatrix", "nsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(tCMatrix_validate, object))

## CSC, symmetric
setClass("nsCMatrix",
	 contains = c("CsparseMatrix", "nsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(sCMatrix_validate, object))

## CSR, general
setClass("ngRMatrix",
	 contains = c("RsparseMatrix", "nsparseMatrix", "generalMatrix"))

## CSR, triangular
setClass("ntRMatrix",
	 contains = c("RsparseMatrix", "nsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(tRMatrix_validate, object))

## CSR, symmetric
setClass("nsRMatrix",
	 contains = c("RsparseMatrix", "nsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(sRMatrix_validate, object))

## Triplet general
setClass("ngTMatrix",
	 contains = c("TsparseMatrix", "nsparseMatrix", "generalMatrix"))

## Triplet, triangular
setClass("ntTMatrix",
	 contains = c("TsparseMatrix", "nsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(tTMatrix_validate, object))

## Triplet, symmetric
setClass("nsTMatrix",
	 contains = c("TsparseMatrix", "nsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(sTMatrix_validate, object))


## ...... Sparse, logical ..............................................

## CSC, general
setClass("lgCMatrix",
	 contains = c("CsparseMatrix", "lsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xgCMatrix_validate, object))

## CSC, triangular
setClass("ltCMatrix",
	 contains = c("CsparseMatrix", "lsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(xtCMatrix_validate, object))

## CSC, symmetric
setClass("lsCMatrix",
	 contains = c("CsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
	 validity = function(object) .Call(xsCMatrix_validate, object))

## CSR, general
setClass("lgRMatrix",
	 contains = c("RsparseMatrix", "lsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xgRMatrix_validate, object))

## CSR, triangular
setClass("ltRMatrix",
	 contains = c("RsparseMatrix", "lsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(xtRMatrix_validate, object))

## CSR, symmetric
setClass("lsRMatrix",
	 contains = c("RsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
	 validity = function(object) .Call(xsRMatrix_validate, object))

## Triplet, general
setClass("lgTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xgTMatrix_validate, object))

## Triplet, triangular
setClass("ltTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(xtTMatrix_validate, object))

## Triplet, symmetric
setClass("lsTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
	 validity = function(object) .Call(xsTMatrix_validate, object))

## Diagonal
setClass("ldiMatrix", contains = c("diagonalMatrix", "lMatrix"))


## ...... Sparse, double ...............................................

## CSC, general
setClass("dgCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xgCMatrix_validate, object))

## CSC, triangular
setClass("dtCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(xtCMatrix_validate, object))

## CSC, symmetric
setClass("dsCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsCMatrix_validate, object))

## CSR, general
setClass("dgRMatrix",
	 contains = c("RsparseMatrix", "dsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xgRMatrix_validate, object))

## CSR, triangular
setClass("dtRMatrix",
	 contains = c("RsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(xtRMatrix_validate, object))

## CSR, symmetric
setClass("dsRMatrix",
	 contains = c("RsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsRMatrix_validate, object))

## Triplet, general
setClass("dgTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xgTMatrix_validate, object))

## Triplet, triangular
setClass("dtTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(xtTMatrix_validate, object))

## Triplet, symmetric
setClass("dsTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsTMatrix_validate, object))

## Diagonal
setClass("ddiMatrix", contains = c("diagonalMatrix", "dMatrix"))

if (FALSE) { # TODO
## CSC, symmetic, positive semidefinite
setClass("dpCMatrix", contains = "dsCMatrix",
	 validity = function(object) TODO("test that 'object' is positive semidefinite"))

## Indicator matrix of a factor
setClass("indicator", contains = "dgCMatrix",
         slots = c(levels = "character"),
         validity = function(object) TODO("test that there exists a factor 'g' such that identical(object, as(g, \"sparseMatrix\")) is TRUE"))
} # TODO


## ...... Sparse, index ................................................

## Row index
setClass("indMatrix", contains = c("sparseMatrix", "generalMatrix"),
	 slots = c(perm = "integer"),
         validity = function(object) .Call(indMatrix_validate, object))

## Permutation
setClass("pMatrix", contains = c("indMatrix"),
	 validity = function(object) .Call(pMatrix_validate, object))

if(FALSE) {
## MJ: We really ought to support a _column_ index, too, for which
##     multiplication on the right is equivalent to selecting columns.
##     Is it too late for us to define indMatrix like below ??
setClass("indMatrix", contains = c("sparseMatrix", "generalMatrix"),
         slots = c(margin = "integer", perm = "integer"),
         prototype = prototype(margin = 1L), # to be valid
         ## indMatrix_validate() also checking for 'margin' in 1:2, etc.
	 validity = function(object) .Call(indMatrix_validate, object))
}


########################################################################
##  2. MatrixFactorization
########################################################################

## ------ The Mother Class "MatrixFactorization" -----------------------

setClass("MatrixFactorization", contains = "VIRTUAL",
         slots = c(Dim = "integer"),
         validity = function(object) .Call(MatrixFactorization_validate, object))


## ------ LU -----------------------------------------------------------

setClass("LU", contains = c("MatrixFactorization", "VIRTUAL"))

## MJ: if BunchKaufman should extend dtrMatrix,
##     then maybe denseLU should extend dgeMatrix ...

setClass("denseLU", contains = "LU",
	 slots = c(Dimnames = "list", x = "numeric", perm = "integer"),
	 validity = function(object) .Call(denseLU_validate, object))

setClass("sparseLU", contains = "LU",
	 slots = c(L = "dtCMatrix", U = "dtCMatrix",
		   p = "integer", q = "integer"),
         validity = function(object) .Call(sparseLU_validate, object))

## unused:
if(FALSE) {
setClass("csn_LU",
         slots = c(L = "dgCMatrix", U = "dgCMatrix", Pinv = "integer"))
setClass("css_LU",
         slots = c(Q = "integer", nz = "integer"))
}


## ------ QR -----------------------------------------------------------

if(FALSE) {
## MJ: It would nice to have symmetry with LU, but then we would need
##     to define methods already available for S3 class 'qr'.  Still ...
setClass("QR", contains = c("MatrixFactorization", "VIRTUAL"))

setClass("denseQR", contains = "QR",
         ## based on S3 class 'qr':
         slots = c(x = "numeric", qraux = "numeric",
                   rank = "integer", pivot = "integer",
                   useLAPACK = "logical"),
         validity = function(object) .Call(denseQR_validate, object))
}

setClass("sparseQR", contains = "MatrixFactorization",
	 slots = c(beta = "numeric", V = "dgCMatrix", R = "dgCMatrix",
		   p = "integer", q = "integer"),
	 validity = function(object) .Call(sparseQR_validate, object))

## unused:
if(FALSE) {
setClass("csn_QR",
         slots = c(L = "dgCMatrix", U = "dgCMatrix", beta = "numeric"))
setClass("css_QR",
         slots = c(Q = "integer", Pinv = "integer",
                   nz = "integer", cp = "integer", parent = "integer"))
}


## ------ Cholesky -----------------------------------------------------

setClass("CholeskyFactorization",
         contains = c("MatrixFactorization", "VIRTUAL"))


## ...... Dense ........................................................

setClass("Cholesky",  contains = c("dtrMatrix", "CholeskyFactorization"),
         validity = function(object) .Call(Cholesky_validate, object))

setClass("pCholesky", contains = c("dtpMatrix", "CholeskyFactorization"),
         validity = function(object) .Call(pCholesky_validate, object))

## unused:
if(FALSE) {
setClass("LDL", contains = c("dtrMatrix", "CholeskyFactorization"))
}


## ...... Sparse .......................................................

## S4 representation of C struct 'cholmod_factor'
setClass("CHMfactor", contains = c("CholeskyFactorization", "VIRTUAL"),
	 slots = c(colcount = "integer", perm = "integer", type = "integer"),
	 validity = function(object) .Call(CHMfactor_validate, object))

## Supernodal CHMfactor
setClass("CHMsuper", contains = c("CHMfactor", "VIRTUAL"),
	 slots = c(pi = "integer", px = "integer", s = "integer",
                   super = "integer"),
	 validity = function(object) .Call(CHMsuper_validate, object))

setClass("nCHMsuper", contains = "CHMsuper") # symbolic factorization
setClass("dCHMsuper", contains = "CHMsuper", slots = c(x = "numeric"))

## Simplicial CHMfactor
setClass("CHMsimpl",
	 contains = c("CHMfactor", "VIRTUAL"),
	 slots = c(p = "integer", i = "integer", nz = "integer",
		   prv = "integer", nxt = "integer"),
	 validity = function(object) .Call(CHMsimpl_validate, object))

setClass("nCHMsimpl", contains = "CHMsimpl") # symbolic factorization
setClass("dCHMsimpl", contains = "CHMsimpl", slots = c(x = "numeric"))


## ------ Bunch-Kaufman ------------------------------------------------

setClass("BunchKaufman", contains = c("dtrMatrix", "MatrixFactorization"),
	 slots = c(perm = "integer"),
	 validity = function(object) .Call(BunchKaufman_validate, object))

setClass("pBunchKaufman", contains = c("dtpMatrix", "MatrixFactorization"),
	 slots = c(perm = "integer"),
	 validity = function(object) .Call(pBunchKaufman_validate, object))


## ------ Schur --------------------------------------------------------

## MJ: My preference would be for signature EValues = "complex" here,
##     _even if_ base::eigen() gives a "numeric" _or_ "complex" result.
##     We can already think of type "complex" as containing types "double"
##     and "integer".  And for backwards compatibility, we could always
##     define a method for initialize() coercing 'EValues' from numeric
##     to complex.

## For eigenvalues:
setClassUnion("number", members = c("numeric", "complex"))

setClass("Schur", contains = "MatrixFactorization",
         slots = c(Q = "Matrix", T = "Matrix", EValues = "number"),
         validity = function(object) .Call(Schur_validate, object))


########################################################################
##  3. sparseVector
########################################################################

## ------ The Mother Class 'sparseVector' ------------------------------

## "longindex" should allow sparseVector of length >= 2^31,
## which is necessary, e.g., when coercing from large sparseMatrix
##
## > setClass("longindex", contains = "numeric")
##
## but we use "numeric" here instead (for simplicity? efficiency?) ...
## note that "numeric" contains "integer" (whether I like it or not) ...

setClass("sparseVector", contains = "VIRTUAL",
         slots = c(length = "numeric", i = "numeric"), # 1-based index!
         prototype = prototype(length = 0),
         validity = function(object) {
             len <- object@length
             if(length(len) != 1L)
                 return("'length' slot does not have length 1")
             if(!is.finite(len))
                 return("'length' slot is not finite")
             if(len < 0)
                 return("'length' slot is negative")
             i <- object@i
             i.len <- length(i)
             if(i.len == 0L)
                 return(TRUE)
             if(i.len > len)
                 return("'i' slot has length greater than 'length' slot")
             i.num <- is.double(i)
             if(i.num)
                 i <- trunc(i)
             i.uns <- is.unsorted(i, strictly = TRUE)
             if(is.na(i.uns))
                 "'i' slot contains NA"
             else if(i.uns || i[1L] < 1 || i[i.len] > len) {
                 m <- if(i.uns)
                          "'i' slot is not strictly increasing"
                      else "'i' slot has elements not in 1:<'length' slot>"
                 if(i.num)
                     paste0(m, " after truncation towards zero")
                 else m
             } else TRUE
         })

## Allow users to do new("[nlidz]sparseVector", i=, x=) with unsorted 'i'
setMethod("initialize", "sparseVector",
          function(.Object, i, x, ...) {
              if(has.x <- !missing(x))
                  x <- x # MJ: why is this necessary?
              if(!missing(i)) {
                  i.uns <- is.unsorted(i, strictly = TRUE)
                  i <-
                      if(is.na(i.uns) || !i.uns)
                          i
                      else {
                          ## we know that there are no NA, and the order of
                          ## ties does not matter (since ties are an error),
                          ## hence it is safe to use "quick" here
                          m <- if(is.integer(length(i))) "radix" else "quick"
                          if(.hasSlot(.Object, "x") && has.x) {
                              s <- sort.int(i, method = m, index.return = TRUE)
                              x <- x[s$ix]
                              s$x
                          } else sort.int(i, method = m)
                      }
              }
              callNextMethod()
          })


## ------ Non-Virtual Subclasses ---------------------------------------

.valid.xsparseVector <- function(object) {
    if(length(object@x) != length(object@i))
        "'i' and 'x' slots do not have equal length"
    else TRUE
}

## No 'x' slot, hence nothing more to validate:
setClass("nsparseVector", contains = "sparseVector")

setClass("lsparseVector", contains = "sparseVector",
	 slots = c(x = "logical"),
	 validity = .valid.xsparseVector)

setClass("isparseVector", contains = "sparseVector",
	 slots = c(x = "integer"),
	 validity = .valid.xsparseVector)

setClass("dsparseVector", contains = "sparseVector",
	 slots = c(x = "numeric"),
	 validity = .valid.xsparseVector)

setClass("zsparseVector", contains = "sparseVector",
	 slots = c(x = "complex"),
	 validity = .valid.xsparseVector)

rm(.valid.xsparseVector)


########################################################################
##  4. Index and more "miscellaneous" classes, but _not_ class unions
########################################################################

## Idea: represent x = c(seq(from1, to1, by1), seq(from2, to2, by2), ...)
##       as list(first = x[1L], rle = rle(diff(x)))
setClass("rleDiff",
         ## MJ: simpler would be slots = c(first=, lengths=, values=) ...
         slots = c(first = "numeric", rle = "rle"),
	 prototype = prototype(first = integer(0L), rle = rle(integer(0L))),
	 validity = function(object) {
	     if(length(object@first) != 1L)
		 "'first' slot does not have length 1"
	     else if(!is.list(rle <- object@rle))
                 "'rle' slot is not a list"
             else if(length(rle) != 2L)
                 "'rle' slot does not have length 2"
             else if(is.null(nms <- names(rle)) ||
                     anyNA(match(nms, c("lengths", "values"))))
                 "'rle' slot does not have names \"lengths\", \"values\""
             else if(!is.numeric(lens <- rle$lengths))
                 "'lengths' is not numeric"
             else if(!is.numeric(vals <- rle$values))
                 "'values' is not numeric"
             else if(length(lens) != length(vals))
                 "'lengths' and 'values' do not have equal length"
             else if(length(lens) == 0L)
                 TRUE
             else if(anyNA(lens))
                 "'lengths' contains NA"
             else if(is.double(lens)) {
                 if(!(all(is.finite(r <- range(lens))) &&
                      all(lens == trunc(lens))))
                     "'lengths' is not integer-valued"
                 else if(r[1L] < 1)
                     "'lengths' is not positive"
                 else TRUE
             } else {
                 if(min(lens) < 1L)
                     "'lengths' is not positive"
                 else TRUE
             }
	 })

## Idea: represent x = c(seq(from1, to1, by1), seq(from2, to2, by2), ...)
##       as rbind(c(from1, from2, ...), c(to1, to2, ...), c(by1, by2, ...))
## MM: (2010-03-04) more efficient than "rleDiff" [TODO: write rleDiff<->seqMat]
## MJ: (2022-09-06) data.frame(from, to, by) could be _handled_ more efficiently
setClass("seqMat", contains = "matrix",
	 prototype = prototype(matrix(integer(0L), nrow = 3L, ncol = 0L)),
	 validity = function(object) {
             if(!is.numeric(object))
                 "matrix is not numeric"
             else if(nrow(object) != 3L)
		 "matrix does not have 3 rows"
             else if(anyNA(object))
                 "matrix contains NA"
             else if(is.double(object) && !(all(is.finite(range(object))) &&
                                            all(object == trunc(object))))
                 "matrix is not integer-valued"
             else {
                 from <- object[1L, ]
                 to   <- object[2L, ]
                 by   <- object[3L, ]
                 if(any((from < to & by <= 0) | (from > to & by >= 0)))
                     "degenerate sequence(s): sign(to - from) != sign(by)"
                 else TRUE
             }
	 })

## Idea: _ab_stract index
## MJ: (2022-09-06) why not just
##     setClassUnion("abIndex", members = c("numeric", "rleDiff", "seqMat")) ?
setClass("abIndex",
         slots = c(kind = "character", x = "numeric", rleD = "rleDiff"),
         prototype = prototype(kind = "int32", x = integer(0L)),
         validity = function(object) {
             ## MJ: should 'rleD' be "empty" if kind != "rleDiff" ?
             if(length(kind <- object@kind) != 1L)
                 "'kind' slot does not have length 1"
             else switch(kind,
                         "int32" =
                             if(is.integer(object@x))
                                 TRUE
                             else "kind=\"int32\" but 'x' slot is not of type \"integer\"",
                         "double" =
                             if(is.double(object@x))
                                 TRUE
                             else "kind=\"double\" but 'x' slot is not of type \"double\"",
                         "rleDiff" =
                             if(length(object@x) == 0L)
                                 TRUE
                             else "kind=\"rleDiff\" but 'x' slot is nonempty",
                         ## otherwise:
                         "'kind' is not \"int32\", \"double\", or \"rleDiff\"")
         })

setClass("determinant",
         ## based on S3 class 'det':
         slots = c(modulus = "numeric", logarithm = "logical",
                   sign = "integer", call = "call"),
         validity = function(object) {
             if(length(logarithm <- object@logarithm) != 1L)
                 "'logarithm' slot does not have length 1"
             else if(is.na(logarithm))
                 "'logarithm' is not TRUE or FALSE"
             else if(length(modulus <- object@modulus) != 1L)
                 "'modulus' slot does not have length 1"
             else if(logarithm && !is.na(modulus) && modulus < 0)
                 "logarithm=FALSE but 'modulus' slot is negative"
             else if(length(sign <- object@sign))
                 "'sign' slot does not have length 1"
             else if(is.na(sign) || (sign != -1L && sign != 1L))
                 "'sign' slot is not -1 or 1"
             else TRUE
         })

## unused:
if(FALSE) {
setClass("logic", contains = "raw") # "raw" rather than "logical"
}


########################################################################
##  5. Class unions
########################################################################

## NB: numeric = { double, integer }
## NB: many of these are _not_ exported, on purpose

## Union of matrix and Matrix:
## * for certain "catch-all" methods; see, e.g., ./products.R
## * note that is(x, "mMatrix") is stricter than length(dim(x)) == 2L,
##   which allows, e.g., class 'table'
setClassUnion("mMatrix",
              members = c("matrix", "Matrix"))

if(FALSE) { # --NOT YET--
## for setMethod("c", "numMatrixLike"), once that works
setClassUnion("numMatrixLike",
              members = c("logical", "numeric", "mMatrix"))
} # --NOT YET--

if(TRUE) {
## MJ: Somewhat surprisingly, these are not actually used anywhere;
##     xsparseVector is only _mentioned_ in ../man/sparseVector-class.Rd.
##     Keeping for now, if only for didactic reasons ...

## Subclasses of Matrix with an 'x' slot:
## NB: the 'x' slot need not contain all of the data (e.g., when diag = "U")
setClassUnion("xMatrix",
              members = c("ndenseMatrix", "lMatrix", "iMatrix",
                          "dMatrix", "zMatrix"))

## Subclasses of sparseVector with an 'x' slot:
setClassUnion("xsparseVector",
              members = c("lsparseVector", "isparseVector",
                          "dsparseVector", "zsparseVector"))
}

## Intersection of denseMatrix and generalMatrix:
## * currently only used in ./diagMatrix.R
setClassUnion("geMatrix",
              members = c("ngeMatrix", "lgeMatrix", "dgeMatrix"))

## Intersection of nsparseMatrix and CsparseMatrix:
## * _should_ be closer to its members than nsparseMatrix and CsparseMatrix
##   but it is _not_
## * a "fix" would be to define nCsparseMatrix as a (non-union) virtual class
##   _and_ have n[gts]CMatrix extend it
setClassUnion("nCsparseMatrix",
              members = c("ngCMatrix", "ntCMatrix", "nsCMatrix"))
setClassUnion("lCsparseMatrix",
              members = c("lgCMatrix", "ltCMatrix", "lsCMatrix"))
setClassUnion("dCsparseMatrix",
              members = c("dgCMatrix", "dtCMatrix", "dsCMatrix"))

if(FALSE) { # --NOT YET--
## CHOLMOD-like sparseMatrix, i.e., excluding diagonalMatrix and indMatrix:
## * would be useful, e.g., in ./products.R for '%&%',
##   but at the moment it affects dispatch too much
setClassUnion("CRTsparseMatrix",
              members = c("CsparseMatrix", "RsparseMatrix", "TsparseMatrix"))
} # --NOT YET--

## Atomic vectors:
## * note that is(<atomic matrix>, "atomicVector") is FALSE
##   even though is.atomic(<atomic matrix>) is TRUE
setClassUnion("atomicVector",
              members = c("logical", "numeric", "complex", "raw", "character"))

## Numeric-like vectors:
## * for methods handling logical and integer as double; see, e.g., ./solve.R
setClassUnion("numLike",
              members = c("logical", "numeric"))

## Index vectors:
## * for 'i' in x[i], x[i, ], x[, i], etc.
## * TODO: include rleDiff
setClassUnion("index",
              members = c("logical", "numeric", "character"))

## Subassignment values:
## * for 'value' in x[i, j] <- value
setClassUnion("replValue",
              members = c("logical", "numeric", "complex", "raw"))
setClassUnion("replValueSp",
              ## MJ: why Matrix but not matrix ??
              members = c("replValue", "sparseVector", "Matrix"))
