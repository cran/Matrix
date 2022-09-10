## README:
##
## Validity methods should assume that methods for superclasses have passed,
## following 'validObject()'. We should _not_ be testing, e.g., length(Dim),
## typeof(Dimnames), etc. repeatedly ...
##
## When checking whether a class is validated correctly, take care to follow
## the 'contains' recursively!!


## To be used in 'initialize' method for "Matrix", or other constructors
## NB: This must be defined _here_ and _not_ be migrated to ./Auxiliaries.R
fixupDN <- function(dn) .Call(R_DimNames_fixup, dn)

## MJ: no longer
if (FALSE) {
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
}


## ====== Virtual Classes ==============================================

## ------ The Mother Class "Matrix" ------------------------------------

## Virtual class of all Matrix objects
setClass("Matrix", contains = "VIRTUAL",
	 slots = c(Dim = "integer", Dimnames = "list"),
	 prototype = prototype(Dim = integer(2L), Dimnames = list(NULL, NULL)),
	 validity = function(object) .Call(Matrix_validate, object))

if(FALSE) {
## Allowing 'Dimnames' to define 'Dim' ... would require changes
## to DimNames_validate() in ../src/Mutils.c and how it is used
## in Matrix_validate()
setMethod("initialize", "Matrix",
          function(.Object, ...) {
              .Object <- callNextMethod()
              if(length(args <- list(...)) &&
                  any(nzchar(anames <- names(args))) &&
                  "Dimnames" %in% anames) {
                  ## Coerce non-character 'Dimnames[[i]]' to character
                  ## and set zero-length 'Dimnames[[i]]' to NULL
                  .Object@Dimnames <- DN <- fixupDN(.Object@Dimnames)
                  if(!("Dim" %in% anames ||
                       is.null(DN[[1L]]) || is.null(DN[[2L]]))) {
                      ## Take 'Dim' from lengths of 'Dimnames'
                      .Object@Dim <- lengths(DN, use.names = FALSE)
                  }
              }
              .Object
          })
}

setMethod("initialize", "Matrix",
          function(.Object, ...) {
              .Object <- callNextMethod()
              if(length(args <- list(...)) &&
                  any(nzchar(anames <- names(args))) &&
                  "Dimnames" %in% anames) {
                  ## Coerce non-character 'Dimnames[[i]]' to character
                  ## and set zero-length 'Dimnames[[i]]' to NULL
                  .Object@Dimnames <- fixupDN(.Object@Dimnames)
              }
              .Object
          })


## ------ Virtual by structure -----------------------------------------

## Virtual class of composite matrices,
## i.e., those for which it makes sense to define a factorization
setClass("compMatrix", contains = c("Matrix", "VIRTUAL"),
	 slots = c(factors = "list"),
         validity = function(object) .Call(compMatrix_validate, object))

## Virtual class of matrices that are not symmetric, triangular, _or diagonal_
setClass("generalMatrix", contains = c("compMatrix", "VIRTUAL"))

## Virtual class of symmetric matrices
setClass("symmetricMatrix", contains = c("compMatrix", "VIRTUAL"),
	 slots = c(uplo = "character"),
	 prototype = prototype(uplo = "U"),
	 validity = function(object) .Call(symmetricMatrix_validate, object))

## Virtual class of triangular matrices
setClass("triangularMatrix", contains = c("Matrix", "VIRTUAL"),
	 slots = c(uplo = "character", diag = "character"),
	 prototype = prototype(uplo = "U", diag = "N"),
	 validity = function(object) .Call(triangularMatrix_validate, object))


## ------ Virtual by kind ----------------------------------------------

## Virtual class of double matrices
setClass("dMatrix", contains = c("Matrix", "VIRTUAL"), slots = c(x = "numeric"),
         validity = function(object) .Call(dMatrix_validate, object))

## Virtual class of logical matrices,
## >>> typically the result of comparisons, e.g., <dMatrix> <relop> <dMatrix>
## >>> hence NA are allowed and distinct from TRUE, in contrast with "nMatrix"
setClass("lMatrix", contains = c("Matrix", "VIRTUAL"), slots = c(x = "logical"),
         validity = function(object) .Call(lMatrix_validate, object))

## Virtual class of nonzero pattern (or simply "pattern") matrices
## NB: only "ndenseMatrix" requires an 'x' slot
setClass("nMatrix", contains = c("Matrix", "VIRTUAL"))

## Virtual class of integer matrices
setClass("iMatrix", contains = c("Matrix", "VIRTUAL"), slots = c(x = "integer"),
         validity = function(object) .Call(iMatrix_validate, object))

## Virtual class of complex matrices
## >>> 'z' as in the names of Lapack routines
setClass("zMatrix", contains = c("Matrix", "VIRTUAL"), slots = c(x = "complex"),
         validity = function(object) .Call(zMatrix_validate, object))


## ------ Virtual Dense ------------------------------------------------

## Virtual class of dense matrices
## NB: included diagonal matrices until 0.999375-11 (2008-07)
## NB: includes "unpacked" _and_ "packed" matrices
setClass("denseMatrix", contains = c("Matrix", "VIRTUAL"))


## ...... Virtual Dense ... by storage .................................

## Virtual class of dense, "unpacked" matrices, s.t. length(.@x) == n*n
setClass("unpackedMatrix", contains = c("denseMatrix", "VIRTUAL"),
         validity = function(object) .Call(unpackedMatrix_validate, object))

## Virtual class of dense, "packed" matrices, s.t. length(.@x) == n*(n+1)/2
setClass("packedMatrix", contains = c("denseMatrix", "VIRTUAL"),
         slots = c(uplo = "character"),
         prototype = prototype(uplo = "U"),
	 validity = function(object) .Call(packedMatrix_validate, object))


## ...... Virtual Dense ... by kind ....................................

## Virtual class of dense, double matrices
setClass("ddenseMatrix", contains = c("dMatrix", "denseMatrix", "VIRTUAL"))

## Virtual class of dense, logical matrices
setClass("ldenseMatrix", contains = c("lMatrix", "denseMatrix", "VIRTUAL"))

## Virtual class of dense, nonzero pattern matrices
setClass("ndenseMatrix", contains = c("nMatrix", "denseMatrix", "VIRTUAL"),
	 slots = c(x = "logical"),
         validity = function(object) .Call(ndenseMatrix_validate, object))

if(FALSE) { # --NOT YET--
setClass("idenseMatrix", contains = c("iMatrix", "denseMatrix", "VIRTUAL"))
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
## NB: includes diagonal matrices since 0.999375-11 (2008-07)
setClass("sparseMatrix", contains = c("Matrix", "VIRTUAL"))


## ...... Virtual Sparse ... by storage ................................

## Virtual class of diagonal matrices
setClass("diagonalMatrix", contains = c("sparseMatrix", "VIRTUAL"),
         slots = c(diag = "character"),
	 prototype = prototype(diag = "N"),
         validity = function(object) .Call(diagonalMatrix_validate, object))

## Virtual class of sparse matrices with triplet representation
setClass("TsparseMatrix", contains = c("sparseMatrix", "VIRTUAL"),
	 slots = c(i = "integer", j = "integer"),
	 validity = function(object) .Call(Tsparse_validate, object))

## Virtual class of compressed sparse column-oriented matrices
setClass("CsparseMatrix", contains = c("sparseMatrix", "VIRTUAL"),
	 slots = c(i = "integer", p = "integer"),
	 prototype = prototype(p = 0L), # to be valid
         validity = function(object) .Call(Csparse_validate, object))

## Virtual class of compressed sparse row-oriented matrices
setClass("RsparseMatrix", contains = c("sparseMatrix", "VIRTUAL"),
	 slots = c(p = "integer", j = "integer"),
	 prototype = prototype(p = 0L), # to be valid
	 validity = function(object) .Call(Rsparse_validate, object))

if(FALSE) { # --NOT YET-- (fails)
## Would be nice, in theory, for new("dgCMatrix", Dim = c(3L, 3L))
setMethod("initialize", "CsparseMatrix",
          function(.Object, ...) {
              .Object <- callNextMethod()
              .Object@p <- integer(.Object@Dim[2L] + 1L)
              .Object
          })
setMethod("initialize", "RsparseMatrix",
          function(.Object, ...) {
              .Object <- callNextMethod()
              .Object@p <- integer(.Object@Dim[1L] + 1L)
              .Object
          })
} # --NOT YET--

## ...... Virtual Sparse ... by kind ...................................

## Virtual class of sparse, double matrices
setClass("dsparseMatrix", contains = c("dMatrix", "sparseMatrix", "VIRTUAL"))

## Virtual class of sparse, logical matrices
setClass("lsparseMatrix", contains = c("lMatrix", "sparseMatrix", "VIRTUAL"))

## Virtual class of sparse, nonzero pattern matrices
## >>> these are the "pattern" matrices from "symbolic analysis" of sparse OPs
setClass("nsparseMatrix", contains = c("nMatrix", "sparseMatrix", "VIRTUAL"))

if(FALSE) { # --NOT YET--
setClass("isparseMatrix", contains = c("iMatrix", "sparseMatrix", "VIRTUAL"))
setClass("zsparseMatrix", contains = c("zMatrix", "sparseMatrix", "VIRTUAL"))
} # --NOT YET--

## ...... Virtual Sparse ... class intersections .......................
##                               {for method dispatch}

if(FALSE) {
## This is "natural" but gives WARNINGs when other packages use "it"
setClass("dCsparseMatrix",
         contains = c("dsparseMatrix", "CsparseMatrix", "VIRTUAL"))
setClass("lCsparseMatrix",
         contains = c("lsparseMatrix", "CsparseMatrix", "VIRTUAL"))
setClass("nCsparseMatrix",
         contains = c("nsparseMatrix", "CsparseMatrix", "VIRTUAL"))
} else {
## These may work better for other packages
## --> setClassUnion() ... below
}


## ====== Proper (Non-Virtual) Classes =================================

## ------ Proper (Non-Virtual) Dense -----------------------------------

## ...... Dense, double ................................................

## General
## NB: always "unpacked"
setClass("dgeMatrix",
         contains = c("unpackedMatrix", "ddenseMatrix", "generalMatrix"))

## Unpacked, triangular
setClass("dtrMatrix",
	 contains = c("unpackedMatrix", "ddenseMatrix", "triangularMatrix"))

## Packed, triangular
setClass("dtpMatrix",
	 contains = c("packedMatrix", "ddenseMatrix", "triangularMatrix"))

## Unpacked, symmetric
setClass("dsyMatrix",
         contains = c("unpackedMatrix", "ddenseMatrix", "symmetricMatrix"))

## Packed, symmetric
setClass("dspMatrix",
	 contains = c("packedMatrix", "ddenseMatrix", "symmetricMatrix"))

## Unpacked, symmetric, positive semidefinite
setClass("dpoMatrix", contains = "dsyMatrix",
	 validity = function(object) .Call(dpoMatrix_validate, object))

## Packed, symmetric, positive semidefinite
setClass("dppMatrix", contains = "dspMatrix",
         validity = function(object) .Call(dppMatrix_validate, object))

## Unpacked correlation matrices
setClass("corMatrix", contains = "dpoMatrix", slots = c(sd = "numeric"),
	 validity = function(object) .Call(corMatrix_validate, object))


## ...... Dense, logical ...............................................

## General
## NB: always "unpacked"
setClass("lgeMatrix",
         contains = c("unpackedMatrix", "ldenseMatrix", "generalMatrix"))

## Unpacked, triangular
setClass("ltrMatrix",
	 contains = c("unpackedMatrix", "ldenseMatrix", "triangularMatrix"))

## Packed, triangular
setClass("ltpMatrix",
	 contains = c("packedMatrix", "ldenseMatrix", "triangularMatrix"))

## Unpacked, symmetric
setClass("lsyMatrix",
         contains = c("unpackedMatrix", "ldenseMatrix", "symmetricMatrix"))

## Packed, symmetric
setClass("lspMatrix",
	 contains = c("packedMatrix", "ldenseMatrix", "symmetricMatrix"))


## ...... Dense, nonzero pattern .......................................

## General
## NB: always "unpacked"
setClass("ngeMatrix",
         contains = c("unpackedMatrix", "ndenseMatrix", "generalMatrix"))

## Unpacked, triangular
setClass("ntrMatrix",
	 contains = c("unpackedMatrix", "ndenseMatrix", "triangularMatrix"))

## Packed, triangular
setClass("ntpMatrix",
	 contains = c("packedMatrix", "ndenseMatrix", "triangularMatrix"))

## Unpacked, symmetric
setClass("nsyMatrix",
         contains = c("unpackedMatrix", "ndenseMatrix", "symmetricMatrix"))

## Packed, symmetric
setClass("nspMatrix",
	 contains = c("packedMatrix", "ndenseMatrix", "symmetricMatrix"))


## ------ Proper (Non-Virtual) Sparse ----------------------------------

## ...... Sparse, double ...............................................

## Diagonal
setClass("ddiMatrix", contains = c("diagonalMatrix", "dMatrix"))

## Triplet, general
setClass("dgTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xTMatrix_validate, object))

## Triplet, triangular
setClass("dtTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(tTMatrix_validate, object))

## NB: We should _not_ have ".tTMatrix" inherit from ".gTMatrix",
## because a ".tTMatrix" could be less than fully stored if diag = "U".
## Methods for ".gTMatrix" applied to such ".tTMatrix" would produce
## incorrect results, even though all slots are present.

## Triplet, symmetric
setClass("dsTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
         ## Likewise only storing one triangle
	 validity = function(object) .Call(tTMatrix_validate, object))

## CSC, general
setClass("dgCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xCMatrix_validate, object))

## CSC, triangular
setClass("dtCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(tCMatrix_validate, object))

## CSC, symmetric
setClass("dsCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
         ## Likewise only storing one triangle
	 validity = function(object) .Call(tCMatrix_validate, object))

## CSR, general
setClass("dgRMatrix",
	 contains = c("RsparseMatrix", "dsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xRMatrix_validate, object))

## CSR, triangular
setClass("dtRMatrix",
	 contains = c("RsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(tRMatrix_validate, object))

## CSR, symmetric
setClass("dsRMatrix",
	 contains = c("RsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
         ## Likewise only storing one triangle
	 validity = function(object) .Call(tRMatrix_validate, object))

if (FALSE) { ## TODO
## Indicator matrix of a factor ... needs more careful definition
setClass("indicators", contains = "dgCMatrix", slots = c(levels = "character"))

## CSC, symmetic, positive semidefinite
setClass("dpCMatrix", contains = "dsCMatrix",
	 validity = function(object) TODO("test for positive semidefinite ??"))
}


## ...... Sparse, logical ..............................................

## Diagonal
setClass("ldiMatrix", contains = c("diagonalMatrix", "lMatrix"))

## Triplet, general
setClass("lgTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xTMatrix_validate, object))

## Triplet, triangular
setClass("ltTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(tTMatrix_validate, object))

## Triplet, symmetric
setClass("lsTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
         ## Likewise only storing one triangle
	 validity = function(object) .Call(tTMatrix_validate, object))

## CSC, general
setClass("lgCMatrix",
	 contains = c("CsparseMatrix", "lsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xCMatrix_validate, object))

## CSC, triangular
setClass("ltCMatrix",
	 contains = c("CsparseMatrix", "lsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(xCMatrix_validate, object))

## CSC, symmetric
setClass("lsCMatrix",
	 contains = c("CsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
         ## Likewise only storing one triangle
	 validity = function(object) .Call(xCMatrix_validate, object))

## CSR, general
setClass("lgRMatrix",
	 contains = c("RsparseMatrix", "lsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xRMatrix_validate, object))

## CSR, triangular
setClass("ltRMatrix",
	 contains = c("RsparseMatrix", "lsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(tRMatrix_validate, object))

## CSR, symmetric
setClass("lsRMatrix",
	 contains = c("RsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
         ## Likewise only storing one triangle
	 validity = function(object) .Call(tRMatrix_validate, object))


## ...... Sparse, nonzero pattern ......................................

## NB: Unlike [^n]sparseMatrix, there is no 'x' slot to validate here.

## Triplet general
setClass("ngTMatrix",
	 contains = c("TsparseMatrix", "nsparseMatrix", "generalMatrix"))

## Triplet, triangular
setClass("ntTMatrix",
	 contains = c("TsparseMatrix", "nsparseMatrix", "triangularMatrix"))

## Triplet, symmetric
setClass("nsTMatrix",
	 contains = c("TsparseMatrix", "nsparseMatrix", "symmetricMatrix"))

## CSC, general
setClass("ngCMatrix",
	 contains = c("CsparseMatrix", "nsparseMatrix", "generalMatrix"))

## CSC, triangular
setClass("ntCMatrix",
	 contains = c("CsparseMatrix", "nsparseMatrix", "triangularMatrix"))

## CSC, symmetric
setClass("nsCMatrix",
	 contains = c("CsparseMatrix", "nsparseMatrix", "symmetricMatrix"))

## CSR, general
setClass("ngRMatrix",
	 contains = c("RsparseMatrix", "nsparseMatrix", "generalMatrix"))

## CSR, triangular
setClass("ntRMatrix",
	 contains = c("RsparseMatrix", "nsparseMatrix", "triangularMatrix"))

## CSR, symmetric
setClass("nsRMatrix",
	 contains = c("RsparseMatrix", "nsparseMatrix", "symmetricMatrix"))


if(FALSE) { # --NOT YET--

## ...... Sparse, integer ..............................................

## Triplet, general
setClass("igTMatrix",
	 contains = c("TsparseMatrix", "isparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xTMatrix_validate, object))

## Triplet, triangular
setClass("itTMatrix",
	 contains = c("TsparseMatrix", "isparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(tTMatrix_validate, object))

## Triplet, symmetric
setClass("isTMatrix",
	 contains = c("TsparseMatrix", "isparseMatrix", "symmetricMatrix"),
         ## Likewise only storing one triangle
	 validity = function(object) .Call(tTMatrix_validate, object))

## CSC, general
setClass("igCMatrix",
	 contains = c("CsparseMatrix", "isparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xCMatrix_validate, object))

## CSC, triangular
setClass("itCMatrix",
	 contains = c("CsparseMatrix", "isparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(tCMatrix_validate, object))

## CSC, symmetric
setClass("isCMatrix",
	 contains = c("CsparseMatrix", "isparseMatrix", "symmetricMatrix"),
         ## Likewise only storing one triangle
	 validity = function(object) .Call(tCMatrix_validate, object))

## CSR, general
setClass("igRMatrix",
	 contains = c("RsparseMatrix", "isparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xRMatrix_validate, object))

## CSR, triangular
setClass("itRMatrix",
	 contains = c("RsparseMatrix", "isparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(tRMatrix_validate, object))

## CSR, symmetric
setClass("isRMatrix",
	 contains = c("RsparseMatrix", "isparseMatrix", "symmetricMatrix"),
         ## Likewise only storing one triangle
	 validity = function(object) .Call(tRMatrix_validate, object))

} # --NOT YET--



##-------------------- index and permutation matrices--------------------------

setClass("indMatrix", slots = c(perm = "integer"),
	 contains = c("sparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(indMatrix_validate, object))

setClass("pMatrix", slots = c(perm = "integer"),
	 contains = c("indMatrix"),
	 validity = function(object) .Call(pMatrix_validate, object))


### Factorization classes ---------------------------------------------

## Mother class:
setClass("MatrixFactorization", slots = c(Dim = "integer"), contains = "VIRTUAL",
	 validity = function(object) .Call(MatrixFactorization_validate, object))

setClass("CholeskyFactorization", contains = "MatrixFactorization", "VIRTUAL")

## -- Those (exceptions) inheriting from "Matrix" : ---

setClass("Cholesky",  contains = c("dtrMatrix", "CholeskyFactorization"))

#unUsed: setClass("LDL", contains = c("dtrMatrix", "CholeskyFactorization"))

setClass("pCholesky", contains = c("dtpMatrix", "CholeskyFactorization"))

## These are currently only produced implicitly from *solve()
setClass("BunchKaufman", contains = c("dtrMatrix", "MatrixFactorization"),
	 slots = c(perm = "integer"),
	 validity = function(object) .Call(BunchKaufman_validate, object))

setClass("pBunchKaufman", contains = c("dtpMatrix", "MatrixFactorization"),
	 slots = c(perm = "integer"),
	 validity = function(object) .Call(pBunchKaufman_validate, object))

## -- the usual ``non-Matrix'' factorizations : ---------

setClass("CHMfactor", # cholmod_factor struct as S4 object
	 contains = c("CholeskyFactorization", "VIRTUAL"),
	 slots = c(colcount = "integer", perm = "integer", type = "integer"),
	 validity = function(object) .Call(CHMfactor_validate, object))

setClass("CHMsuper",		       # supernodal cholmod_factor
	 contains = c("CHMfactor", "VIRTUAL"),
	 slots = c(super = "integer", pi = "integer", px = "integer",
		   s = "integer"),
	 validity = function(object) .Call(CHMsuper_validate, object))

setClass("CHMsimpl",		       # simplicial cholmod_factor
	 contains = c("CHMfactor", "VIRTUAL"),
	 slots = c(p = "integer", i = "integer", nz = "integer",
		   nxt = "integer", prv = "integer"),
	 validity = function(object) .Call(CHMsimpl_validate, object))

setClass("dCHMsuper", contains = "CHMsuper", slots = c(x = "numeric"))

setClass("nCHMsuper", contains = "CHMsuper")

setClass("dCHMsimpl", contains = "CHMsimpl", slots = c(x = "numeric"))

setClass("nCHMsimpl", contains = "CHMsimpl")

##--- LU ---

setClass("LU", contains = c("MatrixFactorization", "VIRTUAL"))

setClass("denseLU", contains = "LU",
	 slots = c(x = "numeric", perm = "integer", Dimnames = "list"),
	 validity = function(object) .Call(LU_validate, object))

setClass("sparseLU", contains = "LU",
	 slots = c(L = "dtCMatrix", U = "dtCMatrix",
		   p = "integer", q = "integer"))

##--- QR ---

setClass("sparseQR", contains = "MatrixFactorization",
	 slots = c(V = "dgCMatrix", beta = "numeric",
		   p = "integer", R = "dgCMatrix", q = "integer"),
	 validity = function(object) .Call(sparseQR_validate, object))

##-- "SPQR" ---> ./spqr.R  for now

## "denseQR" -- ?  (``a version of''  S3 class "qr")

if (FALSE) { ## unused classes
setClass("csn_QR", slots = c(U = "dgCMatrix", L = "dgCMatrix",
                             beta = "numeric"))

setClass("csn_LU", slots = c(U = "dgCMatrix", L = "dgCMatrix",
                             Pinv = "integer"))

setClass("css_QR", slots = c(Pinv = "integer", Q = "integer",
                             parent = "integer", cp = "integer",
                             nz = "integer"))

setClass("css_LU", slots = c(Q = "integer", nz = "integer"))
}

##-- Schur ---

## non-"Matrix" Class 1  --- For Eigen values:
setClassUnion("number", members = c("numeric", "complex"))

setClass("Schur", contains = "MatrixFactorization",
	 slots = c(T = "Matrix", # <- "block-triangular"; maybe triangular
                   Q = "Matrix", EValues = "number"),
	 validity = function(object) {
	     dim <- object@Dim
	     if((n <- dim[1]) != dim[2])
		 "'Dim' slot is not (n,n)"
	     else if(any(dim(object@T) != n))
		 "'dim(T)' is incorrect"
	     else if(any(dim(object@Q) != n))
		 "'dim(Q)' is incorrect"
	     else if(length(object@EValues) != n)
		 "'EValues' is not of correct length"
	     else TRUE
	 })


### Class Union :  no inheritance, but is(*, <class>) :

setClassUnion("mMatrix", members = c("matrix", "Matrix"))
if(FALSE) ## to be used in setMethod("c", "numM...") -- once that works
setClassUnion("numMatrixLike", members = c("logical", "integer", "numeric", "mMatrix"))

## CARE: Sometimes we'd want all those for which 'x' contains all the data.
##       e.g. Diagonal() is "ddiMatrix" with 'x' slot of length 0, not containing 1
##       same for other  diag="U" (e.g. tridiag)
setClassUnion("xMatrix", ## those Matrix classes with an 'x' slot
              c("dMatrix",
                "iMatrix",
                "lMatrix",
                "ndenseMatrix",
                "zMatrix"))

if(TRUE) { ##--- variant of setClass("dCsparse..." ..) etc working better for other pkgs -----
    ## currently *not* (explicitly) exported

## "classical" Cholmod-like sparseMatrix  (not "indMatrix" or "diagonalMatrix"):
if(FALSE)
  setClassUnion("CRTsparseMatrix", members = c("CsparseMatrix", "RsparseMatrix", "TsparseMatrix"))
  ## would be useful e.g. in ./products.R for "%&%" --- but it changes the method ordering
  ## changing too much (for now)

## These should be *closer* to their members than both {dln]sparse* and Csparse* -- but they are *NOT*
## Could "fix" this be adding these as virtual classes and have the dgC* etc contain *these*
setClassUnion("dCsparseMatrix", members = c("dgCMatrix", "dtCMatrix", "dsCMatrix"))
setClassUnion("lCsparseMatrix", members = c("lgCMatrix", "ltCMatrix", "lsCMatrix"))
setClassUnion("nCsparseMatrix", members = c("ngCMatrix", "ntCMatrix", "nsCMatrix"))

## dense general
setClassUnion("geMatrix", members = c("dgeMatrix", "lgeMatrix", "ngeMatrix"))

## dput(intersect(names(getClass("nsparseMatrix")@subclasses),
##                names(getClass("TsparseMatrix")@subclasses)))
if(FALSE)  ##-- not yet ---------
setClassUnion("nTsparseMatrix", members = c("ngTMatrix", "ntTMatrix", "nsTMatrix"))

}


## --------------------- non-"Matrix" Classes --------------------------------

## --- "General" (not Matrix at all) ----

## e.g. for "Arith" methods, NB: --> see "numericVector" below (incl "integer")
setClassUnion("numLike", members = c("numeric", "logical"))

##setClassUnion("numIndex", members = "numeric")

## Note "rle" is a sealed oldClass (and "virtual" as w/o prototype)
setClass("rleDiff", slots = c(first = "numLike", rle = "rle"),
	 prototype = prototype(first = integer(),
			       rle = rle(integer())),
	 validity = function(object) {
	     if(length(object@first) != 1)
		 return("'first' must be of length one")
	     rl <- object@rle
	     if(!is.list(rl) || length(rl) != 2 ||
		!identical(sort(names(rl)), c("lengths", "values")))
		 return("'rle' must be a list (lengths = *, values = *)")
	     if(length(lens <- rl$lengths) != length(vals <- rl$values))
		 return("'lengths' and 'values' differ in length")
	     if(any(lens <= 0))
		 return("'lengths' must be positive")
	     TRUE
	 })

### 2010-03-04 -- thinking about *implementing* some 'abIndex' methodology,
### I conclude that the following structure would probably be even more
### efficient than the "rleDiff" one :
### IDEA: Store subsequences in a numeric matrix of three rows, where
### ----- one column = [from, to, by]  defining a sub seq()ence

## for now, at least use it, and [TODO!] define  "seqMat" <--> "abIndex" coercions:
setClass("seqMat", contains = "matrix",
	 prototype = prototype(matrix(0, nrow = 3, ncol=0)),
	 validity = function(object) {
	     if(!is.numeric(object)) return("is not numeric")
	     d <- dim(object)
	     if(length(d) != 3 || d[1] != 3)
		 return("not a	 3 x n	matrix")
	     if(any(object != floor(object)))
		 return("some entries are not integer valued")
	     TRUE
	 })

setClass("abIndex", # 'ABSTRact Index'
         slots = c(kind = "character",
                   ## one of ("int32", "double", "rleDiff")
                                        # i.e., numeric or "rleDiff"
                   x = "numLike", # for  numeric [length 0 otherwise]
                   rleD = "rleDiff"),  # "rleDiff" result
         prototype = prototype(kind = "int32", x = integer(0)),# rleD = ... etc
         validity = function(object) {
            switch(object@kind,
                   "int32" = if(!is.integer(object@x))
                   return("'x' slot must be integer when kind is 'int32'")
                   ,
                   "double" = if(!is.double(object@x))
                   return("'x' slot must be double when kind is 'double'")
                   ,
                   "rleDiff" = {
                       if(length(object@x))
                   return("'x' slot must be empty when kind is 'rleDiff'")
                   },
                   ## otherwise
                   return("'kind' must be one of (\"int32\", \"double\", \"rleDiff\")")
                   )
            TRUE
         })

## for 'i' in x[i] or A[i,] : (numeric = {double, integer})
## TODO: allow "abIndex" as well !
setClassUnion("index", members =  c("numeric", "logical", "character"))

## "atomic vectors" (-> ?is.atomic ) -- but note that is.atomic(<matrix>) !
## ---------------  those that we want to convert from old-style "matrix"
setClassUnion("atomicVector", ## "double" is not needed, and not liked by some
	      members = c("logical", "integer", "numeric",
			  "complex", "raw", "character"))

## NB: --> see "numLike" above
if(FALSE) # not used anywhere
setClassUnion("numericVector", members = c("logical", "integer", "numeric"))




## --- Matrix - related (but not "Matrix" nor "Decomposition/Factorization):

### Sparse Vectors ---- here use 1-based indexing ! -----------

## 'longindex' should allow sparseVectors of "length" > 2^32,
## which is necessary e.g. when converted from large sparse matrices
## setClass("longindex", contains = "numeric")
## but we use "numeric" instead, for simplicity (efficiency?)
setClass("sparseVector",
         slots = c(length = "numeric", i = "numeric"), contains = "VIRTUAL",
         ##                     "longindex"    "longindex"
         ## note that "numeric" contains "integer" (if I like it or not..)
	 prototype = prototype(length = 0),
         validity = function(object) {
	     n <- object@length
	     if(anyNA(i <- object@i))	 "'i' slot has NAs"
	     else if(any(!is.finite(i))) "'i' slot is not all finite"
	     else if(any(i < 1))	 "'i' must be >= 1"
	     else if(n == 0 && length(i))"'i' must be empty when the object length is zero"
	     else if(any(i > n)) sprintf("'i' must be in 1:%d", n)
	     else if(is.unsorted(i, strictly=TRUE))
		 "'i' must be sorted strictly increasingly"
             else TRUE
         })

##' initialization -- ensuring that  'i' is sorted (and 'x' alongside)
setMethod("initialize", "sparseVector",
          function(.Object, i, x, ...) {
              has.x <- !missing(x)
              if(!missing(i)) {
                  i <- ## (be careful to assign in all cases)
                      if(is.unsorted(i, strictly=TRUE)) {
                          if(is(.Object, "xsparseVector") && has.x) {
                              si <- sort.int(i, index.return=TRUE)
                              x <- x[si$ix]
                              si$x
                          } else sort.int(i, method = "quick")
                      } else i
              }
              if(has.x)
                  x <- x
              callNextMethod()
          })

.validXspVec <- function(object) {
    ## n <- object@length
    if(length(object@i) != length(object@x))
        "'i' and 'x' differ in length"
    else TRUE
}
setClass("dsparseVector",
	 slots = c(x = "numeric"), contains = "sparseVector",
	 validity = .validXspVec)
setClass("isparseVector",
	 slots = c(x = "integer"), contains = "sparseVector",
	 validity = .validXspVec)
setClass("lsparseVector",
	 slots = c(x = "logical"), contains = "sparseVector",
	 validity = .validXspVec)
setClass("zsparseVector",
	 slots = c(x = "complex"), contains = "sparseVector",
	 validity = .validXspVec)
## nsparse has no new slot: 'i' just contains the locations!
setClass("nsparseVector", contains = "sparseVector")

setClassUnion("xsparseVector", ## those sparseVector's with an 'x' slot
              c("dsparseVector",
                "isparseVector",
                "lsparseVector",
                "zsparseVector"))

## for 'value' in  x[..] <- value hence for all "contents" of our Matrices:
setClassUnion("replValue",   members = c("numeric", "logical", "complex", "raw"))
setClassUnion("replValueSp", members = c("replValue", "sparseVector", "Matrix"))


setClass("determinant",
	 slots = c(modulus = "numeric",
		   logarithm = "logical",
		   sign = "integer",
		   call = "call"))

## --- New "logic" class -- currently using "raw" instead of "logical"
## LOGIC setClass("logic", contains = "raw")
