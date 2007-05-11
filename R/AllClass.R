.onLoad <- function(lib, pkg) {
    if(is.null(getOption("max.print")))
	options(max.print = 10000)#-> show() of large matrices
}

## --- New "logic" class -- currently using "raw" instead of "logical"
## LOGIC setClass("logic", contains = "raw")

## ------------- Virtual Classes ----------------------------------------

## Mother class of all Matrix objects
setClass("Matrix",
	 representation(Dim = "integer", Dimnames = "list", "VIRTUAL"),
	 prototype = prototype(Dim = integer(2), Dimnames = list(NULL,NULL)),
	 validity = function(object) {
	     Dim <- object@Dim
	     if (length(Dim) != 2)
		 return("Dim slot must be of length 2")
	     if (any(Dim < 0))
		 return("Dim slot must contain non-negative values")
	     Dn <- object@Dimnames
	     if (!is.list(Dn) || length(Dn) != 2)
		 return("'Dimnames' slot must be list of length 2")
	     lDn <- sapply(Dn, length)
	     if (lDn[1] > 0 && lDn[1] != Dim[1])
		 return("length(Dimnames[[1]])' must match Dim[1]")
	     if (lDn[2] > 0 && lDn[2] != Dim[2])
		 return("length(Dimnames[[2]])' must match Dim[2]")
	     ## 'else'	ok :
	     TRUE
	 })

## The class of composite matrices - i.e. those for which it makes sense to
## create a factorization
setClass("compMatrix",	representation(factors = "list", "VIRTUAL"),
	 contains = "Matrix")

## Virtual classes of Matrices determined by above/below diagonal relationships

setClass("generalMatrix", representation = "VIRTUAL", contains = "compMatrix")

setClass("symmetricMatrix",
	 representation(uplo = "character", "VIRTUAL"),
	 contains = "compMatrix",
	 prototype = prototype(uplo = "U"),
	 validity = function(object) .Call(symmetricMatrix_validate, object))

setClass("triangularMatrix",
	 representation(uplo = "character", diag = "character", "VIRTUAL"),
	 contains = "Matrix",
	 prototype = prototype(uplo = "U", diag = "N"),
	 validity = function(object) .Call(triangularMatrix_validate, object))


## Virtual class of numeric matrices
setClass("dMatrix",
	 representation(x = "numeric", "VIRTUAL"), contains = "Matrix",
	 validity = function(object)
	 .Call(dMatrix_validate, object))

## Virtual class of integer matrices
setClass("iMatrix",
	 representation(x = "integer", "VIRTUAL"), contains = "Matrix")

## Virtual class of logical matrices
setClass("lMatrix",
## LOGIC representation(x = "logic", "VIRTUAL"), contains = "Matrix")
         representation(x = "logical", "VIRTUAL"), contains = "Matrix")

## Virtual class of nonzero pattern matrices
setClass("nMatrix", representation("VIRTUAL"), contains = "Matrix")
## aka 'pattern' matrices -- have no x slot

## Virtual class of complex matrices
setClass("zMatrix", # letter 'z' is as in the names of Lapack subroutines
	 representation(x = "complex", "VIRTUAL"), contains = "Matrix")

## Virtual class of dense matrices (including "packed")
setClass("denseMatrix", representation("VIRTUAL"),
	 contains = "Matrix")

## Virtual class of dense, numeric matrices
setClass("ddenseMatrix", representation("VIRTUAL"),
	 contains = c("dMatrix", "denseMatrix"))

## Virtual class of dense, logical matrices
setClass("ldenseMatrix", representation("VIRTUAL"),
	 contains = c("lMatrix", "denseMatrix"))

## Virtual class of dense, nonzero pattern matrices - rarely used, for completeness
setClass("ndenseMatrix", representation(x = "logical", "VIRTUAL"),
	 contains = c("nMatrix", "denseMatrix"))

## diagonal: has 'diag' slot;  diag = "U"  <--> have identity matrix
setClass("diagonalMatrix", representation(diag = "character", "VIRTUAL"),
	 contains = "denseMatrix",
	 validity = function(object) {
	     d <- object@Dim
	     if(d[1] != (n <- d[2])) return("matrix is not square")
	     lx <- length(object@x)
	     if(object@diag == "U" && lx != 0)
		 return("diag = \"U\" (identity matrix) requires empty 'x' slot")
	     if(object@diag == "N" && lx != n)
		 return("diagonal matrix has 'x' slot of length != 'n'")
	     TRUE
	 },
	 prototype = prototype(diag = "N")
	 )

## virtual SPARSE ------------

setClass("sparseMatrix", representation("VIRTUAL"), contains = "Matrix")

## sparse matrices in Triplet representation (dgT, lgT, ..):
setClass("TsparseMatrix", representation(i = "integer", j = "integer", "VIRTUAL"),
	 contains = "sparseMatrix",
	 validity = function(object) .Call(Tsparse_validate, object)
         )

setClass("CsparseMatrix", representation(i = "integer", p = "integer", "VIRTUAL"),
	 contains = "sparseMatrix",
	 prototype = prototype(p = 0:0),# to be valid
         validity = function(object) .Call(Csparse_validate, object)
         )

setClass("RsparseMatrix", representation(p = "integer", j = "integer", "VIRTUAL"),
	 contains = "sparseMatrix",
	 prototype = prototype(p = 0:0),# to be valid
	 ## TODO:
	 ## , validity = function(object) .Call(Rsparse_validate, object)
         )

setClass("dsparseMatrix", representation("VIRTUAL"),
	 contains = c("dMatrix", "sparseMatrix"))

setClass("lsparseMatrix", representation("VIRTUAL"),
	 contains = c("lMatrix", "sparseMatrix"))

## these are the "pattern" matrices for "symbolic analysis" of sparse OPs:
setClass("nsparseMatrix", representation("VIRTUAL"),
	 contains = c("nMatrix", "sparseMatrix"))

## ------------------ Proper (non-virtual) Classes ----------------------------

##----------------------  DENSE	 -----------------------------------------

## numeric, dense, general matrices
setClass("dgeMatrix", contains = c("ddenseMatrix", "generalMatrix"),
	 ## checks that length( @ x) == prod( @ Dim):
	 validity =
	 function(object) .Call(dgeMatrix_validate, object)
	 )
## i.e. "dgeMatrix" cannot be packed, but "ddenseMatrix" can ..

## numeric, dense, non-packed, triangular matrices
setClass("dtrMatrix",
	 contains = c("ddenseMatrix", "triangularMatrix"),
	 validity =
	 function(object) .Call(dtrMatrix_validate, object)
	 )

## numeric, dense, packed, triangular matrices
setClass("dtpMatrix",
	 contains = c("ddenseMatrix", "triangularMatrix"),
	 validity =
	 function(object) .Call(dtpMatrix_validate, object)
	 )

## numeric, dense, non-packed symmetric matrices
setClass("dsyMatrix",
         contains = c("ddenseMatrix", "symmetricMatrix"),
	 validity =
	 function(object) .Call(dsyMatrix_validate, object)
	 )

## numeric, dense, packed symmetric matrices
setClass("dspMatrix",
	 contains = c("ddenseMatrix", "symmetricMatrix"),
	 validity =
	 function(object) .Call(dspMatrix_validate, object)
	 )

## numeric, dense, non-packed, positive-definite, symmetric matrices
setClass("dpoMatrix", contains = "dsyMatrix",
	 validity = function(object) .Call(dpoMatrix_validate, object)
	 )

## numeric, dense, packed, positive-definite, symmetric matrices
setClass("dppMatrix", contains = "dspMatrix",
	 validity = function(object) .Call(dppMatrix_validate, object)

)
##----- logical dense Matrices -- e.g. as result of <ddenseMatrix>  COMPARISON

## logical, dense, general matrices
setClass("lgeMatrix", contains = c("ldenseMatrix", "generalMatrix"),
	 ## checks that length( @ x) == prod( @ Dim):
	 validity = function(object) stopifnot(length(object@x) == prod(object@Dim))
	 )
## i.e. "lgeMatrix" cannot be packed, but "ldenseMatrix" can ..

## logical, dense, non-packed, triangular matrices
setClass("ltrMatrix",
	 contains = c("ldenseMatrix", "triangularMatrix"))

## logical, dense, packed, triangular matrices
setClass("ltpMatrix",
	 contains = c("ldenseMatrix", "triangularMatrix"))

## logical, dense, non-packed symmetric matrices
setClass("lsyMatrix",
	 contains = c("ldenseMatrix", "symmetricMatrix"))

## logical, dense, packed symmetric matrices
setClass("lspMatrix",
	 contains = c("ldenseMatrix", "symmetricMatrix"),
	 validity = function(object) .Call(dspMatrix_validate, object)
	 ## "dsp" and "lsp" have the same validate
	 )

##----- nonzero pattern dense Matrices -- "for completeness"

## logical, dense, general matrices
setClass("ngeMatrix", contains = c("ndenseMatrix", "generalMatrix"),
	 ## checks that length( @ x) == prod( @ Dim):
	 validity = function(object)
         stopifnot(length(object@x) == prod(object@Dim)))
## i.e. "ngeMatrix" cannot be packed, but "ndenseMatrix" can ..

## logical, dense, non-packed, triangular matrices
setClass("ntrMatrix",
	 contains = c("ndenseMatrix", "triangularMatrix"))

## logical, dense, packed, triangular matrices
setClass("ntpMatrix",
	 contains = c("ndenseMatrix", "triangularMatrix"))

## logical, dense, non-packed symmetric matrices
setClass("nsyMatrix",
	 contains = c("ndenseMatrix", "symmetricMatrix"))

## logical, dense, packed symmetric matrices
setClass("nspMatrix",
	 contains = c("ndenseMatrix", "symmetricMatrix"),
	 validity = function(object) .Call(dspMatrix_validate, object)
	 ## "dsp" and "nsp" have the same validate
	 )


## 'diagonalMatrix' already has validity checking
## diagonal, numeric matrices;	    "d*" has 'x' slot :
setClass("ddiMatrix", contains = c("diagonalMatrix", "ddenseMatrix"))# or "dMatrix"
## diagonal, logical matrices; "ldense*" has 'x' slot :
setClass("ldiMatrix", contains = c("diagonalMatrix", "ldenseMatrix"))

setClass("corMatrix", representation(sd = "numeric"), contains = "dpoMatrix",
	 validity = function(object) {
	     ## assuming that 'dpoMatrix' validity check has already happened:
	     n <- object@Dim[2]
	     if(length(sd <- object@sd) != n)
		 return("'sd' slot must be of length 'dim(.)[1]'")
	     if(any(sd < 0))
		 return("'sd' slot has negative entries")
	     TRUE
	 })


##-------------------- S P A R S E (non-virtual) --------------------------

##---------- numeric sparse matrix classes --------------------------------

## numeric, sparse, triplet general matrices
setClass("dgTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xTMatrix_validate, object)
	 )

## Should not have dtTMatrix inherit from dgTMatrix because a dtTMatrix could
## be less than fully stored if diag = "U".  Methods for the dgTMatrix
## class would not produce correct results even though all the slots
## are present.

## numeric, sparse, triplet triangular matrices
setClass("dtTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(tTMatrix_validate, object)
	 )

## numeric, sparse, triplet symmetric matrices(also only store one triangle)
setClass("dsTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
	 validity = function(object) .Call(tTMatrix_validate, object)
	 )

## numeric, sparse, sorted compressed sparse column-oriented general matrices
setClass("dgCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xCMatrix_validate, object)
	 )

## see comments for dtTMatrix above
## numeric, sparse, sorted compressed sparse column-oriented triangular matrices
setClass("dtCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(tCMatrix_validate, object)
	 )

## see comments for dsTMatrix above
## numeric, sparse, sorted compressed sparse column-oriented symmetric matrices
setClass("dsCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
	 validity = function(object) .Call(tCMatrix_validate, object)
	 )

## numeric, sparse, sorted compressed sparse row-oriented general matrices
setClass("dgRMatrix",
	 contains = c("RsparseMatrix", "dsparseMatrix", "generalMatrix"),
	 ##TODO: validity = function(object) .Call(dgRMatrix_validate, object)
	 )

## numeric, sparse, sorted compressed sparse row-oriented triangular matrices
setClass("dtRMatrix",
	 contains = c("RsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 ##TODO: validity = function(object) .Call(dtRMatrix_validate, object)

	 )

## numeric, sparse, sorted compressed sparse row-oriented symmetric matrices
setClass("dsRMatrix",
	 contains = c("RsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
	 ##TODO: validity = function(object) .Call(dsRMatrix_validate, object)
	 )

##---------- logical sparse matrix classes --------------------------------

## these classes are used in symbolic analysis to determine the
## locations of non-zero entries

## logical, sparse, triplet general matrices
setClass("lgTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xTMatrix_validate, object)
	 )

## logical, sparse, triplet triangular matrices
setClass("ltTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(xTMatrix_validate, object)
	 )

## logical, sparse, triplet symmetric matrices
setClass("lsTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
	 validity = function(object) .Call(xTMatrix_validate, object)
	 )

## logical, sparse, sorted compressed sparse column-oriented general matrices
setClass("lgCMatrix",
	 contains = c("CsparseMatrix", "lsparseMatrix", "generalMatrix"),
	 validity = function(object) .Call(xCMatrix_validate, object)
	 )

## logical, sparse, sorted compressed sparse column-oriented triangular matrices
setClass("ltCMatrix",
	 contains = c("CsparseMatrix", "lsparseMatrix", "triangularMatrix"),
	 validity = function(object) .Call(xCMatrix_validate, object)
	 )

## logical, sparse, sorted compressed sparse column-oriented symmetric matrices
setClass("lsCMatrix",
	 contains = c("CsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
	 validity = function(object) .Call(xCMatrix_validate, object)
	 )

## logical, sparse, sorted compressed sparse row-oriented general matrices
setClass("lgRMatrix",
	 representation(j = "integer", p = "integer"),
	 contains = c("RsparseMatrix", "lsparseMatrix", "generalMatrix"),
	 ##TODO: validity = function(object) .Call(lgRMatrix_validate, object)
	 )

## logical, sparse, sorted compressed sparse row-oriented triangular matrices
setClass("ltRMatrix",
	 contains = c("RsparseMatrix", "lsparseMatrix", "triangularMatrix"),
	 ##TODO: validity = function(object) .Call(ltRMatrix_validate, object)
	 )

## logical, sparse, sorted compressed sparse row-oriented symmetric matrices
setClass("lsRMatrix",
	 contains = c("RsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
	 ##TODO: validity = function(object) .Call(lsRMatrix_validate, object)
	 )

##---------- nonzero pattern sparse matrix classes ---------------------------

## these classes are used in symbolic analysis to determine the
## locations of non-zero entries

## nonzero pattern, sparse, triplet general matrices
setClass("ngTMatrix",
	 contains = c("TsparseMatrix", "nsparseMatrix", "generalMatrix")
         ## validity: Tsparse_validate should be enough
	 )

## nonzero pattern, sparse, triplet triangular matrices
setClass("ntTMatrix",
	 contains = c("TsparseMatrix", "nsparseMatrix", "triangularMatrix"),
         ## validity: Tsparse_ and triangular*_validate should be enough
	 )

## nonzero pattern, sparse, triplet symmetric matrices
setClass("nsTMatrix",
	 contains = c("TsparseMatrix", "nsparseMatrix", "symmetricMatrix"),
         ## validity: Tsparse_ and symmetric*_validate should be enough
	 )

## nonzero pattern, sparse, sorted compressed column-oriented general matrices
setClass("ngCMatrix",
	 contains = c("CsparseMatrix", "nsparseMatrix", "generalMatrix"),
         ## validity: Csparse_validate should be enough
	 )

## nonzero pattern, sparse, sorted compressed column-oriented triangular matrices
setClass("ntCMatrix",
	 contains = c("CsparseMatrix", "nsparseMatrix", "triangularMatrix"),
         ## validity: Csparse_ and triangular*_validate should be enough
	 )

## nonzero pattern, sparse, sorted compressed column-oriented symmetric matrices
setClass("nsCMatrix",
	 contains = c("CsparseMatrix", "nsparseMatrix", "symmetricMatrix"),
         ## validity: Csparse_ and symmetric*_validate should be enough
	 )

## nonzero pattern, sparse, sorted compressed row-oriented general matrices
setClass("ngRMatrix",
	 contains = c("RsparseMatrix", "nsparseMatrix", "generalMatrix"),
	 )

## nonzero pattern, sparse, sorted compressed row-oriented triangular matrices
setClass("ntRMatrix",
	 contains = c("RsparseMatrix", "nsparseMatrix", "triangularMatrix"),
	 )

## nonzero pattern, sparse, sorted compressed row-oriented symmetric matrices
setClass("nsRMatrix",
	 contains = c("RsparseMatrix", "nsparseMatrix", "symmetricMatrix"),
	 )

##-------------------- permutation ----------------------------------------

setClass("pMatrix", representation(perm = "integer"),
	 contains = c("sparseMatrix", "generalMatrix"),
	 validity = function(object) {
	     d <- object@Dim
	     if (d[2] != (n <- d[1])) return("pMatrix must be square")
	     perm <- object@perm
	     if (length(perm) != n)
		 return(paste("length of 'perm' slot must be", n))
	     if(n > 0 &&
		!(all(range(perm) == c(1, n)) && length(unique(perm)) == n))
		 return("'perm' slot is not a valid permutation")
	     TRUE
	 })


### Factorization classes ---------------------------------------------

## Mother class:
setClass("MatrixFactorization", representation(Dim = "integer", "VIRTUAL"))

## -- Those (exceptions) inheriting from "Matrix" : ---

setClass("Cholesky",  contains = c("dtrMatrix", "MatrixFactorization"))

#unUsed: setClass("LDL", contains = c("dtrMatrix", "MatrixFactorization"))

setClass("pCholesky", contains = c("dtpMatrix", "MatrixFactorization"))

## These are currently only produced implicitly from *solve()
setClass("BunchKaufman",
	 contains = c("dtrMatrix", "MatrixFactorization"),
	 representation(perm = "integer"),
	 validity =
	 function(object) .Call(BunchKaufman_validate, object)
	 )

setClass("pBunchKaufman",
	 contains = c("dtpMatrix", "MatrixFactorization"),
	 representation(perm = "integer"),
	 validity =
	 function(object) .Call(pBunchKaufman_validate, object))

## -- the usual ``non-Matrix'' factorizations : ---------

setClass("CHMfactor",		 # cholmod_factor struct as S4 object
	 contains = "MatrixFactorization",
	 representation(colcount = "integer", perm = "integer",
			type = "integer", "VIRTUAL"),
	 validity = function(object) .Call(CHMfactor_validate, object))

setClass("CHMsuper",		       # supernodal cholmod_factor
	 contains = "CHMfactor",
	 representation(super = "integer", pi = "integer", px = "integer",
			s = "integer", "VIRTUAL"),
	 validity = function(object) .Call(CHMsuper_validate, object))

setClass("CHMsimpl",		       # simplicial cholmod_factor
	 contains = "CHMfactor",
	 representation(p = "integer", i = "integer", nz = "integer",
			nxt = "integer", prv = "integer", "VIRTUAL"),
	 validity = function(object) .Call(CHMsimpl_validate, object))

setClass("dCHMsuper", contains = "CHMsuper", representation(x = "numeric"))

setClass("nCHMsuper", contains = "CHMsuper")

setClass("dCHMsimpl", contains = "CHMsimpl", representation(x = "numeric"))

setClass("nCHMsimpl", contains = "CHMsimpl")

##--- LU ---

setClass("LU", contains = "MatrixFactorization", representation("VIRTUAL"))

setClass("denseLU", contains = "LU",
	 representation(x = "numeric", perm = "integer"),
	 validity = function(object) .Call(LU_validate, object))

setClass("sparseLU", contains = "LU",
	 representation(L = "dgCMatrix", U = "dgCMatrix",
			p = "integer", q = "integer"))

##--- QR ---

setClass("sparseQR", contains = "MatrixFactorization",
	 representation(V = "dgCMatrix", beta = "numeric",
			p = "integer", R = "dgCMatrix", q = "integer"))

## "denseQR" -- ?  (``a version of''  S3 class "qr")

if (FALSE) { ## unused classes
setClass("csn_QR", representation(U = "dgCMatrix", L = "dgCMatrix",
                                  beta = "numeric"))

setClass("csn_LU", representation(U = "dgCMatrix", L = "dgCMatrix",
                                  Pinv = "integer"))

setClass("css_QR", representation(Pinv = "integer", Q = "integer",
                                  parent = "integer", cp = "integer",
                                  nz = "integer"))

setClass("css_LU", representation(Q = "integer", nz = "integer"))
}


### Class Union :  no inheritance, but is(*, <class>) :

## Definition  Packed := dense with length( . @x) < prod( . @Dim)
##	       ~~~~~~
## REPLACED the following with	isPacked() in ./Auxiliaries.R :
## setClassUnion("packedMatrix",
##		 members = c("dspMatrix", "dppMatrix", "dtpMatrix",
##		  "lspMatrix", "ltpMatrix", "diagonalMatrix"))


## --------------------- non-"Matrix" Classes --------------------------------

## --- "General" (not Matrix at all) ----

## for 'i' in x[i] or A[i,] : (numeric = {double, integer})
setClassUnion("index", members =  c("numeric", "logical", "character"))

## "atomic vectors" (-> ?is.atomic ) --
## ---------------  those that we want to convert from old-style "matrix"
setClassUnion("atomicVector", ## numeric = {integer, double} but all 3 should *directly* be atomic
	      members = c("logical", "integer", "double", "numeric",
			  "complex", "raw", "character"))

## --- Matrix - related (but not "Matrix" nor "Decomposition/Factorization):

### for 'value' in  x[..] <- value hence for all "contents" of our Matrices:
setClassUnion("replValue", members =  c("numeric", "logical", "complex", "raw"))


setClass("determinant",
	 representation(modulus = "numeric",
			logarithm = "logical",
			sign = "integer",
			call = "call"))

