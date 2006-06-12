.onLoad <- function(lib, pkg) {
    if(is.null(getOption("max.print")))
	options(max.print = 10000)#-> show() of large matrices
}

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
	 contains = "compMatrix")

setClass("triangularMatrix",
	 representation(uplo = "character", diag = "character", "VIRTUAL"),
	 contains = "Matrix",
	 validity = function(object)
	 .Call(triangularMatrix_validate, object)
	 )


## Virtual class of numeric matrices
setClass("dMatrix",
	 representation(x = "numeric", "VIRTUAL"), contains = "Matrix",
	 validity = function(object)
	 .Call(dMatrix_validate, object))

## Virtual class of integer matrices
setClass("iMatrix",
	 representation(x = "integer", "VIRTUAL"), contains = "Matrix")

## Virtual class of logical matrices
setClass("lMatrix", representation("VIRTUAL"), contains = "Matrix")
## Note that logical sparse matrices do not need an x slot so the x
## slot is part of the ldenseMatrix class

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
setClass("ldenseMatrix", representation(x = "logical", "VIRTUAL"),
	 contains = c("lMatrix", "denseMatrix"))

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
	 contains = "sparseMatrix")

setClass("CsparseMatrix", representation(i = "integer", p = "integer", "VIRTUAL"),
	 contains = "sparseMatrix")

setClass("RsparseMatrix", representation(p = "integer", j = "integer", "VIRTUAL"),
	 contains = "sparseMatrix")

setClass("dsparseMatrix", representation("VIRTUAL"),
	 contains = c("dMatrix", "sparseMatrix"))

setClass("lsparseMatrix", representation("VIRTUAL"),
	 contains = c("lMatrix", "sparseMatrix"))

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
	 prototype = prototype(uplo = "U", diag = "N"),
	 validity =
	 function(object) .Call(dtrMatrix_validate, object)
	 )

## numeric, dense, packed, triangular matrices
setClass("dtpMatrix",
	 contains = c("ddenseMatrix", "triangularMatrix"),
	 prototype = prototype(uplo = "U", diag = "N"),
	 validity =
	 function(object) .Call(dtpMatrix_validate, object)
	 )

## numeric, dense, non-packed symmetric matrices
setClass("dsyMatrix",
	 ## FIXME?
	 ##> 'ddense*' before 'dge*' so it can use d* or ddense* methods
	 ##> WITHOUT a coerce to dge* (losing symmetry)
	 ##> gives error in crossprod() dispatch
	 ##> contains = c("ddenseMatrix", "dgeMatrix", "symmetricMatrix"),
	 contains = c("ddenseMatrix", "symmetricMatrix"),
	 prototype = prototype(uplo = "U"),
	 validity =
	 function(object) .Call(dsyMatrix_validate, object)
	 )

## numeric, dense, packed symmetric matrices
setClass("dspMatrix",
	 prototype = prototype(uplo = "U"),
	 contains = c("ddenseMatrix", "symmetricMatrix"),
	 validity =
	 function(object) .Call(dspMatrix_validate, object)
	 )

## numeric, dense, non-packed, positive-definite, symmetric matrices
setClass("dpoMatrix", contains = "dsyMatrix",
	 validity =
	 function(object) .Call(dpoMatrix_validate, object)
	 )

## numeric, dense, packed, positive-definite, symmetric matrices
setClass("dppMatrix", contains = "dspMatrix",
	 validity =
	 function(object) .Call(dppMatrix_validate, object)
	 )

##----- logical dense Matrices -- e.g. as result of <ddenseMatrix>  COMPARISON

## numeric, dense, general matrices
setClass("lgeMatrix", contains = c("ldenseMatrix", "generalMatrix"),
	 ## checks that length( @ x) == prod( @ Dim):
	 validity = function(object) stopifnot(length(object@x) == prod(object@Dim))
	 )
## i.e. "lgeMatrix" cannot be packed, but "ldenseMatrix" can ..

## numeric, dense, non-packed, triangular matrices
setClass("ltrMatrix",
	 contains = c("ldenseMatrix", "triangularMatrix"),
	 prototype = prototype(uplo = "U", diag = "N"))

## numeric, dense, packed, triangular matrices
setClass("ltpMatrix",
	 contains = c("ldenseMatrix", "triangularMatrix"),
	 prototype = prototype(uplo = "U", diag = "N")
	 ## validity: ldense*, triangular*  should suffice
	 )

## numeric, dense, non-packed symmetric matrices
setClass("lsyMatrix",
	 contains = c("ldenseMatrix", "symmetricMatrix"),
	 prototype = prototype(uplo = "U")
	 ##, validity = function(object) .Call(lsyMatrix_validate, object)
	 )

## numeric, dense, packed symmetric matrices
setClass("lspMatrix",
	 contains = c("ldenseMatrix", "symmetricMatrix"),
	 prototype = prototype(uplo = "U"),
	 validity = function(object)
	 .Call(dspMatrix_validate, object)
	 ## "dsp" and "lsp" have the same validate
	 )

## 'diagonalMatrix' already has validity checking
## diagonal, numeric matrices;	    "d*" has 'x' slot :
setClass("ddiMatrix", contains = c("diagonalMatrix", "ddenseMatrix"))# or "dMatrix"
## diagonal, logical matrices; "ldense*" has 'x' slot :
setClass("ldiMatrix", contains = c("diagonalMatrix", "ldenseMatrix"))


##-------------------- S P A R S E (non-virtual) --------------------------

##---------- numeric sparse matrix classes --------------------------------

## numeric, sparse, triplet general matrices
setClass("dgTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "generalMatrix"),
	 validity =
	 function(object) .Call(dgTMatrix_validate, object)
	 )

## Should not have dtTMatrix inherit from dgTMatrix because a dtTMatrix could
## be less than fully stored if diag = "U".  Methods for the dgTMatrix
## class would not produce correct results even though all the slots
## are present.

## numeric, sparse, triplet triangular matrices
setClass("dtTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 prototype = prototype(uplo = "U", diag = "N"),
	 validity =
	 function(object) .Call(dtTMatrix_validate, object)
	 )

## Should not have dsTMatrix inherit from dgTMatrix because a dsTMatrix
## is not fully stored.	 Methods for the dgTMatrix class would not
## produce correct results even though all the slots are present.

## numeric, sparse, triplet symmetric matrices
setClass("dsTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
	 prototype = prototype(uplo = "U"),
	 validity =
	 function(object) .Call(dsTMatrix_validate, object)
	 )

## numeric, sparse, sorted compressed sparse column-oriented general matrices
setClass("dgCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "generalMatrix"),
	 prototype = prototype(p = 0:0),# to be valid
	 validity =
	 function(object) .Call(dgCMatrix_validate, object)
	 )

## see comments for dtTMatrix above
## numeric, sparse, sorted compressed sparse column-oriented triangular matrices
setClass("dtCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 prototype = prototype(p = 0:0, uplo = "U", diag = "N"),# to be valid
	 validity =
	 function(object) .Call(tsc_validate, object)
	 )

## see comments for dsTMatrix above
## numeric, sparse, sorted compressed sparse column-oriented symmetric matrices
setClass("dsCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
	 prototype = prototype(p = 0:0, uplo = "U"),# to be valid
	 validity =
	 function(object) .Call(dsCMatrix_validate, object)
	 )

## numeric, sparse, sorted compressed sparse row-oriented general matrices
setClass("dgRMatrix",
	 contains = c("RsparseMatrix", "dsparseMatrix", "generalMatrix"),
	 prototype = prototype(p = 0:0),
	 ##TODO: validity = function(object) .Call(dgRMatrix_validate, object)
	 )

## numeric, sparse, sorted compressed sparse row-oriented triangular matrices
setClass("dtRMatrix",
	 contains = c("RsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 prototype = prototype(p = 0:0, uplo = "U", diag = "N"),# to be valid
	 ##TODO: validity = function(object) .Call(dtRMatrix_validate, object)

	 )

## numeric, sparse, sorted compressed sparse row-oriented symmetric matrices
setClass("dsRMatrix",
	 contains = c("RsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
	 prototype = prototype(p = 0:0, uplo = "U"),# to be valid
	 ##TODO: validity = function(object) .Call(dsRMatrix_validate, object)
	 )

##---------- logical sparse matrix classes --------------------------------

## these classes are used in symbolic analysis to determine the
## locations of non-zero entries

## logical, sparse, triplet general matrices
setClass("lgTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "generalMatrix"),
	 validity =
	 function(object) .Call(lgTMatrix_validate, object)
	 )

## logical, sparse, triplet triangular matrices
setClass("ltTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "triangularMatrix"),
	 prototype = prototype(uplo = "U", diag = "N"),
	 validity =
	 function(object) .Call(ltTMatrix_validate, object)
	 )

## logical, sparse, triplet symmetric matrices
setClass("lsTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
	 validity =
	 function(object) .Call(lsTMatrix_validate, object)
	 )

## logical, sparse, sorted compressed sparse column-oriented general matrices
setClass("lgCMatrix",
	 contains = c("CsparseMatrix", "lsparseMatrix", "generalMatrix"),
	 prototype = prototype(p = 0:0),# to be valid
	 validity =
	 function(object) .Call(lgCMatrix_validate, object)
	 )

## logical, sparse, sorted compressed sparse column-oriented triangular matrices
setClass("ltCMatrix",
	 contains = c("CsparseMatrix", "lsparseMatrix", "triangularMatrix"),
	 prototype = prototype(p = 0:0, uplo = "U", diag = "N"),# to be valid
	 validity =
	 function(object) .Call(ltCMatrix_validate, object)
	 )

## logical, sparse, sorted compressed sparse column-oriented symmetric matrices
setClass("lsCMatrix",
	 contains = c("CsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
	 prototype = prototype(p = 0:0, uplo = "U"),# to be valid
	 validity =
	 function(object) .Call(lsCMatrix_validate, object)
	 )

## logical, sparse, sorted compressed sparse row-oriented general matrices
setClass("lgRMatrix",
	 representation(j = "integer", p = "integer"),
	 contains = c("RsparseMatrix", "lsparseMatrix", "generalMatrix"),
	 validity =
	 function(object) .Call(lgRMatrix_validate, object)
	 )

## logical, sparse, sorted compressed sparse row-oriented triangular matrices
setClass("ltRMatrix",
	 contains = c("RsparseMatrix", "lsparseMatrix", "triangularMatrix"),
	 validity =
	 function(object) .Call(ltRMatrix_validate, object)
	 )

## logical, sparse, sorted compressed sparse row-oriented symmetric matrices
setClass("lsRMatrix",
	 contains = c("RsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
	 validity =
	 function(object) .Call(lsRMatrix_validate, object)
	 )

### Factorization classes ---------------------------------------------

setClass("Cholesky", contains = "dtrMatrix")

setClass("LDL", contains = "dtrMatrix")

setClass("correlation", representation(sd = "numeric"), contains = "dpoMatrix")

setClass("pCholesky", contains = "dtpMatrix")

setClass("BunchKaufman",
	 representation(perm = "integer"),
	 contains = "dtrMatrix",
	 validity =
	 function(object) .Call(BunchKaufman_validate, object)
	 )

setClass("pBunchKaufman",
	 representation(perm = "integer"),
	 contains = "dtpMatrix",
	 validity =
	 function(object) .Call(pBunchKaufman_validate, object)
	 )

setClass("dCholCMatrix",
	 representation(perm = "integer", Parent = "integer", D = "numeric"),
	 contains = "dtCMatrix",
	 validity =
	 function(object) .Call(dCholCMatrix_validate, object)
	 )

setClass("lCholCMatrix",
	 representation(perm = "integer", Parent = "integer"),
	 contains = "ltCMatrix",
	 validity =
	 function(object) .Call(lCholCMatrix_validate, object)
	 )

setClass("CHMfactor",		 # cholmod_factor struct as S4 object
	 representation(colcount = "integer", perm = "integer",
                        type = "integer", "VIRTUAL"),
	 validity =
	 function(object) .Call(CHMfactor_validate, object)
	 )

setClass("CHMsuper",		       # supernodal cholmod_factor
	 representation(super = "integer", pi = "integer", px = "integer",
			s = "integer", "VIRTUAL"),
	 contains = "CHMfactor",
	 validity =
	 function(object) .Call(CHMsuper_validate, object))

setClass("CHMsimpl",		       # simplicial cholmod_factor
	 representation(p = "integer", i = "integer",
			nz = "integer", nxt = "integer", prv = "integer", "VIRTUAL"),
	 contains = "CHMfactor",
	 validity =
	 function(object) .Call(CHMsuper_validate, object))

setClass("dCHMsuper", representation(x = "numeric"), contains = "CHMsuper")

setClass("lCHMsuper", contains = "CHMsuper")

setClass("dCHMsimpl", representation(x = "numeric"), contains = "CHMsimpl")

setClass("lCHMsimpl", contains = "CHMsimpl")


##-------------------- permutation ----------------------------------------

setClass("pMatrix", representation(perm = "integer"),
	 contains = "sparseMatrix",
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

### Class Union :  no inheritance, but is(*, <class>) :

## Definition  Packed := dense with length( . @x) < prod( . @Dim)
##	       ~~~~~~
## REPLACED the following with	isPacked() in ./Auxiliaries.R :
## setClassUnion("packedMatrix",
##		 members = c("dspMatrix", "dppMatrix", "dtpMatrix",
##		  "lspMatrix", "ltpMatrix", "diagonalMatrix"))


## --------------------- non-"Matrix" Classes --------------------------------

## --- "General" (not Matrix at all) ----

## for 'i' in x[i] or A[i,] :
setClassUnion("index", members =  c("numeric", "logical", "character"))


## --- Matrix - related ----

setClass("determinant",
	 representation(modulus = "numeric",
			logarithm = "logical",
			sign = "integer",
			call = "call"))

setClass("LU",
	 representation(x = "numeric", perm = "integer"),
	 validity = function(object) .Call(LU_validate, object)
	 )

## -------------------- lmer-related Classes --------------------------------

setOldClass("data.frame")
setOldClass("family")
setOldClass("logLik")
setOldClass("terms")

## mixed effects representation
setClass("mer",
	 representation(## original data
			flist = "list", # list of grouping factors
			Zt = "dgCMatrix",  # sparse representation of Z'
			X = "matrix",	   # X
			y = "numeric",	   # y
			wts = "numeric",   # weights
			wrkres = "numeric",# working residuals (copy of y for LMMs)
			method = "character", # parameter estimation method
			useScale = "logical", # should scale factor be included
			family = "family", # glm family
			call = "call",	   # call to model-fitting function
			## invariants derived from data structure
			cnames = "list",   # column names of model matrices
			nc = "integer",	   # dimensions of blocks in Omega
			Gp = "integer",	   # Pointers to groups of rows in Zt
			## quantities that vary when Z, X or y are updated
			XtX = "dpoMatrix", # X'X
			ZtZ = "dsCMatrix", # Z'Z
			ZtX = "dgeMatrix", # Z'X
			Zty = "numeric",   # Z'y
			Xty = "numeric",   # X'y
			## primary slots that vary during the optimization
			## When Omega is updated, these are updated
			Omega = "list", # list of relative precision matrices
			## Cholesky factor of inflated [Z:X:y]'[Z:X:y]
			L = "dCHMsuper", # sparse Cholesky factor of Z'Z + Omega
			RZX = "dgeMatrix",
			RXX = "dtrMatrix",
			rZy = "numeric",
			rXy = "numeric",
			devComp = "numeric", # Components of deviance
			deviance = "numeric", # Current deviance (ML and REML)
			## Secondary slots only evaluated when requested.
			fixef = "numeric",
			ranef = "numeric",
			RZXinv = "dgeMatrix",
			bVar = "list",
			gradComp = "list",
			## status indicator
			status = "logical"
			)
	)

## Representation of a linear or generalized linear mixed effects model
setClass("lmer",
	 representation(assign = "integer", frame = "data.frame",
			terms = "terms"),
	 contains = "mer")

setClass("summary.mer", # the "mer" result ``enhanced'' :
	 representation(
			isG   = "logical",
			methTitle = "character",
			logLik= "logLik",
			ngrps = "integer",
			sigma = "numeric", # scale, non-negative number
			coefs = "matrix",
			vcov = "dpoMatrix",
			REmat = "matrix",
			AICtab= "data.frame"
			),
	 contains = "mer")

setClass("summary.lmer", contains = c("summary.mer", "lmer"))

setClass("ranef.lmer", contains = "list")

setClass("coef.lmer", contains = "list")

setClass("pedigree", representation =
	 list(sire = "integer", dam = "integer", label = "character"),
	 validity = function(object) {
	     n <- length(sire <- object@sire)
	     if (length(dam <- object@dam) != n)
		 return("sire and dam slots must be the same length")
	     if (length(object@label) != n)
		 return("'label' slot must have the same length as 'sire' and 'dam'")
	     if(n == 0) return(TRUE)
	     animal <- 1:n
	     snmiss <- !is.na(sire)
	     dnmiss <- !is.na(dam)
	     if (any(sire[snmiss] >= animal[snmiss]) ||
		 any(dam[dnmiss] >= animal[dnmiss]))
		 return("the sire and dam must precede the offspring")
             if (any(sire[snmiss] < 1 | sire[snmiss] > n |
                     dam[dnmiss] < 1 | dam[dnmiss] > n))
                 return(paste("Non-missing sire or dam must be in [1,",
                              n, "]", sep = ''))
	     TRUE
	 })
