.onLoad <- function(lib, pkg) {
    if(is.null(getOption("max.print")))
        options(max.print = 10000)#-> show() of large matrices
}

## ------------- Virtual Classes ----------------------------------------

# Mother class of all Matrix objects
setClass("Matrix", representation(Dim = "integer", Dimnames = "list",
                                  "VIRTUAL"),
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
             ## 'else'  ok :
             TRUE
         })

# Virtual class of numeric matrices
setClass("dMatrix",
         representation(x = "numeric", "VIRTUAL"), contains = "Matrix")

# Virtual class of integer matrices
setClass("iMatrix",
         representation(x = "integer", "VIRTUAL"), contains = "Matrix")

# Virtual class of logical matrices
setClass("lMatrix", representation("VIRTUAL"), contains = "Matrix")
## Note that logical sparse matrices do not need an x slot so the x
## slot is part of the ldenseMatrix class

# Virtual class of complex matrices
setClass("zMatrix", # letter 'z' is as in the names of Lapack subroutines
         representation(x = "complex", "VIRTUAL"), contains = "Matrix")

# Virtual class of dense matrices
setClass("denseMatrix", representation("VIRTUAL"), contains = "Matrix")

# Virtual class of dense, numeric matrices
setClass("ddenseMatrix",
         representation(rcond = "numeric", factors = "list", "VIRTUAL"),
         contains = c("dMatrix", "denseMatrix"))

# Virtual class of dense, logical matrices
setClass("ldenseMatrix",
         representation(x = "logical", factors = "list", "VIRTUAL"),
         contains = c("lMatrix", "denseMatrix"))

## virtual SPARSE ------------

setClass("sparseMatrix", representation("VIRTUAL"), contains = "Matrix")

## general Triplet Matrices (dgT, lgT, ..):
setClass("gTMatrix", representation(i = "integer", j = "integer", "VIRTUAL"),
         contains = "sparseMatrix")

setClass("dsparseMatrix", representation("VIRTUAL"),
         contains = c("dMatrix", "sparseMatrix"))

setClass("lsparseMatrix", representation("VIRTUAL"),
         contains = c("lMatrix", "sparseMatrix"))

## ------------------ Proper (non-virtual) Classes ----------------------------

##----------------------  DENSE  -----------------------------------------

# numeric, dense, general matrices
setClass("dgeMatrix", contains = "ddenseMatrix",
         ## checks that length( @ x) == prod( @ Dim):
         validity = function(object) .Call("dgeMatrix_validate", object)
         )
## i.e. "dgeMatrix" cannot be packed, but "ddenseMatrix" can ..

# numeric, dense, non-packed, triangular matrices
setClass("dtrMatrix",
         representation(uplo = "character", diag = "character"),
         contains = "dgeMatrix",
         prototype = prototype(uplo = "U", diag = "N"),
         validity = function(object) .Call("dtrMatrix_validate", object)
         )

# numeric, dense, packed, triangular matrices
setClass("dtpMatrix",
         representation(uplo = "character", diag = "character"),
         contains = "ddenseMatrix",
         prototype = prototype(uplo = "U", diag = "N"),
         validity = function(object) .Call("dtpMatrix_validate", object)
         )

# numeric, dense, non-packed symmetric matrices
setClass("dsyMatrix",
         representation(uplo = "character"),
         prototype = prototype(uplo = "U"),
         contains = "dgeMatrix",
         validity = function(object) .Call("dsyMatrix_validate", object)
         )

# numeric, dense, packed symmetric matrices
setClass("dspMatrix",
         representation(uplo = "character"),
         prototype = prototype(uplo = "U"),
         contains = "ddenseMatrix",
         validity = function(object) .Call("dspMatrix_validate", object)
         )

# numeric, dense, non-packed, positive-definite, symmetric matrices
setClass("dpoMatrix", contains = "dsyMatrix",
         validity = function(object) .Call("dpoMatrix_validate", object)
         )

# numeric, dense, packed, positive-definite, symmetric matrices
setClass("dppMatrix", contains = "dspMatrix",
         validity = function(object) .Call("dppMatrix_validate", object)
         )

##-------------------- S P A R S E (non-virtual) --------------------------

##---------- numeric sparse matrix classes --------------------------------

# numeric, sparse, triplet general matrices
setClass("dgTMatrix",
         representation(factors = "list"),
         contains = c("gTMatrix", "dsparseMatrix"),
         validity = function(object) .Call("dgTMatrix_validate", object)
         )

# numeric, sparse, triplet triangular matrices
setClass("dtTMatrix",
         representation(uplo = "character", diag = "character"),
         contains = "dgTMatrix",
         validity = function(object) .Call("dtTMatrix_validate", object)
         )

# numeric, sparse, triplet symmetric matrices
setClass("dsTMatrix",
         representation(uplo = "character"),
         contains = "dgTMatrix",
         validity = function(object) .Call("dsTMatrix_validate", object)
         )

# numeric, sparse, sorted compressed sparse column-oriented general matrices
setClass("dgCMatrix",
         representation(i = "integer", p = "integer", factors = "list"),
         contains = "dsparseMatrix",
         validity = function(object) .Call("dgCMatrix_validate", object)
         )

# numeric, sparse, sorted compressed sparse column-oriented triangular matrices
setClass("dtCMatrix",
         representation(uplo = "character", diag = "character"),
         contains = "dgCMatrix",
         validity = function(object) .Call("tsc_validate", object)
         )

# numeric, sparse, sorted compressed sparse column-oriented symmetric matrices
setClass("dsCMatrix",
         representation(uplo = "character"),
         contains = "dgCMatrix",
         validity = function(object) .Call("dsCMatrix_validate", object)
         )

# numeric, sparse, sorted compressed sparse row-oriented general matrices
setClass("dgRMatrix",
         representation(j = "integer", p = "integer", factors = "list"),
         contains = "dsparseMatrix",
         ##TODO: validity = function(object) .Call("dgRMatrix_validate", object)
         )

# numeric, sparse, sorted compressed sparse row-oriented triangular matrices
setClass("dtRMatrix",
         representation(uplo = "character", diag = "character"),
         contains = "dgRMatrix",
         ##TODO: validity = function(object) .Call("dtRMatrix_validate", object)
         )

# numeric, sparse, sorted compressed sparse row-oriented symmetric matrices
setClass("dsRMatrix",
         representation(uplo = "character"),
         contains = "dgRMatrix",
         ##TODO: validity = function(object) .Call("dsRMatrix_validate", object)
         )

##---------- logical sparse matrix classes --------------------------------

## these classes are used in symbolic analysis to determine the
## locations of non-zero entries

# logical, sparse, triplet general matrices
setClass("lgTMatrix",
         contains = c("gTMatrix", "lsparseMatrix"),
         validity = function(object) .Call("lgTMatrix_validate", object)
         )

# logical, sparse, triplet triangular matrices
setClass("ltTMatrix",
         representation(uplo = "character", diag = "character"),
         contains = "lgTMatrix",
         validity = function(object) .Call("ltTMatrix_validate", object)
         )

# logical, sparse, triplet symmetric matrices
setClass("lsTMatrix",
         representation(uplo = "character"),
         contains = "lgTMatrix",
         validity = function(object) .Call("lsTMatrix_validate", object)
         )

# logical, sparse, sorted compressed sparse column-oriented general matrices
setClass("lgCMatrix",
         representation(i = "integer", p = "integer"),
         contains = "lsparseMatrix",
         validity = function(object) .Call("lgCMatrix_validate", object)
         )

# logical, sparse, sorted compressed sparse column-oriented triangular matrices
setClass("ltCMatrix",
         representation(uplo = "character", diag = "character"),
         contains = "lgCMatrix",
         validity = function(object) .Call("ltCMatrix_validate", object)
         )

# logical, sparse, sorted compressed sparse column-oriented symmetric matrices
setClass("lsCMatrix",
         representation(uplo = "character"),
         contains = "lgCMatrix",
         validity = function(object) .Call("lsCMatrix_validate", object)
         )

# logical, sparse, sorted compressed sparse row-oriented general matrices
setClass("lgRMatrix",
         representation(j = "integer", p = "integer"),
         contains = "lsparseMatrix",
         validity = function(object) .Call("lgRMatrix_validate", object)
         )

# logical, sparse, sorted compressed sparse row-oriented triangular matrices
setClass("ltRMatrix",
         representation(uplo = "character", diag = "character"),
         contains = "lgRMatrix",
         validity = function(object) .Call("ltRMatrix_validate", object)
         )

# logical, sparse, sorted compressed sparse row-oriented symmetric matrices
setClass("lsRMatrix",
         representation(uplo = "character"),
         contains = "lgRMatrix",
         validity = function(object) .Call("lsRMatrix_validate", object)
         )

## Compressed sparse column matrix in blocks

setClass("dgBCMatrix",
         representation(p = "integer", i = "integer", x = "array"))

## Factorization classes

setClass("Cholesky", contains = "dtrMatrix")

setClass("pCholesky", contains = "dtpMatrix")

setClass("BunchKaufman", representation(perm = "integer"), contains = "dtrMatrix",
         validity = function(object) .Call("BunchKaufman_validate", object));

setClass("pBunchKaufman", representation(perm = "integer"), contains = "dtpMatrix",
         validity = function(object) .Call("pBunchKaufman_validate", object));

setClass("dCholCMatrix",
         representation(perm = "integer", Parent = "integer", D = "numeric"),
         contains = "dtCMatrix",
         validity = function(object) .Call("dCholCMatrix_validate", object))

setClass("lCholCMatrix",
         representation(perm = "integer", Parent = "integer"),
         contains = "ltCMatrix",
         validity = function(object) .Call("lCholCMatrix_validate", object))

##-------------------- permutation ----------------------------------------

setClass("pMatrix", representation(perm = "integer"), contains = "Matrix",
         validity = function(object) {
             dd <- object@Dim
             n <- dd[1]
             perm <- object@perm
             if (dd[2] != n) return("pMatrix must be symmetric")
             if (length(perm) != n)
                 return(paste("length of 'perm' slot must be", n))
             if (!(all(range(perm) == c(1, n)) && length(unique(perm)) == n))
                 return("'perm' slot is not a valid permutation")
             TRUE
         })

## --------------------- non-"Matrix" Classes --------------------------------

setClass("determinant",
         representation(modulus ="numeric",
                        logarithm = "logical",
                        sign = "integer",
                        call = "call"))

setClass("LU", representation(x = "numeric",
                              perm = "integer"),
         validity = function(object) .Call("LU_validate", object))

## Deprecated:
                       # positive-definite symmetric matrices as matrices
setClass("pdmatrix", contains="matrix")

#                        # factors of positive-definite symmetric matrices
# setClass("pdfactor", representation("matrix", logDet = "numeric"))

                       # correlation matrices and standard deviations
setClass("corrmatrix", representation("matrix", stdDev = "numeric"))

## -------------------- lmer-related Classes --------------------------------

setOldClass("data.frame")
setOldClass("family")
setOldClass("logLik")
setOldClass("terms")

setClass("VarCorr",
         representation(scale="numeric",
                        reSumry="list",
                        useScale="logical"),
         prototype = list(scale = 1.0, useScale = TRUE))

## mixed effects representation
setClass("mer",
         representation(
                        flist = "list", # list of grouping factors
                        perm = "list",  # list of permutations of levels (0-based)
                        Parent = "list",# list of Parent arrays for ZZpO
                        D = "list",     # list of diagonal factors (upper triangle)
                        bVar = "list",  # list of conditional variance factors (upper triangle)
                        L = "list",     # list of blocks of L
                        ZZpO = "list",  # list of diagonal blocks of Z'Z+Omega
                        Omega = "list", # list of relative precision matrices
                        method = "character", # parameter estimation method
                        RXX = "matrix", # Augmented RXX component or its inverse
                        RZX = "matrix", # Augmented RZX component or its inverse
                        XtX = "matrix", # Original X'X matrix
                        ZtZ = "list",   # list of blocks of Z'Z
                        ZtX = "matrix", # Original Z'X matrix
                        cnames = "list",# column names of model matrices
                        devComp = "numeric", # Components of deviance
                        deviance = "numeric", # Current deviance (ML and REML)
                        nc = "integer", # number of columns in (augmented)
                                        # model matrices and number of observations
                        Gp = "integer", # Pointers to groups of rows in RZX
                        status = "logical"
                        ),
         validity = function(object) {
             .Call("lmer_validate", object, PACKAGE = "Matrix")
         })

## Representation of a linear or generalized linear mixed effects model
setClass("lmer",
         representation(assign = "integer", call = "call",
                        family = "family", fitted = "numeric",
                        fixed = "numeric", frame = "data.frame",
                        logLik = "logLik", residuals = "numeric",
                        terms = "terms"),
         contains = "mer")

## Representation of a generalized linear mixed effects model
##setClass("glmer",
##         representation(family = "family", glmmll = "numeric", fixed = "numeric"),
##         contains = "lmer")

setClass("summary.lmer",
         representation(useScale="logical",
                        showCorrelation="logical"),
         contains = "lmer")

setClass("lmer.ranef",
         representation(varFac = "list", stdErr = "numeric"),
         contains = "list")

setClass("lmer.ranef.confint", contains = "list")

setClass("lmer.coef",
         representation(varFac = "list", stdErr = "numeric"),
         contains = "list")
