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
setClass("lMatrix",
         representation(x = "logical", "VIRTUAL"), contains = "Matrix")

# Virtual class of complex matrices
setClass("zMatrix", # letter 'z' is as in the names of Lapack subroutines
         representation(x = "complex", "VIRTUAL"), contains = "Matrix")

# Virtual class of dense, numeric matrices
setClass("ddenseMatrix",
         representation(rcond = "numeric", factors = "list", "VIRTUAL"),
         contains = "dMatrix")

## virtual SPARSE ------------

setClass("sparseMatrix", contains = "Matrix")# "VIRTUAL"

setClass("dsparseMatrix", contains = c("dMatrix", "sparseMatrix"))# "VIRTUAL"

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

# numeric, sparse, triplet general matrices
setClass("dgTMatrix",
         representation(i = "integer", j = "integer", factors = "list"),
         contains = "dsparseMatrix",
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
         validity = function(object) .Call("dgRMatrix_validate", object)
         )

# numeric, sparse, sorted compressed sparse row-oriented triangular matrices
setClass("dtRMatrix",
         representation(uplo = "character", diag = "character"),
         contains = "dgRMatrix",
         validity = function(object) .Call("dtRMatrix_validate", object)
         )

# numeric, sparse, sorted compressed sparse row-oriented symmetric matrices
setClass("dsRMatrix",
         representation(uplo = "character"),
         contains = "dgRMatrix",
         validity = function(object) .Call("dsRMatrix_validate", object)
         )

## Compressed sparse column matrix in blocks
setClass("dgBCMatrix",
         representation(p = "integer", i = "integer", x = "array"))

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
