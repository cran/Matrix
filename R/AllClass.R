# Ensure that the methods package is available
.onLoad <- function(lib, pkg) {
    require("methods", character = TRUE, quietly = TRUE)
}

# Virtual class of all Matrix objects
setClass("Matrix")

# general double-precision matrices
setClass("geMatrix",
         representation(x = "numeric", Dim = "integer", 
                        rcond = "numeric", factorization = "list"),
         prototype = prototype(x = numeric(0), Dim = as.integer(c(0,0)),
                               rcond = numeric(0),
                               factorization = list()),
         contains = "Matrix",
         validity = function(object) {
             .Call("geMatrix_validate", object)
         })

# Dense, non-packed, triangular matrices 
setClass("trMatrix",
         representation(uplo = "character", diag = "character"),
         contains = "geMatrix",
         prototype = prototype(uplo = "U", diag = "N"),
         validity = function(object) {
             .Call("trMatrix_validate", object)
         })

# Dense, packed, triangular matrices
# setClass("tpMatrix", representation(...))

# Dense, non-packed symmetric matrices
setClass("syMatrix",
         representation(uplo = "character"),
         prototype = prototype(uplo = "U",
           x = numeric(0), Dim = as.integer(c(0,0)),
           rcond = numeric(0), factorization = list()),
         contains = "geMatrix",
         validity = function(object) {
             .Call("syMatrix_validate", object)
         })

# Dense, packed, symmetric matrices
# setClass("spMatrix", representation(...))

# Dense, non-packed, postive-definite, symmetric matrices
setClass("poMatrix", contains = "syMatrix",
         prototype = prototype(uplo = "U",
           x = numeric(0), Dim = as.integer(c(0,0)),
           rcond = numeric(0), factorization = list())
         )

# Sparse general matrix in sorted compressed sparse column format
setClass("cscMatrix",
         representation(p = "integer", i = "integer", x = "numeric",
                        Dim = "integer", factorization = "list"),
         prototype = prototype(p = as.integer(0), i = integer(0),
                        x = numeric(0), Dim = as.integer(c(0, 0))),
         validity = function(object)
                    .Call("csc_validate", object)
         )

# Sparse triangular matrix in sorted compressed sparse column format
setClass("tscMatrix",
         representation(uplo = "character", diag = "character"),
         prototype = prototype(p = as.integer(0), i = integer(0),
                        x = numeric(0), Dim = as.integer(c(0, 0)),
                        uplo = 'L', diag = 'N'),
         contains = "cscMatrix",
         validity = function(object)
                    .Call("tsc_validate", object))

# Sparse symmetric matrix in compressed sparse column format.
# Only one triangle is stored, uplo indicates if it is the lower or upper
setClass("sscMatrix",
         representation = representation(uplo = "character"),
         prototype = prototype(p = as.integer(0), i = integer(0),
                        x = numeric(0), Dim = as.integer(c(0, 0)),
                        uplo = 'L'),
         contains = "cscMatrix",
         validity = function(object)
                    .Call("sscMatrix_validate", object))

# Sparse general matrix in triplet format
setClass("tripletMatrix",
         representation(i = "integer", j = "integer", x = "numeric",
                        Dim = "integer"),
         prototype = prototype(i = integer(0), j = integer(0),
         x = numeric(0), Dim = as.integer(c(0,0))),
         validity = function(object)
                    .Call("triplet_validate", object))

setClass("determinant",
         representation(modulus ="numeric",
                        logarithm = "logical",
                        sign = "integer",
                        call = "call"))

setClass("LU", representation(x = "numeric", Dim = "integer",
                              pivot = "integer"),
         validity = function(object)
                    .Call("LU_validate", object))

setClass("Cholesky", contains = "trMatrix")

setClass("sscChol",
         representation = representation(perm = "integer", Parent = "integer",
         D = "numeric"),
         contains = "tscMatrix",
         prototype = prototype(p = as.integer(0), i = integer(0),
                        x = numeric(0), Dim = as.integer(c(0, 0)),
                        uplo = 'L', perm = integer(0), Parent = integer(0),
                        D = numeric(0)),
         validity = function(object)
           .Call("sscChol_validate", object))

setClass("sscCrosstab", representation =
         representation(Gp = "integer", perm = "integer"),
         contains = "sscMatrix",
         validity = function(object)
           .Call("sscCrosstab_validate", object))

setClass("ssclme", representation =
         representation(
                        D = "numeric",  # Diagonal of D in LDL'
                        DIsqrt = "numeric", # inverse square root of D
                        Dim = "integer", # Dimensions of Z'Z and LDL'
                        Gp = "integer", # Pointers to groups of columns of Z
                        Li = "integer", # Row indices of L
                        Lp = "integer", # Column pointers of L
                        Lx = "numeric", # Non-zero, off-diagonals of L
                        Omega = "list", # List of symmetric matrices
                        Parent = "integer", # Elimination tree of L
                        RXX = "matrix", # Augmented RXX component or inverse
                        RZX = "matrix", # Augmented RZX component or inverse
                        XtX = "matrix", # Original X'X matrix
                        ZtX = "matrix", # Original Z'X matrix
                        bVar = "list",  # Diagonal blocks on (Z'Z+W)^{-1}
                        deviance = "numeric", # Current deviance (ML and REML)
                        devComp = "numeric", # Components of deviance
                        i = "integer",  # Row indices of Z'Z
                        nc = "integer", # number of columns in model matrices
                        p = "integer",  # Pointers to columns of Z'Z
                        status = "logical", # record if factored, if inverted
                        x = "numeric"   # Non-zeroes in upper triangle of Z'Z
                        ),
         validity = function(object)
         .Call("ssclme_validate", object))

                       # positive-definite symmetric matrices as matrices
setClass("pdmatrix", contains="matrix")

                       # factors of positive-definite symmetric matrices
setClass("pdfactor", representation("matrix", logDet = "numeric"))

                       # correlation matrices and standard deviations
setClass("corrmatrix", representation("matrix", stdDev = "numeric"))

## Compressed sparse column matrix in blocks
setClass("cscBlocked", representation(p = "integer", i = "integer", x = "array"))

## Block/block sparse symmetric matrices
setClass("bbSparseSy", representation(x = "list", uplo = "character"))

## Block/block cross tabulation
setClass("bbCrosstab", contains = "bbSparseSy")

## Block/block sparse triangular matrices
setClass("bbSparseTr", representation(x = "list", uplo = "character",
                                      diag = "character"))

## Block/block L matrix
setClass("bbLmat", representation(Linv = "list"), contains = "bbSparseTr")
         
## Representation of a linear mixed effects model
setClass("lmeRep",
         representation(
                        D = "list",     # list of diagonal factors (lower triangle)
                        L = "list",     # list of blocks of L
                        Linv = "list",  # list of blocks of L^{-1}
                        Omega = "list", # list of relative precision matrices
                        RXX = "matrix", # Augmented RXX component or its inverse
                        RZX = "matrix", # Augmented RZX component or its inverse
                        XtX = "matrix", # Original X'X matrix
                        ZZx = "list",   # list of blocks of Z'Z
                        ZtX = "matrix", # Original Z'X matrix
                        cnames = "list",# column names of model matrices
                        devComp = "numeric", # Components of deviance
                        deviance = "numeric", # Current deviance (ML and REML)
                        levels = "list",# names of levels of grouping factors
                        nc = "integer", # number of columns in (augmented)
                                        # model matrices and number of observations
                        status = "logical",
                        call = "call"   # omit this after debugging phase
                        ),
         validity = function(object)
         .Call("lmeRep_validate", object))

