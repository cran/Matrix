# Ensure that the methods package is available, initialize symbols
.onLoad <- function(lib, pkg) {
    require("methods", character = TRUE, quietly = TRUE)
    .Call("Matrix_init", PACKAGE = "Matrix")
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
             .Call("geMatrix_validate", object, PACKAGE="Matrix")
         })

# Dense, non-packed, triangular matrices 
setClass("trMatrix",
         representation(uplo = "character", diag = "character"),
         contains = "geMatrix",
         prototype = prototype(uplo = "U", diag = "N"),
         validity = function(object) {
             .Call("trMatrix_validate", object, PACKAGE="Matrix")
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
             .Call("syMatrix_validate", object, PACKAGE="Matrix")
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
                    .Call("csc_validate", object, PACKAGE="Matrix")
         )

# Sparse triangular matrix in sorted compressed sparse column format
setClass("tscMatrix",
         representation(uplo = "character", diag = "character"),
         prototype = prototype(p = as.integer(0), i = integer(0),
                        x = numeric(0), Dim = as.integer(c(0, 0)),
                        uplo = 'L', diag = 'N'),
         contains = "cscMatrix",
         validity = function(object)
                    .Call("tsc_validate", object, PACKAGE="Matrix"))

# Sparse symmetric matrix in compressed sparse column format.
# Only one triangle is stored, uplo indicates if it is the lower or upper
setClass("sscMatrix",
         representation = representation(uplo = "character"),
         prototype = prototype(p = as.integer(0), i = integer(0),
                        x = numeric(0), Dim = as.integer(c(0, 0)),
                        uplo = 'L'),
         contains = "cscMatrix",
         validity = function(object)
                    .Call("sscMatrix_validate", object, PACKAGE="Matrix"))

# Sparse general matrix in triplet format
setClass("tripletMatrix",
         representation(i = "integer", j = "integer", x = "numeric",
                        Dim = "integer"),
         prototype = prototype(i = integer(0), j = integer(0),
         x = numeric(0), Dim = as.integer(c(0,0))),
         validity = function(object)
                    .Call("triplet_validate", object, PACKAGE="Matrix"))

setClass("determinant",
         representation(modulus ="numeric",
                        logarithm = "logical",
                        sign = "integer",
                        call = "call"))

setClass("LU", representation(x = "numeric", Dim = "integer",
                              pivot = "integer"),
         validity = function(object)
                    .Call("LU_validate", object, PACKAGE = "Matrix"))

setClass("Cholesky", contains = "trMatrix")

setClass("sscChol",
         representation = representation(perm = "integer", iperm = "integer"),
         contains = "tscMatrix",
         prototype = prototype(p = as.integer(0), i = integer(0),
                        x = numeric(0), Dim = as.integer(c(0, 0)),
                        uplo = 'L', perm = integer(0), iperm = integer(0)),
         validity = function(object)
           .Call("sscChol_validate", object, PACKAGE = "Matrix"))

setClass("sscCrosstab", representation =
         representation(Gp = "integer", perm = "integer"),
         contains = "sscMatrix",
         validity = function(object)
           .Call("sscCrosstab_validate", object, PACKAGE = "Matrix"))

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
         .Call("ssclme_validate", object, PACKAGE = "Matrix"))

setClass("pdMat",      # parameterized positive-definite matrices
         representation(form="formula",    # a model-matrix formula
                        Names="character", # column (and row) names
                        param="numeric",   # parameter vector
                        Ncol="integer",    # number of columns
                        factor="matrix",   # factor of the pos-def matrix
                        logDet="numeric"   # logarithm of the absolute value
                        ## of the determinant of the factor (i.e. half
                        ## the logarithm of the determinant of the matrix)
                        ),
         prototype(form=formula(NULL),
                   Names=character(0),
                   param=numeric(0),
                   Ncol=as.integer(0),
                   factor=matrix(numeric(0),0,0),
                   logDet=numeric(0))
         )

#setClass("pdSymm", contains="pdMat")    # general symmetric pd matrices

#setClass("pdScalar", contains="pdMat") # special case of positive scalars
setClass("pdLogChol", contains="pdMat") # default parameterization
setClass("pdNatural", contains="pdMat") # log sd and logistic of correlation
#setClass("pdMatrixLog", contains="pdSymm") # matrix logarithm parameterization

setClass("pdDiag", contains="pdMat")    # diagonal pd matrices

setClass("pdIdent", contains="pdMat")   # positive multiple of the identity

setClass("pdCompSymm", contains="pdMat") # compound symmetric pd matrices

#setClass("pdBlocked",                   # block-diagonal pd matrices
#         representation("pdMat", components = "list"))

                       # positive-definite symmetric matrices as matrices
setClass("pdmatrix", contains="matrix")

                       # factors of positive-definite symmetric matrices
setClass("pdfactor", representation("matrix", logDet = "numeric"))

                       # correlation matrices and standard deviations
setClass("corrmatrix", representation("matrix", stdDev = "numeric"))

