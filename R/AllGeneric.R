# Define those generics that we need, if they don't exist

if (!isGeneric("expand")) {
    setGeneric("expand", function(x, ...) standardGeneric("expand"))
}

if (!isGeneric("tcrossprod")) {
    setGeneric("tcrossprod", function(x) standardGeneric("tcrossprod"))
}

if (!isGeneric("isSymmetric")) {
    setGeneric("isSymmetric", function(object, ...) standardGeneric("isSymmetric"))
}

if (!isGeneric("isNested")) {
    setGeneric("isNested", function(object, ...) standardGeneric("isNested"))
}

if (!isGeneric("facmul")) {
    setGeneric("facmul",
               function(x, factor, y, transpose, left, ...)
               standardGeneric("facmul"))
}

if (!isGeneric("lu")) {
    setGeneric("lu", function(x, ...) standardGeneric("lu"))
}

if (!isGeneric("norm")) {
    setGeneric("norm", function(x, type, ...) standardGeneric("norm"))
}

if (!isGeneric("rcond")) {
    setGeneric("rcond", function(x, type, ...) standardGeneric("rcond"))
}

if (!isGeneric("Schur")) {
    setGeneric("Schur", function(x, vectors, ...) standardGeneric("Schur"))
}

if (!isGeneric("unpack")) {
    setGeneric("unpack", function(x, ...) standardGeneric("unpack"))
}

if (!isGeneric("%p%")) {
    setGeneric("%p%", function(a, b) standardGeneric("%p%"))
}

if (!isGeneric("coef<-")) {
    setGeneric("coef<-", function(object, ..., value)
               standardGeneric("coef<-"))
}

## The generics pdFactor, pdMatrix, corFactor and corMatrix will be deprecated

if (!isGeneric("pdFactor")) {
    setGeneric("pdFactor", function(object) standardGeneric("pdFactor"))
}

if (!isGeneric("pdMatrix")) {
    setGeneric("pdMatrix", function(object) standardGeneric("pdMatrix"))
}

if (!isGeneric("corFactor")) {
    ## extractor for transpose inverse square root factor of corr matrix
    setGeneric("corFactor", function(object, ...) standardGeneric("corFactor"))
}

if (!isGeneric("corMatrix")) {
    ## extractor for correlation matrix or the transpose inverse
    ## square root matrix
    setGeneric("corMatrix", function(object, ...) standardGeneric("corMatrix"))
}

# if (!isGeneric("isInitialized")) {
#     setGeneric("isInitialized",
#                function(object) standardGeneric("isInitialized"),
#                valueClass = "logical")
# }

if (!isGeneric("matrix<-")) {
    setGeneric("matrix<-",
               function(object, value) standardGeneric("matrix<-"))
}
