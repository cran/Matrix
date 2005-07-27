#### Define those generics that we need, if they don't exist


if (!isGeneric("expand"))
    setGeneric("expand", function(x, ...) standardGeneric("expand"))

if (!isGeneric("tcrossprod"))
    setGeneric("tcrossprod", function(x) standardGeneric("tcrossprod"))

if (!isGeneric("isSymmetric"))
    setGeneric("isSymmetric", function(object, ...)
               standardGeneric("isSymmetric"))

if (!isGeneric("isNested"))
    setGeneric("isNested", function(object, ...) standardGeneric("isNested"))

if (!isGeneric("facmul"))
    setGeneric("facmul",
               function(x, factor, y, transpose, left, ...)
               standardGeneric("facmul"))

if (!isGeneric("lu"))
    setGeneric("lu", function(x, ...) standardGeneric("lu"))

if (!isGeneric("norm"))
    setGeneric("norm", function(x, type, ...) standardGeneric("norm"))

if (!isGeneric("rcond"))
    setGeneric("rcond", function(x, type, ...) standardGeneric("rcond"))

if (!isGeneric("Schur"))
    setGeneric("Schur", function(x, vectors, ...) standardGeneric("Schur"))

if (!isGeneric("unpack"))
    setGeneric("unpack", function(x, ...) standardGeneric("unpack"))

##- if (!isGeneric("%p%"))
##-     setGeneric("%p%", function(a, b) standardGeneric("%p%"))

if (!isGeneric("expm"))
    setGeneric("expm", function(x) standardGeneric("expm"))

## ----------------------- lmer-related Generics ---------------------------

setGeneric("lmer",
           function(formula, data, family,
                    method = c("REML", "ML", "PQL", "Laplace", "AGQ"),
                    control = list(),
                    subset, weights, na.action, offset,
                    model = TRUE, x = FALSE, y = FALSE,...)
           standardGeneric("lmer"))

if (!isGeneric("LMEoptimize<-")) {
    setGeneric("LMEoptimize<-", function(x, ..., value)
               standardGeneric("LMEoptimize<-"))
}

if (!isGeneric("fixef")) {
    setGeneric("fixef", function(object, ...) standardGeneric("fixef"))
}

fixed.effects <- function(object, ...) {
    ## fixed.effects was an alternative name for fixef
    .Deprecated("fixef")
    mCall = match.call()
    mCall[[1]] = as.name("fixef")
    eval(mCall, parent.frame())
}

if (!isGeneric("ranef")) {
    setGeneric("ranef", function(object, ...)
               standardGeneric("ranef"))
}

random.effects <- function(object, ...) {
    ## random.effects was an alternative name for ranef
    .Deprecated("ranef")
    mCall = match.call()
    mCall[[1]] = as.name("ranef")
    eval(mCall, parent.frame())
}

if (!isGeneric("BIC")) {
    setGeneric("BIC", function(object, ...) standardGeneric("BIC"))
}

setMethod("BIC", "logLik",
          function(object, ...)
          -2 * (c(object) - attr(object, "df") * log(attr(object, "nobs"))/2)
          )

if (!isGeneric("VarCorr")) {
    setGeneric("VarCorr", function(x, ...) standardGeneric("VarCorr"))
}

if (!isGeneric("gradient")) {           # not exported
    setGeneric("gradient", function(x, ...) standardGeneric("gradient"))
}

if (!isGeneric("getFixDF")) {           # not exported
    setGeneric("getFixDF", function(object, ...) standardGeneric("getFixDF"))
}

if (!isGeneric("mcmcsamp")) {
    setGeneric("mcmcsamp", function(obj, nsamp = 1, verbose =
    FALSE, ...) standardGeneric("mcmcsamp"))
}
