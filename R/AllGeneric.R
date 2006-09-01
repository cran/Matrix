#### Define those generics that we need, if they don't exist;
#### not all will be exported

if (!isGeneric("expand"))
    setGeneric("expand", function(x, ...) standardGeneric("expand"))

## tcrossprod() is now in R's base

if (!isGeneric("isDiagonal"))
    setGeneric("isDiagonal", function(object, ...)
               standardGeneric("isDiagonal"))

## isSymmetric is "S3 generic" in R's  base/R/eigen.R

if (!isGeneric("isTriangular"))
    setGeneric("isTriangular", function(object, ...) ## 'upper = NA'
               standardGeneric("isTriangular"))

if (!isGeneric("facmul"))
    setGeneric("facmul",
               function(x, factor, y, transpose, left, ...)
               standardGeneric("facmul"))

if (!isGeneric("lu"))
    setGeneric("lu", function(x, ...) standardGeneric("lu"))

if (!isGeneric("chol"))
    setGeneric("chol", def = function(x, pivot= FALSE,...) standardGeneric("chol"),
               useAsDefault= function(x, pivot= FALSE,...) base::chol(x, pivot, ...))


if (!isGeneric("qr"))
    setGeneric("qr", def =   function(x, tol=1e-7,...) standardGeneric("qr"),
               useAsDefault= function(x, tol=1e-7,...) base::qr(x, tol, ...))

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

if (!isGeneric("writeHB"))
    setGeneric("writeHB", function(obj, file, ...)
               standardGeneric("writeHB"))

if (!isGeneric("writeMM"))
    setGeneric("writeMM", function(obj, file, ...)
               standardGeneric("writeMM"))

## if (!isGeneric("qqmath"))
##     setGeneric("qqmath", function(x, data, ...)
##                standardGeneric("qqmath"))

if (!isGeneric("tril"))
    setGeneric("tril", function(x, k = 0, ...)
               standardGeneric("tril"))

if (!isGeneric("triu"))
    setGeneric("triu", function(x, k = 0, ...)
               standardGeneric("triu"))

if (!isGeneric("band"))
    setGeneric("band", function(x, k1, k2, ...)
               standardGeneric("band"))

if (!isGeneric("Cholesky"))
    setGeneric("Cholesky",
               function(A, perm = TRUE, LDL = TRUE, super = FALSE, ...)
               standardGeneric("Cholesky"))


## ----------------------- lmer-related Generics ---------------------------

## Hmm: If this does not match *exactly* the "formula" - method in ./lmer.R
## ---  the  match.call() in there may give a very different result
setGeneric("lmer",
	   function(formula, data, family = gaussian,
		    method = c("REML", "ML", "PQL", "Laplace", "AGQ"),
		    control = list(), start, subset, weights, na.action,
		    offset, contrasts = NULL, model = TRUE,
		    ...)
	   standardGeneric("lmer"))

if (!isGeneric("isNested"))
    setGeneric("isNested", function(x, ...) standardGeneric("isNested"))

if (!isGeneric("LMEoptimize<-")) {
    setGeneric("LMEoptimize<-", function(x, ..., value)
               standardGeneric("LMEoptimize<-"))
}

if (!isGeneric("fixef")) {
    setGeneric("fixef", function(object, ...) standardGeneric("fixef"))
}

if (!isGeneric("denomDF")) {
    setGeneric("denomDF", function(x, ...) standardGeneric("denomDF"))
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

if (!isGeneric("postVar")) {            # posterior variances
    setGeneric("postVar", function(object, ...)
               standardGeneric("postVar"))
}

if (!isGeneric("gradient")) {           # not exported
    setGeneric("gradient", function(x, ...) standardGeneric("gradient"))
}

if (!isGeneric("getFixDF")) {           # not exported
    setGeneric("getFixDF", function(object, ...) standardGeneric("getFixDF"))
}

if (!isGeneric("mcmcsamp")) {
    setGeneric("mcmcsamp", function(object, n = 1, verbose = FALSE, ...)
	       standardGeneric("mcmcsamp"))
}

if (!exists("simulate", mode = "function")) {
    setGeneric("simulate", function(object, nsim = 1, seed = NULL, ...)
               standardGeneric("simulate"))
}

###---- Group Generics ----
## The following are **WORKAROUND** s currently needed for all non-Primitives:

##  "Math"
setGeneric("log", group="Math")
setGeneric("gamma", group="Math")
setGeneric("lgamma", group="Math")

## "Math2"
setGeneric("round",  group="Math2")
setGeneric("signif", group="Math2")

## "Summary" --- this needs some hoop jumping that may become unnecessary
##               in a future version of R (>= 2.3.x):

.max_def <- function(x, ..., na.rm = FALSE) base::max(x, ..., na.rm = na.rm)
.min_def <- function(x, ..., na.rm = FALSE) base::min(x, ..., na.rm = na.rm)
.range_def <- function(x, ..., na.rm = FALSE) base::range(x, ..., na.rm = na.rm)
.prod_def <- function(x, ..., na.rm = FALSE) base::prod(x, ..., na.rm = na.rm)
.sum_def <- function(x, ..., na.rm = FALSE) base::sum(x, ..., na.rm = na.rm)
.any_def <- function(x, ..., na.rm = FALSE) base::any(x, ..., na.rm = na.rm)
.all_def <- function(x, ..., na.rm = FALSE) base::all(x, ..., na.rm = na.rm)

setGeneric("max", function(x, ..., na.rm = FALSE) standardGeneric("max"),
           useAsDefault = .max_def, group = "Summary")
setGeneric("min", function(x, ..., na.rm = FALSE) standardGeneric("min"),
           useAsDefault = .min_def, group="Summary")
setGeneric("range", function(x, ..., na.rm = FALSE) standardGeneric("range"),
           useAsDefault = .range_def, group="Summary")
setGeneric("prod", function(x, ..., na.rm = FALSE) standardGeneric("prod"),
           useAsDefault = .prod_def, group="Summary")
setGeneric("sum", function(x, ..., na.rm = FALSE) standardGeneric("sum"),
           useAsDefault = .sum_def, group="Summary")
setGeneric("any", function(x, ..., na.rm = FALSE) standardGeneric("any"),
           useAsDefault = .any_def, group="Summary")
setGeneric("all", function(x, ..., na.rm = FALSE) standardGeneric("all"),
           useAsDefault = .all_def, group="Summary")
