#### Define those generics that we need, if they don't exist;
#### not all will be exported

    setGeneric("expand", function(x, ...) standardGeneric("expand"))

## tcrossprod() is now in R's base

    setGeneric("isDiagonal", function(object, ...)
               standardGeneric("isDiagonal"))

## isSymmetric is "S3 generic" in R's  base/R/eigen.R

    setGeneric("isTriangular", function(object, ...) ## 'upper = NA'
               standardGeneric("isTriangular"))

    setGeneric("facmul",
               function(x, factor, y, transpose, left, ...)
               standardGeneric("facmul"))

    setGeneric("lu", function(x, ...) standardGeneric("lu"))

    setGeneric("chol", def = function(x, pivot= FALSE,...) standardGeneric("chol"),
               useAsDefault= function(x, pivot= FALSE,...) base::chol(x, pivot, ...))


    setGeneric("qr", def =   function(x, tol=1e-7,...) standardGeneric("qr"),
               useAsDefault= function(x, tol=1e-7,...) base::qr(x, tol, ...))

    setGeneric("norm", function(x, type, ...) standardGeneric("norm"))

    setGeneric("rcond", function(x, type, ...) standardGeneric("rcond"))

    setGeneric("Schur", function(x, vectors, ...) standardGeneric("Schur"))

    setGeneric("unpack", function(x, ...) standardGeneric("unpack"))

##-     setGeneric("%p%", function(a, b) standardGeneric("%p%"))

    setGeneric("expm", function(x) standardGeneric("expm"))

    setGeneric("writeHB", function(obj, file, ...)
               standardGeneric("writeHB"))

    setGeneric("writeMM", function(obj, file, ...)
               standardGeneric("writeMM"))

##     setGeneric("qqmath", function(x, data, ...)
##                standardGeneric("qqmath"))

    setGeneric("tril", function(x, k = 0, ...)
               standardGeneric("tril"))

    setGeneric("triu", function(x, k = 0, ...)
               standardGeneric("triu"))

    setGeneric("band", function(x, k1, k2, ...)
               standardGeneric("band"))

    setGeneric("Cholesky",
               function(A, perm = TRUE, LDL = TRUE, super = FALSE, ...)
               standardGeneric("Cholesky"))


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

## Add '...' so our methods can add  'sparseResult':
setGeneric("colSums",
	   def = function(x, na.rm = FALSE, dims = 1, ...) standardGeneric("colSums"),
	   useAsDefault = function(x,  na.rm = FALSE, dims = 1, ...)
	   base::colSums(x, na.rm=na.rm, dims=dims))
setGeneric("colMeans",
	   def = function(x, na.rm = FALSE, dims = 1, ...) standardGeneric("colMeans"),
	   useAsDefault = function(x,  na.rm = FALSE, dims = 1, ...)
	   base::colMeans(x, na.rm=na.rm, dims=dims))
setGeneric("rowSums",
	   def = function(x, na.rm = FALSE, dims = 1, ...) standardGeneric("rowSums"),
	   useAsDefault = function(x,  na.rm = FALSE, dims = 1, ...)
	   base::rowSums(x, na.rm=na.rm, dims=dims))
setGeneric("rowMeans",
	   def = function(x, na.rm = FALSE, dims = 1, ...) standardGeneric("rowMeans"),
	   useAsDefault = function(x,  na.rm = FALSE, dims = 1, ...)
	   base::rowMeans(x, na.rm=na.rm, dims=dims))
