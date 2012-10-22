#### Define those generics that we need, if they don't exist;
#### not all will be exported

    setGeneric("expand", function(x, ...) standardGeneric("expand"))

## these two are *not* exported (yet) :

    setGeneric("isDiagonal", function(object, ...)
               standardGeneric("isDiagonal"))

    setGeneric("isTriangular", function(object, ...) ## 'upper = NA'
               standardGeneric("isTriangular"))

## isSymmetric is "S3 generic" in R's  base/R/eigen.R

    setGeneric("facmul",
               function(x, factor, y, transpose, left, ...)
               standardGeneric("facmul"))

setGeneric("BunchKaufman", function(x, ...) standardGeneric("BunchKaufman"))
setGeneric("lu", function(x, ...) standardGeneric("lu"))

##NB ## do not redefine the "base signature"

##NB setGeneric("chol", def = function(x, pivot= FALSE,...) standardGeneric("chol"),
##NB            useAsDefault= function(x, pivot= FALSE,...) base::chol(x, pivot, ...))

##NB     setGeneric("qr", def =   function(x, tol=1e-7,...) standardGeneric("qr"),
##NB                useAsDefault= function(x, tol=1e-7,...) base::qr(x, tol, ...))

    setGeneric("Schur", function(x, vectors, ...) standardGeneric("Schur"))

setGeneric("unpack", function(x, ...) standardGeneric("unpack"))
setGeneric("pack", function(x, ...) standardGeneric("pack"))

##-     setGeneric("%p%", function(a, b) standardGeneric("%p%"))

    setGeneric("expm", function(x) standardGeneric("expm"))

    setGeneric("writeMM", function(obj, file, ...)
               standardGeneric("writeMM"))

    setGeneric("tril", function(x, k = 0, ...)
               standardGeneric("tril"))

    setGeneric("triu", function(x, k = 0, ...)
               standardGeneric("triu"))

    setGeneric("band", function(x, k1, k2, ...)
               standardGeneric("band"))

    setGeneric("Cholesky",
	       function(A, perm = TRUE, LDL = !super, super = FALSE,
			Imult = 0, ...)
               standardGeneric("Cholesky"))

setGeneric("symmpart", function(x) standardGeneric("symmpart"))
setGeneric("skewpart", function(x) standardGeneric("skewpart"))

## A version of coercion to  "symmetric" which does *NOT* check,
## but just takes the ## upper (or lower) values and
## ``declares'' the symmetric:
setGeneric("forceSymmetric",
	   function(x, uplo) standardGeneric("forceSymmetric"))

setGeneric("nnzero", function(x, na.counted = NA) standardGeneric("nnzero"),
	   signature = "x")

setGeneric("updown", function(update, C, L) standardGeneric("updown"))

if(as.numeric(R.version$`svn rev`) < 60620)
setGeneric("toeplitz", function(x, ...) standardGeneric("toeplitz"),
           useAsDefault= function(x, ...) stats::toeplitz(x))
## and an entry in ../man/sparseVector-class.Rd
