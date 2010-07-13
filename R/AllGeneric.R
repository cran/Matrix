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

if(getRversion() < "2.11.0" || R.version$`svn rev` < 51018)
    setGeneric("norm", function(x, type, ...) standardGeneric("norm"))
## it is *implicit* generic in standard R from 2.11.0 (~= rev 51018)

if(getRversion() < "2.10.0" || R.version$`svn rev` < 49870) {
  setGeneric("rcond", function(x, norm, ...) standardGeneric("rcond"),
	   signature = c("x", "norm"),
	   useAsDefault = function(x, norm, ...) base::rcond(x, norm=norm, ...))
}## else:  rcond() has implicit generic in newer versions of R

    setGeneric("Schur", function(x, vectors, ...) standardGeneric("Schur"))

    setGeneric("unpack", function(x, ...) standardGeneric("unpack"))

##-     setGeneric("%p%", function(a, b) standardGeneric("%p%"))

    setGeneric("expm", function(x) standardGeneric("expm"))

##     setGeneric("writeHB", function(obj, file, ...)
##                standardGeneric("writeHB"))

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
##' <description>
##' Updates the mean vector mu given the linear predictor
##' gamma. Evaluate the residuals and the weighted sum of squared
##' residuals.
##' <details>
##' Note that the offset is added to the linear predictor before
##' calculating mu.
##' 
##' The sqrtXwt matrix can be updated but the sqrtrwt should not be in
##' that the weighted sum of squared residuals should be calculated
##' relative to fixed weights.  Reweighting is done in a separate call.
##' @title Update the fitted mean response
##' @param respM a response module
##' @param gamma the value of the linear predictor before adding the offset
##' @param ... 
##' @return updated respM
setGeneric("updateMu", function(respM, gamma, ...)
           standardGeneric("updateMu")) 

##' <description>
##' Update the weights, sqrtrwt and sqrtXwt
##' <details>
##' @title Update the residual and X weights
##' @param respM a response module
##' @param ... 
##' @return updated response module
setGeneric("updateWts", function(respM, ...)
           standardGeneric("updateWts"))

if (FALSE) {                            # don't need this generic in R
##' <description>
##' Set new values of the coefficients.  Can be called with a single
##' vector argument and with a pair of vectors, representing a base and
##' an increment, plus a step factor.
##' <details>
##' @title set new values of the coefficients
##' @param predM a predictor module
##' @param base coefficient base value
##' @param incr increment
##' @param step step factor, defaults to 0 in which case incr is ignored
##' @param ... 
##' @return predM
setGeneric("setCoef", function(predM, base, incr, step = 0, ...) standardGeneric("setCoef"))
}

##' <description>
##' Update any internal structures associated with sqrtXwt and the
##' weighted residuals.  The "V" matrix is evaluated from X using the
##' sqrtXwt matrix and a Vtr vector is calculated.
##' <details>
##' @title 
##' @param predM a predictor module
##' @param sqrtXwt the sqrtXwt matrix
##' @param wtres the vector of weighted residuals
##' @param ... 
##' @return updated predM
setGeneric("reweight", function(predM, sqrtXwt, wtres, ...)
           standardGeneric("reweight"))

if (FALSE) {                       # don't nee this generic in R
##' <description>
##' Return the gamma vector
##' <details>
##' @title 
##' @param predM a predictor module
##' @param ... 
##' @return X %*% coef
setGeneric("gammaInc", function(predM, ...)
           standardGeneric("gammaInc"))
}

##' <description>
##' Solve for the coefficients, usually in the form of
##' coef <- solve(predM@fac, predM@Vtr, system = "A")
##' <details>
##' The squared length of the intermediate solution is attached as an
##' attribute of the returned value. 
##' @title solve for the coefficients or coefficient increment
##' @param predM 
##' @param ... 
##' @return coefficient vector or increment
setGeneric("solveCoef", function(predM, ...)
           standardGeneric("solveCoef"))
