### pdLogChol - a general positive definite structure parameterized
###   by the non-zero elements of the Cholesky factor.  The logarithms of
###   the diagonal elements are the first Ncol elements of the parameter
###   vector

setGeneric('pdLogChol',
           function(value, form=formula(NULL), nam = character(), data=list(),
                    ...)
           standardGeneric('pdLogChol'))

setMethod('pdLogChol',
          signature(value = 'formula', form = 'missing',
                    nam = 'missing', data = 'missing'),
          function(value, form, nam, data, ...) {
              new('pdLogChol', form = value)
          })

## Methods for the pdLogChol class

setReplaceMethod("coef",
                 signature(object = "pdLogChol", value = "numeric"),
                 function(object, value) {
                     .Call("pdLogChol_coefGets", object, value,
                           PACKAGE = "Matrix")
                 })

setAs("pdLogChol", "pdmatrix",
      function(from) new("pdmatrix", crossprod(from@factor)),
      function(from, value) {
          nc <- ncol(value)
          if (!identical(nc, dim(value)[2]))
              stop("value must be a square matrix")
          if (length(from@param) < 1) {
              from@Ncol <- nc
          }
          if (from@Ncol != nc)
              stop("can not change length of an initialized pdMat object")
          Names <- dimnames(value)[[2]]
          if (!is.null(Names))
              from@Names <- Names
          fact <- .Call("nlme_Chol", as(value, "matrix"), PACKAGE="Matrix")
          from@factor <- fact
          from@logDet = sum(log(diag(fact)))
          from@param <- c(log(diag(fact)), fact[col(fact) > row(fact)])
          from
      })

setMethod("solve", signature(a="pdLogChol", b="missing"),
          function(a, b) {
              if (!isInitialized(a))
                  stop(paste("Uninitialized", class(a), "object"))
              as(a, "pdmatrix") <- crossprod(t(solve(a@factor)))
              a
          })

setMethod("summary", signature(object="pdLogChol"),
          function(object, structName, noCorrelation, ...) {
              if (missing(structName)) structName =
                   "General positive-definite, Log-Cholesky parametrization"
              if (missing(noCorrelation)) noCorrelation = FALSE
              callNextMethod()
          })

