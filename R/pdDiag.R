### pdDiag - diagonal structure parameterized by the logarithm of
###   the square root of the diagonal terms.

setGeneric('pdDiag',
           function(value, form, nam, data, ...)
           standardGeneric('pdDiag'))

setMethod("pdDiag",
          signature(value = 'formula', form = 'missing',
                    nam = 'missing', data = 'missing'),
          function(value, form, nam, data, ...) {
              new('pdDiag', form = value)
          })

## Methods for the pdDiag class

setReplaceMethod("coef",
                 signature(object = "pdDiag", value = "numeric"),
                 function(object, value) {
                     .Call("pdDiag_coefGets", object, value,
                           PACKAGE = "Matrix")
                 })

setAs("pdDiag", "pdmatrix",
      function(from) {
          if (!isInitialized(from))
              stop(paste("Uninitialized", class(from), "object"))
          value <- diag(exp(2 * from@param), ncol = from@Ncol)
          nam <- from@Names
          if (length(nam) == length(from@param)) {
              dimnames(value) <- list(nam, nam)
          }
          new("pdmatrix", value)
      },
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
          coef(from) <- log(diag(value))/2
          from
      })

setMethod("solve", signature(a="pdDiag", b="missing"),
          function(a, b) {
              if (!isInitialized(a))
                  stop(paste("Uninitialized", class(a), "object"))
              coef(a) <- -a@param
              a
          })
