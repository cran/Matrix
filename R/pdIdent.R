### pdIdent: multiple of the identity matrix - the parameter is
###  half the logarithm of the multiple.

## constructor for pdIdent

setGeneric('pdIdent',
           function(value, form, nam, data, ...)
           standardGeneric('pdIdent'))

setMethod("pdIdent",
          signature(value = 'formula', form = 'missing',
                    nam = 'missing', data = 'missing'),
          function(value, form, nam, data, ...) {
              new('pdIdent', form = value)
          })

## methods for the pdIdent class

setReplaceMethod("coef",
                 signature(object = "pdIdent", value = "numeric"),
                 function(object, value) {
                     if (length(value) != 1)
                         stop("coef for pdIdent class must be of length 1")
                     if (length(object@Ncol) < 1 || object@Ncol < 1)
                         stop(paste("Changing parameter of uninitialized",
                                    class(object), "object"))
                     object@param <- value
                     diag(object@factor) <- exp(value)
                     object@logDet <- object@Ncol * object@param
                     object
                 })

setAs('pdIdent', 'pdmatrix',
      function(from) {
          if (!isInitialized(from))
              stop("Uninitialized pdIdent object")
          Ncol <- from@Ncol
          value <- diag(exp(2 * from@param), Ncol)
          nam <- from@Names
          if (length(nam) == Ncol) {
              dimnames(value) <- list(nam, nam)
          }
          new("pdmatrix", value)
      },
      function(from, value) {
          nc <- ncol(value)
          if (!identical(nc, dim(value)[1]))
              stop("value must be a square matrix")
          if (length(from@param) < 1) {
              from@Ncol <- nc
              from@factor <- diag(nrow = nc, ncol = nc)
          }
          if (from@Ncol != nc)
              stop("can not change length of an initialized pdMat object")
          if (length(colnames(value)) != nc) {
              if (length(from@Names) == 0)
                  from@Names <- paste("V", 1:nc, sep = "")
          } else {
              from@Names <- colnames(value)
          }
          from@param <- mean(log(diag(value)))/2
          diag(from@factor) <-  exp(from@param)
          from@logDet <- nc * from@param
          from
      })

setAs("pdIdent", "corrmatrix",
      function(from) {
          if (!isInitialized(from))
              stop(paste("Uninitialized", class(from), "object"))
          Ncol <- from@Ncol
          val <- diag(Ncol)
          stdDev <- rep(exp(from@param), Ncol)
          if (length(nm <- from@Names) == 0) {
              nm <- paste("V", 1:Ncol, sep = "")
          }
          dimnames(val) <- list(nm, nm)
          names(stdDev) <- nm
          new("corrmatrix", val, stdDev = stdDev)
      })

setMethod("solve", signature(a="pdIdent", b="missing"),
          function(a, b)
      {
          if (!isInitialized(a))
              stop(paste("Uninitialized", class(a), "object"))
          a@param <- -a@param
          diag(a@factor) <- exp(a@param)
          a@logDet <- a@Ncol * a@param
          a
      })
