### pdNatural - a general positive definite structure parameterized
###   by the log of the square root of the diagonal elements and the
###   generalized logit of the correlations. This is NOT an unrestricted
###   parametrization

setGeneric('pdNatural',
           function(value, form=formula(NULL), nam = character(),
                    data=list(), ...)
           standardGeneric('pdNatural'))

setMethod('pdNatural',
          signature(value = 'formula', form = 'missing',
                    nam = 'missing', data = 'missing'),
          function(value, form, nam, data, ...) {
              new('pdNatural', form = value)
          })

setMethod('pdNatural',
          signature(value = 'pdMat', form = 'missing',
                    nam = 'missing', data = 'missing'),
          function(value, form, nam, data) {
              val <- new('pdNatural', form = value@form, Names = value@Names)
              as(val, "corrmatrix") <- as(value, "corrmatrix")
              val
          })


setReplaceMethod("coef",
                 signature(object = "pdNatural", value = "numeric"),
                 function(object, value) {
                     npar <- length(value)
                     lenPar <- length(object@param)
                     if (npar != lenPar) {
                         Ncol <- round((sqrt(8*length(value) + 1) - 1)/2)
                         np <- (Ncol * (Ncol + 1))/2
                         if (np != npar)
                             stop(paste("coef for a pdNatural object cannot have",
                                        "length", npar))
                         if (lenPar <= 0 && length(object@Names) > 0) {
                             lenPar <- length(object@Names)
                             lenPar <- (lenPar * (lenPar+1))/2
                         }
                         if (lenPar && lenPar != npar)
                             stop("coef for a pdNatural object has inconsistent length")
                         object@Ncol <- as.integer(Ncol)
                         object@factor <- matrix(0., Ncol, Ncol)
                     }
                     Ncol <- object@Ncol
                     corr <- diag(nrow = Ncol, ncol = Ncol)
                     expz <- exp(value[-(1:Ncol)])
                     corr[upper.tri(corr)] <- (expz - 1)/(1 + expz)
                     corr[lower.tri(corr)] <- t(corr)[lower.tri(corr)]
                     corrFact <- La.chol(corr)
                     object@factor <- t(exp(value[1:Ncol]) * t(corrFact))
                     object@param <- value
                     object@logDet <-
                         sum(log(diag(corrFact)))+sum(value[1:Ncol])
                     object
                 })

setAs("pdNatural", "corrmatrix",
      function(from) {
          if (!isInitialized(from))
              stop(paste("Uninitialized", class(from), "object"))
          value <- .Call("pdNatural_corrmatrix", from, PACKAGE = "Matrix")
          if (length(from@Names) == from@Ncol)
              dimnames(value) <- list(from@Names, from@Names)
          new("corrmatrix", value, stdDev = exp(from@param[1:from@Ncol]))
      },
      function(from, value) {
          nc <- ncol(value)
          if (!identical(nc, dim(value)[1]))
              stop("value must be a square matrix")
          if (length(from@param) < 1) {
              from@Ncol <- nc
          }
          if (from@Ncol != nc)
              stop("can not change length of an initialized pdMat object")
          Names <- dimnames(value)[[2]]
          if (!is.null(Names))
              from@Names <- Names
          rho <- value[upper.tri(value)]
          coef(from) <- c(log(value@stdDev), log((1.+rho)/(1.-rho)))
          from
      })

setAs("pdNatural", "pdmatrix",
      function(from) {
          if (!isInitialized(from))
              stop(paste("Uninitialized", class(from), "object"))
          new("pdmatrix", .Call("pdNatural_pdmatrix", from, PACKAGE="Matrix"))
      },
      function(from, value) {
          stdDev <- sqrt(diag(value))
          value <- array(t(value/stdDev)/stdDev, dim(value))
          colNames <- colnames(value)
          diag(value) <- rep(1.0, from@Ncol)
          if (!is.null(colNames)) {
              dimnames(value) <- list(colNames, colNames)
              names(stdDev) <- colNames
          }
          as(from, "corrmatrix") <-
              new("corrmatrix", value, stdDev = stdDev)
          from
      })

setMethod("solve", signature(a = "pdNatural", b = "missing"),
          function(a, b) {
              if (!isInitialized(a))
                  stop(paste("Unitialized", class(a), "object"))
              Ncol <- a@Ncol
              if (Ncol > 1) {
                  as(a, "pdmatrix") <- solve(as(a, "pdmatrix"))
              } else {
                  coef(a) <- -coef(a)
              }
              a
          })

#setMethod("LMEgradient",
#          signature(x="pdNatural", A="matrix", nlev="numeric"),
#          function(x, A, nlev)
#          .Call("pdNatural_LMEgradient", x, A, nlev, PACKAGE="lme4")
#          )

#setMethod("summary", signature(object="pdNatural"),
#          function(object,
#                   structName = "General positive-definite, Natural parametrization",
#                   ...) {
#              callNextMethod()
#          })
