## General methods for the virtual class pdMat

setMethod("formula", "pdMat", function(x, ...) x@form,
          valueClass = "formula")

setMethod("dim", "pdMat", function(x) c(x@Ncol, x@Ncol),
          valueClass = "integer")

setMethod("names", "pdMat", function(x) x@Names,
          valueClass = "character")

setMethod("coef", signature(object="pdMat"),
          function(object, ...) object@param,
          valueClass = "numeric")

setMethod("isInitialized", "pdMat",
          function(object) length(object@param) != 0,
          valueClass = "logical")

setReplaceMethod("names", "pdMat",
                 function(x, value) {
                     x@Names <- value
                     x
                 })

## The generics Names and Names<- will be deprecated in nlme_4.0
#setMethod("Names", "ANY", function(object, ...) names(object))
#setReplaceMethod("Names", "ANY",
#                 function(object, value) "names<-"(object, value))

## In general we use coersion to the pdmatrix class as the pdMatrix method

setMethod("pdMatrix", signature(object="pdMat"),
          function(object) as(object, "pdmatrix"))

setAs("pdMat", "matrix",
      function(from) as(as(from, "pdmatrix"), "matrix"),
      function(from, value) {
          as(from, "pdmatrix") <- new("pdmatrix", value)
          from
      })

## In general we use coersion to the pdfactor class as the pdFactor method

setMethod("pdFactor", signature(object="pdMat"),
          function(object) as(object, "pdfactor"),
          valueClass = "pdfactor")

## This coersion should be overridden in classes with a simpler way of
## obtaining the factor

setAs("pdMat", "pdfactor",
      function(from) {
          if (!isInitialized(from))
              stop(paste("Uninitialized", class(from), "object"))
          new("pdfactor", from@factor, logDet = from@logDet)
      },
      function(from, value) {
          as(from, "pdmatrix") <- crossprod(value)
          from
      })

setAs("matrix", "pdmatrix",
      function(from) new("pdmatrix", from))

setAs("pdfactor", "pdmatrix",
      function(from) new("pdmatrix", crossprod(from)))

setAs("pdmatrix", "pdfactor",
      function(from) {
          val <- new("pdfactor",
                     .Call("nlme_Chol", as(from, "pdmatrix"), PACKAGE="Matrix"))
          val@logDet <- sum(log(diag(val)))
          val
      },
      function(from, value) {
          as(from, "pdmatrix") <- crossprod(value)
          from
      })

setReplaceMethod("matrix", signature(object="pdMat", value="matrix"),
                 function(object, value) {
                     as(object, "pdmatrix") <- value
                 })

#setMethod("logDet", signature(object="pdDiag", covariate="missing"),
#          function(object, covariate, ...) {
#              if (!isInitialized(object))
#                  stop(paste("Uninitialized", class(object), "object"))
#              object@logDet
#          },
#          valueClass = "numeric")

## In general we use coersion to the corrmatrix class as the corMatrix method

setMethod("corMatrix", signature(object="pdMat"),
          function(object) as(object, "corrmatrix"),
          valueClass = "corrmatrix")

## This should be overridden in classes with a simpler way of obtaining
## the correlation matrix

setAs("pdMat", "corrmatrix",
      function(from) {
          Var  <- as(from, "pdmatrix")
          stdDev <- sqrt(diag(Var))
          colNames <- colnames(Var)
          if (is.null(colNames))
              colNames <- rownames(Var, do.NULL = FALSE, prefix="V")
          names(stdDev) <- colNames
          value <- array(t(Var/stdDev)/stdDev, dim(Var),
                         list(colNames, colNames))
          diag(value) <- rep(1.0, from@Ncol)
          new("corrmatrix", value, stdDev = stdDev)
      })

## gradient of the positive-definite matrix with respect to the parameters
## This method uses finite differences.  It should be overridden in
## explicit pdMat classes for which the gradient can be calculated explicitly

#setMethod("pdgradient", "pdMat",
#          function(x) {
#              m0 <- as(x, "pdmatrix")
#              pars <- coef(x)
#              val <- array(0., c(dim(m0), length(pars)))
#              eps <- sqrt(.Machine$double.eps)
#              ind <- seq(along = pars)
#              for (i in ind) {
#                  coef(x) <- pars + eps * pars * (ind == i)
#                  val[,,i] <- (as(x, "pdmatrix") - m0)/(pars[i] * eps)
#              }
#              val
#          })

setMethod("show", "pdMat",
          function(object) {
              if (isInitialized(object)) {
                  cat("positive definite matrix of class",
                      class(object), "\n")
                  print(as(object, "pdmatrix"))
              } else {
                  cat(paste("Uninitialized", class(object), "object\n"))
              }
          })

setMethod("solve", signature(a="pdMat", b="missing"),
          function(a, b, ...) {
              if (!isInitialized(a))
                  stop(paste("Uninitialized", class(a),
                             "object", deparse(substitute(a))))
              as(a, "matrix") <- solve(as(a, "matrix"))
              a
          })

setMethod("summary", signature(object="pdMat"),
          function(object, structName, noCorrelation , ...) {
              cl = as.character(class(object))
              if (missing(structName))
                  structName = cl
              if (missing(noCorrelation))
                  noCorrelation = switch(cl, pdDiag = , pdIdent = TRUE, FALSE)
              if (isInitialized(object)) {
                  new("summary.pdMat", cor = as(object, "corrmatrix"),
                      structName = structName,
                      noCorrelation = noCorrelation,
                      formula = formula(object))
              } else {
                  object
              }
          })

setAs("pdmatrix", "corrmatrix",
      function(from) {
          ss = sqrt(diag(from))
          new("corrmatrix", t(from/ss)/ss, stdDev = ss)
      })

## This function is for testing only.  It will eventually be removed.
##
pdgradNumeric <- function(x)
{
    m0 <- as(x, "pdmatrix")
    pars <- coef(x)
    val <- array(0., c(dim(m0), length(pars)))
    eps <- sqrt(.Machine$double.eps)
    ind <- seq(along = pars)
    for (i in ind) {
        coef(x) <- pars + eps * pars * (ind == i)
        val[,,i] <- (as(x, "pdmatrix") - m0)/(pars[i] * eps)
    }
    val
}

