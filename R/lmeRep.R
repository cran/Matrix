setReplaceMethod("coef", signature(object = "lmeRep", value = "numeric"),
                 function(object, unconst = FALSE, ..., value)
                 .Call("lmeRep_coefGets", object, as.double(value),
                       unconst, PACKAGE = "Matrix"))

setAs("ssclme", "sscMatrix",
      function(from)
      new("sscMatrix", i = from@i, p = from@p, Dim = from@Dim, x = from@x,
          uplo = "U"))

setAs("ssclme", "tscMatrix",
      function(from)
      new("tscMatrix", i = from@Li, p = from@Lp, Dim = from@Dim, x = from@Lx,
          uplo = "L", diag = "U"))

setMethod("coef", signature(object = "lmeRep"),
          function(object, unconst = FALSE, ...) {
              .Call("lmeRep_coef", object, unconst, PACKAGE = "Matrix")
          })

setMethod("deviance", "lmeRep",
          function(object, REML = FALSE, ...) {
              chol(object)
              object@deviance[[ifelse(REML, "REML", "ML")]]
          })

setMethod("chol", signature(x = "lmeRep"),
          function(x, pivot = FALSE, LINPACK = pivot)
          .Call("lmeRep_factor", x, PACKAGE = "Matrix")
          )

setMethod("solve", signature(a = "lmeRep", b = "missing"),
          function(a, b, ...)
          .Call("lmeRep_invert", a, PACKAGE = "Matrix")
          )


