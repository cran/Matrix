setReplaceMethod("coef", signature(object = "ssclme", value = "numeric"),
                 function(object, unconst = FALSE, ..., value)
                 .Call("ssclme_coefGets", object, as.double(value),
                       unconst, PACKAGE = "Matrix"))

setAs("ssclme", "sscMatrix",
      function(from)
      new("sscMatrix", i = from@i, p = from@p, Dim = from@Dim, x = from@x,
          uplo = "U"))

setAs("ssclme", "tscMatrix",
      function(from)
      new("tscMatrix", i = from@Li, p = from@Lp, Dim = from@Dim, x = from@Lx,
          uplo = "L", diag = "U"))

setMethod("coef", signature(object = "ssclme"),
          function(object, unconst = FALSE, ...) {
              .Call("ssclme_coef", object, unconst, PACKAGE = "Matrix")
          })

