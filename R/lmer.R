setReplaceMethod("coef", signature(object = "lmer", value = "numeric"),
                 function(object, unconst = FALSE, ..., value)
                 .Call("lmer_coefGets", object, as.double(value),
                       unconst))

setMethod("coef", signature(object = "lmer"),
          function(object, unconst = FALSE, ...) {
              .Call("lmer_coef", object, unconst)
          })

setMethod("deviance", "lmer",
          function(object, REML = FALSE, ...) {
              chol(object)
              object@deviance[[ifelse(REML, "REML", "ML")]]
          })

setMethod("chol", signature(x = "lmer"),
          function(x, pivot = FALSE, LINPACK = pivot)
          .Call("lmer_factor", x)
          )

setMethod("solve", signature(a = "lmer", b = "missing"),
          function(a, b, ...)
          .Call("lmer_invert", a)
          )


