### pdCompSymm: Compound symmetry structure

setGeneric('pdCompSymm',
           function(value, form, nam, data, ...)
           standardGeneric('pdCompSymm'))

setMethod("pdCompSymm",
          signature(value = 'formula', form = 'missing',
                    nam = 'missing', data = 'missing'),
          function(value, form, nam, data, ...) {
              new('pdCompSymm', form = value)
          })

setMethod("isInitialized", "pdCompSymm",
          function(object) (length(object@ncol) != 0))

setReplaceMethod("coef",
                 signature(object = "pdCompSymm", value = "numeric"),
                 function(object, value) {
                     .Call("pdCompSymm_coefGets", object, value,
                           PACKAGE = "Matrix")
                 })

setAs('pdCompSymm', 'pdmatrix',
      function(from) {
          if (length(from@Ncol) == 1 && from@Ncol >= 1 &&
              length(from@param) == 2) {
              value <- exp(2 * from@param) * diag(from@Ncol)
              nam <- from@Names
              if (length(nam) == ncol(value)) {
                  dimnames(value) <- list(nam, nam)
              }
              return(value)
          } else {
              stop("Uninitialized pdCompSymm object")
          }
      })
