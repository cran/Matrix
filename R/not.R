## METHODS FOR GENERIC: ! (not)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("!", c(x = "Matrix"),
          function(x) !.M2kind(x, "l"))

setMethod("!", c(x = "sparseVector"),
          function(x) !.V2kind(x, "l"))

setMethod("!", c(x = "ndenseMatrix"),
          function(x) {
              if(.M.shape(x) == "t")
                  x <- .M2gen(x)
              x@x <- { y <- x@x; if(anyNA(y)) !(y | is.na(y)) else !y }
              x
          })
setMethod("!", c(x = "ldenseMatrix"),
          function(x) {
              if(.M.shape(x) == "t")
                  x <- .M2gen(x)
              x@x <- !x@x
              x
          })

setMethod("!", c(x = "nsparseMatrix"),
          function(x) {
              x <- .sparse2dense(if(.M.shape(x) == "t") .M2gen(x) else x)
              x@x <- !x@x
              x
          })
setMethod("!", c(x = "lsparseMatrix"),
          function(x) {
              x <- .sparse2dense(if(.M.shape(x) == "t") .M2gen(x) else x)
              x@x <- !x@x
              x
          })

setMethod("!", c(x = "ndiMatrix"),
          function(x) {
              if(x@diag == "N" && anyNA(y <- x@x))
                  x@x <- y | is.na(y)
              x <- .diag2dense(x, ".", "g")
              x@x <- !x@x
              x
          })
setMethod("!", c(x = "ldiMatrix"),
          function(x) {
              x <- .diag2dense(x, ".", "g")
              x@x <- !x@x
              x
          })

setMethod("!", c(x = "nsparseVector"),
          function(x) !.V2v(x))
setMethod("!", c(x = "lsparseVector"),
          function(x) !.V2v(x))

setMethod("!", c(x = "indMatrix"),
          function(x) {
              x <- .ind2dense(x)
              x@x <- !x@x
              x
          })
