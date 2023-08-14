## METHODS FOR CLASS: triangularMatrix (virtual)
## triangular matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("isTriangular", signature(object = "triangularMatrix"),
          function(object, upper = NA, ...) {
              if(is.na(upper))
                  `attr<-`(TRUE, "kind", object@uplo)
              else
                  object@uplo == (if(upper) "U" else "L") || isDiagonal(object)
          })

## NB: [dz]t.Matrix should _not_ use this method as it does not
## tolerate numerical fuzz
setMethod("isSymmetric", signature(object = "triangularMatrix"),
          function(object, checkDN = TRUE, ...) {
              if(checkDN) {
                  ca <- function(check.attributes = TRUE, ...) check.attributes
                  if(ca(...) && !isSymmetricDN(object@Dimnames))
                      return(FALSE)
              }
              isDiagonal(object)
          })
