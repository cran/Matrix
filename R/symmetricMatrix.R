## METHODS FOR CLASS: symmetricMatrix (virtual)
## Hermitian {incl. real, symmetric} matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dimnames", signature(x = "symmetricMatrix"),
          function(x) symmDN(x@Dimnames))

setMethod("isSymmetric", signature(object = "symmetricMatrix"),
          function(object, ...) TRUE)

setMethod("isTriangular", signature(object = "symmetricMatrix"),
          function(object, upper = NA, ...) {
              if(!isDiagonal(object))
                  FALSE
              else if(is.na(upper))
                  `attr<-`(TRUE, "kind", "U")
              else TRUE
          })
