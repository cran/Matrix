setAs("trMatrix", "geMatrix",
      function(from) {
          .Call("trMatrix_as_geMatrix", from, PACKAGE="Matrix")
      })

setAs("trMatrix", "matrix",
      function(from) {
          .Call("trMatrix_as_matrix", from, PACKAGE="Matrix")
      })

setMethod("%*%", signature(x = "trMatrix", y = "geMatrix"),
          function(x, y)
          .Call("trMatrix_geMatrix_mm", x, y, PACKAGE = "Matrix"))

setMethod("%*%", signature(x = "geMatrix", y = "trMatrix"),
          function(x, y)
          .Call("trMatrix_geMatrix_mm_R", y, x, PACKAGE = "Matrix"))

setMethod("crossprod", signature(x = "trMatrix", y = "missing"),
          function(x, y = NULL) crossprod(as(x, "geMatrix")),
          valueClass = "poMatrix")

setMethod("determinant", signature(x = "trMatrix", logarithm = "missing"),
          function(x, logarithm, ...) determinant(x, TRUE))

setMethod("determinant", signature(x = "trMatrix", logarithm = "logical"),
          function(x, logarithm, ...) {
              dg = diag(x)
              if (logarithm) {
                  modulus = sum(log(abs(dg)))
                  sgn = prod(sign(dg))
              } else {
                  modulus = prod(dg)
                  sgn = sign(modulus)
                  modulus = abs(modulus)
              }
              attr(modulus, "logarithm") = logarithm
              val = list(modulus = modulus, sign = sgn)
              class(val) = "det"
              val
          })

setMethod("norm", signature(x = "trMatrix", type = "character"),
          function(x, type, ...)
          .Call("trMatrix_norm", x, type, PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("norm", signature(x = "trMatrix", type = "missing"),
          function(x, type, ...)
          .Call("trMatrix_norm", x, "O", PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("rcond", signature(x = "trMatrix", type = "character"),
          function(x, type, ...)
          .Call("trMatrix_rcond", x, type, PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("rcond", signature(x = "trMatrix", type = "missing"),
          function(x, type, ...)
          .Call("trMatrix_rcond", x, "O", PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("solve", signature(a = "trMatrix", b="missing"),
          function(a, b, ...)
          .Call("trMatrix_solve", a, PACKAGE = "Matrix"),
          valueClass = "trMatrix")

setMethod("solve", signature(a = "trMatrix", b="matrix"),
          function(a, b, ...)
          .Call("trMatrix_matrix_solve", a, b, PACKAGE = "Matrix"),
          valueClass = "matrix")

setMethod("t", signature(x = "trMatrix"),
          function(x) {
              val = new("trMatrix", Dim = x@Dim, x = t(as(x, "matrix"))@x)
              if (x@uplo == "U") val@uplo = "L"
              if (x@diag == "U") val@diag = "U"
              val
          }, valueClass = "trMatrix")
