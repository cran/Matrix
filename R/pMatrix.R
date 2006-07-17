#### Permutation Matrices -- Coercion and Methods

setAs("integer", "pMatrix",
      function(from) {
          n <- length(from)
          nn <- names(from)
          new("pMatrix", Dim = rep.int(n, 2), Dimnames = list(nn,nn),
              perm = from)
      })

setAs("pMatrix", "matrix",
      function(from) {
	  fp <- from@perm
	  r <- diag(nrow = length(fp))[fp,]
	  if(.has.DN(from)) dimnames(r) <- from@Dimnames
	  r
      })

## coerce to 0/1 sparse matrix, i.e. sparse logical :
setAs("pMatrix", "lgTMatrix",
      function(from) {
          d <- from@Dim
	  new("lgTMatrix", i = seq(length = d[1]) - 1:1, j = from@perm - 1:1,
              Dim = d, Dimnames = from@Dimnames)
      })

setAs("pMatrix", "TsparseMatrix",
      function(from) as(from, "lgTMatrix"))

setMethod("solve", signature(a = "pMatrix", b = "missing"),
          function(a, b) {
              bp <- ap <- a@perm
              bp[ap] <- seq(along = ap)
              new("pMatrix", perm = bp, Dim = a@Dim,
                  Dimnames = rev(a@Dimnames))
          }, valueClass = "pMatrix")

setMethod("t", signature(x = "pMatrix"), function(x) solve(x))

setMethod("%*%", signature(x = "matrix", y = "pMatrix"),
	  function(x, y) x[ , y@perm], valueClass = "matrix")

setMethod("%*%", signature(x = "pMatrix", y = "matrix"),
	  function(x, y) y[x@perm ,], valueClass = "matrix")

setMethod("%*%", signature(x = "pMatrix", y = "pMatrix"),
	  function(x, y) {
              stopifnot(identical(d <- x@Dim, y@Dim))
              n <- d[1]
              ## FIXME: dimnames dealing: as with S3 matrix's  %*%
              x@perm <- x@perm[y@perm]
              x
          })

setMethod("%*%", signature(x = "Matrix", y = "pMatrix"),
	  function(x, y) x[, y@perm])

setMethod("%*%", signature(x = "pMatrix", y = "Matrix"),
          function(x, y) y[x@perm , ])
