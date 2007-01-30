#### Permutation Matrices -- Coercion and Methods

## The typical   'constructor' : coerce from  'index'
setAs("integer", "pMatrix",
      function(from) {
          n <- length(from)
          nn <- names(from)
          new("pMatrix", Dim = rep.int(n, 2), Dimnames = list(nn,nn),
              perm = from)
      })

setAs("numeric", "pMatrix",
      function(from)
	  if(all(from == (i <- as.integer(from)))) as(i, "pMatrix")
	  else stop("coercion to 'pMatrix' only works from integer numeric"))

setAs("pMatrix", "matrix",
      function(from) {
	  fp <- from@perm
	  r <- diag(nrow = length(fp))[fp,]
	  if(.has.DN(from)) dimnames(r) <- from@Dimnames
	  r
      })

## coerce to 0/1 sparse matrix, i.e. sparse pattern
setAs("pMatrix", "ngTMatrix",
      function(from) {
          d <- from@Dim
	  new("ngTMatrix", i = seq_len(d[1]) - 1:1, j = from@perm - 1:1,
              Dim = d, Dimnames = from@Dimnames)
      })

setAs("pMatrix", "TsparseMatrix", function(from) as(from, "ngTMatrix"))
setAs("pMatrix", "nMatrix",	  function(from) as(from, "ngTMatrix"))
setAs("pMatrix", "lMatrix", function(from) as(as(from, "nMatrix"),"lMatrix"))

setAs("pMatrix", "CsparseMatrix",
      function(from) as(as(from, "ngTMatrix"), "CsparseMatrix"))


setMethod("solve", signature(a = "pMatrix", b = "missing"),
	  function(a, b) {
	      bp <- ap <- a@perm
	      bp[ap] <- seq_along(ap)
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


.pMat.nosense <- function (x, i, j, ..., value)
    stop('partially replacing "pMatrix" entries is not sensible')
setReplaceMethod("[", signature(x = "pMatrix", i = "index"), .pMat.nosense)
setReplaceMethod("[", signature(x = "pMatrix", i = "missing", j = "index"),
		 .pMat.nosense) ##   explicit  ^^^^^^^^^^^^ for disambiguation
