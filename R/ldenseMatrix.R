#### "ldenseMatrix" - virtual class of logical dense matrices
####  ------------
#### Contains  lge*;  ltr*, ltp*;  lsy*, lsp*;	 ldi*

## Logical -> Double {of same structure}:

setAs("lgeMatrix", "dgeMatrix", l2d_Matrix)
setAs("lsyMatrix", "dsyMatrix", l2d_Matrix)
setAs("lspMatrix", "dspMatrix", l2d_Matrix)
setAs("ltrMatrix", "dtrMatrix", l2d_Matrix)
setAs("ltpMatrix", "dtpMatrix", l2d_Matrix)

### NOTA BENE: Much of this is *very* parallel to ./ndenseMatrix.R
###						  ~~~~~~~~~~~~~~~~

## all need be coercable to "lgeMatrix":

setAs("lsyMatrix", "lgeMatrix",	 function(from)
      .Call(lsyMatrix_as_lgeMatrix, from, 0:0))
setAs("ltrMatrix", "lgeMatrix",	 function(from)
      .Call(ltrMatrix_as_lgeMatrix, from, 0:0))
setAs("ltpMatrix", "lgeMatrix",
      function(from) as(as(from, "ltrMatrix"), "lgeMatrix"))
setAs("lspMatrix", "lgeMatrix",
      function(from) as(as(from, "lsyMatrix"), "lgeMatrix"))
## and the reverse
setAs("lgeMatrix", "ltpMatrix",
      function(from) as(as(from, "ltrMatrix"), "ltpMatrix"))
setAs("lgeMatrix", "lspMatrix",
      function(from) as(as(from, "lsyMatrix"), "lspMatrix"))


## packed <->  non-packed :

setAs("lspMatrix", "lsyMatrix",
      function(from)
      .Call(lspMatrix_as_lsyMatrix, from, 0:0))

setAs("lsyMatrix", "lspMatrix",
      function(from)
      .Call(lsyMatrix_as_lspMatrix, from, 0:0))

setAs("ltpMatrix", "ltrMatrix",
      function(from)
      .Call(ltpMatrix_as_ltrMatrix, from, 0:0))

setAs("ltrMatrix", "ltpMatrix",
      function(from)
      .Call(ltrMatrix_as_ltpMatrix, from, 0:0))



### -> symmetric :

if(FALSE) ## not sure if this is a good idea ... -- FIXME?
setIs("lgeMatrix", "lsyMatrix",
      test = function(obj) isSymmetric(obj),
      replace = function(obj, value) { ## copy all slots
	  for(n in slotNames(obj)) slot(obj, n) <- slot(value, n)
      })

### Alternative (at least works):
setAs("lgeMatrix", "lsyMatrix",
      function(from) {
	  if(isSymmetric(from))
	      new("lsyMatrix", x = from@x, Dim = from@Dim,
		  Dimnames = from@Dimnames, factors = from@factors)
	  else stop("not a symmetric matrix")
      })

setAs("lgeMatrix", "ltrMatrix",
      function(from) {
	  if(isT <- isTriangular(from))
	      new("ltrMatrix", x = from@x, Dim = from@Dim,
		  Dimnames = from@Dimnames, uplo = attr(isT, "kind"))
	  ## TODO: also check 'diag'
	  else stop("not a triangular matrix")
      })


###  ldense* <-> "matrix" :

## 1) "lge* :
setAs("lgeMatrix", "matrix",
      function(from) array(from@x, dim = from@Dim, dimnames = from@Dimnames))

setAs("matrix", "lgeMatrix",
      function(from) {
	  new("lgeMatrix",
	      x = as.logical(from),
	      Dim = as.integer(dim(from)),
	      Dimnames = .M.DN(from))
      })

## 2) base others on "lge*":

setAs("matrix", "lsyMatrix",
      function(from) as(as(from, "lgeMatrix"), "lsyMatrix"))
setAs("matrix", "lspMatrix",
      function(from) as(as(from, "lsyMatrix"), "lspMatrix"))
setAs("matrix", "ltrMatrix",
      function(from) as(as(from, "lgeMatrix"), "ltrMatrix"))
setAs("matrix", "ltpMatrix",
      function(from) as(as(from, "ltrMatrix"), "ltpMatrix"))

## Useful if this was called e.g. for as(*, "lsyMatrix"), but it isn't
setAs("matrix", "ldenseMatrix", function(from) as(from, "lgeMatrix"))

setAs("ldenseMatrix", "matrix", ## uses the above l*M. -> lgeM.
      function(from) as(as(from, "lgeMatrix"), "matrix"))

## dense |-> compressed :

setAs("lgeMatrix", "lgTMatrix",
      function(from) {
	  ## Non'zeros':
	  nF <- nonFALSE(from@x)## == nz.NA(from@x, na. = TRUE)
	  ## cheap but not so efficient:
	  d <- dim(from)
	  ij <- which(array(nF, dim = d), arr.ind = TRUE) - 1:1
	  new("lgTMatrix", i = ij[,1], j = ij[,2], x = from@x[nF],
	      Dim = d, Dimnames = from@Dimnames,
	      factors = from@factors)
      })

setAs("lgeMatrix", "lgCMatrix",
      function(from) as(as(from, "lgTMatrix"), "lgCMatrix"))

setMethod("as.logical", signature(x = "ldenseMatrix"),
	  function(x, ...) as(x, "lgeMatrix")@x)

###----------------------------------------------------------------------


setMethod("t", signature(x = "lgeMatrix"), t_geMatrix)
setMethod("t", signature(x = "ltrMatrix"), t_trMatrix)
setMethod("t", signature(x = "lsyMatrix"), t_trMatrix)
setMethod("t", signature(x = "ltpMatrix"),
	  function(x) as(callGeneric(as(x, "ltrMatrix")), "ltpMatrix"))
setMethod("t", signature(x = "lspMatrix"),
	  function(x) as(callGeneric(as(x, "lsyMatrix")), "lspMatrix"))

setMethod("!", "ltrMatrix",
	  function(e1) {
	      e1@x <- !e1@x
	      ## And now we must fill one triangle with '!FALSE' results :

	      ## TODO: the following should be .Call using
	      ##	a variation of make_array_triangular:
	      r <- as(e1, "lgeMatrix")
	      n <- e1@Dim[1]
	      coli <- rep(1:n, each=n)
	      rowi <- rep(1:n, n)
	      Udiag <- e1@diag == "U"
	      log.i <-
		  if(e1@uplo == "U") {
		      if(Udiag) rowi >= coli else rowi > coli
		  } else {
		      if(Udiag) rowi <= coli else rowi < coli
		  }
	      r@x[log.i] <- TRUE
	      r
	  })

setMethod("!", "ltpMatrix", function(e1) !as(e1, "ltrMatrix"))

## for the other ldense* ones
setMethod("!", "lgeMatrix",
	  function(e1) { e1@x <- !e1@x ; e1 })
## FIXME : this loses symmetry "lsy" and "lsp":
setMethod("!", "ldenseMatrix",
	  function(e1) !as(e1, "lgeMatrix"))


setMethod("|", signature(e1="lgeMatrix", e2="lgeMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      e1@x <- e1@x | e2@x
	      e1
	  })
setMethod("&", signature(e1="lgeMatrix", e2="lgeMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      e1@x <- e1@x & e2@x
	      e1
	  })

setMethod("|", signature(e1="ldenseMatrix", e2="ldenseMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      as(e1, "lgeMatrix") | as(e2, "lgeMatrix")
	  })

setMethod("&", signature(e1="ldenseMatrix", e2="ldenseMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      as(e1, "lgeMatrix") & as(e2, "lgeMatrix")
	  })


setMethod("as.vector", signature(x = "ldenseMatrix", mode = "missing"),
	  function(x) as(x, "lgeMatrix")@x)
