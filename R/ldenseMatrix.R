#### "ldenseMatrix" - virtual class of logical dense matrices
####  ------------
#### Contains  lge*;  ltr*, ltp*;  lsy*, lsp*;   ldi*

## Logical -> Double {of same structure}:

setAs("lgeMatrix", "dgeMatrix", l2d_Matrix)
setAs("lsyMatrix", "dsyMatrix", l2d_Matrix)
setAs("lspMatrix", "dspMatrix", l2d_Matrix)
setAs("ltrMatrix", "dtrMatrix", l2d_Matrix)
setAs("ltpMatrix", "dtpMatrix", l2d_Matrix)

## all need be coercable to "lgeMatrix":

setAs("lsyMatrix", "lgeMatrix",  function(from)
      .Call(lsyMatrix_as_lgeMatrix, from))
setAs("ltrMatrix", "lgeMatrix",  function(from)
      .Call(ltrMatrix_as_lgeMatrix, from))
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
      .Call(lspMatrix_as_lsyMatrix, from))

setAs("lsyMatrix", "lspMatrix",
      function(from)
      .Call(lsyMatrix_as_lspMatrix, from))

setAs("ltpMatrix", "ltrMatrix",
      function(from)
      .Call(ltpMatrix_as_ltrMatrix, from))

setAs("ltrMatrix", "ltpMatrix",
      function(from)
      .Call(ltrMatrix_as_ltpMatrix, from))



### -> symmetric :

if(FALSE) ## cannot easily work around R bug  -- FIXME --
setIs("lgeMatrix", "lsyMatrix",
### BUG in R: this fails, because isSymmetric() is namespace hidden and NOT found
##B      test = function(obj) isSymmetric(obj),
##B and this fails too:
##B      test = function(obj) Matrix:::isSymmetric(obj),
      replace = function(obj, value) {
          ## copy all slots
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
	      Dimnames =
	      if(!is.null(dn <- dimnames(from))) dn else list(NULL,NULL)
	      )
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
          ##  cheap but not so efficient:
          ij <- which(as(from,"matrix"), arr.ind = TRUE) - 1:1
          new("lgTMatrix", i = ij[,1], j = ij[,2],
              Dim = from@Dim, Dimnames = from@Dimnames,
              factors = from@factors)
      })

setAs("lgeMatrix", "lgCMatrix",
      function(from) as(as(from, "lgTMatrix"), "lgCMatrix"))

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

## for the other ldense* ones:
setMethod("!", "ldenseMatrix",
          function(e1) { e1@x <- !e1@x ; e1 })

setMethod("as.vector", signature(x = "ldenseMatrix", mode = "missing"),
          function(x) as(x, "lgeMatrix")@x)
