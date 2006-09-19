#### "ndenseMatrix" - virtual class of nonzero pattern dense matrices
####  ------------
#### Contains  nge*;  ntr*, ntp*;  nsy*, nsp*;   ndi*

## Nonzero Pattern -> Double {of same structure}:

setAs("ngeMatrix", "dgeMatrix", n2d_Matrix)
setAs("nsyMatrix", "dsyMatrix", n2d_Matrix)
setAs("nspMatrix", "dspMatrix", n2d_Matrix)
setAs("ntrMatrix", "dtrMatrix", n2d_Matrix)
setAs("ntpMatrix", "dtpMatrix", n2d_Matrix)

### NOTA BENE: Much of this is *very* parallel to ./ldenseMatrix.R
###						  ~~~~~~~~~~~~~~~~

## all need be coercable to "ngeMatrix":

setAs("nsyMatrix", "ngeMatrix",  function(from)
      .Call(lsyMatrix_as_lgeMatrix, from, 1:1))
setAs("ntrMatrix", "ngeMatrix",  function(from)
      .Call(ltrMatrix_as_lgeMatrix, from, 1:1))
setAs("ntpMatrix", "ngeMatrix",
      function(from) as(as(from, "ntrMatrix"), "ngeMatrix"))
setAs("nspMatrix", "ngeMatrix",
      function(from) as(as(from, "nsyMatrix"), "ngeMatrix"))
## and the reverse
setAs("ngeMatrix", "ntpMatrix",
      function(from) as(as(from, "ntrMatrix"), "ntpMatrix"))
setAs("ngeMatrix", "nspMatrix",
      function(from) as(as(from, "nsyMatrix"), "nspMatrix"))


## packed <->  non-packed :

setAs("nspMatrix", "nsyMatrix",
      function(from)
      .Call(lspMatrix_as_lsyMatrix, from, 1:1))

setAs("nsyMatrix", "nspMatrix",
      function(from)
      .Call(lsyMatrix_as_lspMatrix, from, 1:1))

setAs("ntpMatrix", "ntrMatrix",
      function(from)
      .Call(ltpMatrix_as_ltrMatrix, from, 1:1))

setAs("ntrMatrix", "ntpMatrix",
      function(from)
      .Call(ltrMatrix_as_ltpMatrix, from, 1:1))



### -> symmetric :

if(FALSE) ## not sure if this is a good idea ... -- FIXME?
setIs("ngeMatrix", "nsyMatrix",
      test = function(obj) isSymmetric(obj),
      replace = function(obj, value) { ## copy all slots
          for(n in slotNames(obj)) slot(obj, n) <- slot(value, n)
      })

### Alternative (at least works):
setAs("ngeMatrix", "nsyMatrix",
      function(from) {
	  if(isSymmetric(from))
	      new("nsyMatrix", x = from@x, Dim = from@Dim,
		  Dimnames = from@Dimnames, factors = from@factors)
	  else stop("not a symmetric matrix")
      })

setAs("ngeMatrix", "ntrMatrix",
      function(from) {
	  if(isT <- isTriangular(from))
	      new("ntrMatrix", x = from@x, Dim = from@Dim,
		  Dimnames = from@Dimnames, uplo = attr(isT, "kind"))
          ## TODO: also check 'diag'
	  else stop("not a triangular matrix")
      })


###  ldense* <-> "matrix" :

## 1) "nge* :
setAs("ngeMatrix", "matrix",
      function(from) array(from@x, dim = from@Dim, dimnames = from@Dimnames))

setAs("matrix", "ngeMatrix",
      function(from) {
	  new("ngeMatrix",
	      x = as.logical(from),
	      Dim = as.integer(dim(from)),
	      Dimnames = .M.DN(from))
      })

## 2) base others on "nge*":

setAs("matrix", "nsyMatrix",
      function(from) as(as(from, "ngeMatrix"), "nsyMatrix"))
setAs("matrix", "nspMatrix",
      function(from) as(as(from, "nsyMatrix"), "nspMatrix"))
setAs("matrix", "ntrMatrix",
      function(from) as(as(from, "ngeMatrix"), "ntrMatrix"))
setAs("matrix", "ntpMatrix",
      function(from) as(as(from, "ntrMatrix"), "ntpMatrix"))

## Useful if this was called e.g. for as(*, "nsyMatrix"), but it isn't
setAs("matrix", "ndenseMatrix", function(from) as(from, "ngeMatrix"))

setAs("ndenseMatrix", "matrix", ## uses the above l*M. -> lgeM.
      function(from) as(as(from, "ngeMatrix"), "matrix"))

## dense |-> compressed :

setAs("ngeMatrix", "ngTMatrix",
      function(from) {
          ##  cheap but not so efficient:
          ij <- which(as(from,"matrix"), arr.ind = TRUE) - 1:1
          new("ngTMatrix", i = ij[,1], j = ij[,2],
              Dim = from@Dim, Dimnames = from@Dimnames,
              factors = from@factors)
      })

setAs("ngeMatrix", "ngCMatrix",
      function(from) as(as(from, "ngTMatrix"), "ngCMatrix"))

setMethod("as.logical", signature(x = "ndenseMatrix"),
	  function(x, ...) as(x, "ngeMatrix")@x)

###----------------------------------------------------------------------


setMethod("t", signature(x = "ngeMatrix"), t_geMatrix)
setMethod("t", signature(x = "ntrMatrix"), t_trMatrix)
setMethod("t", signature(x = "nsyMatrix"), t_trMatrix)
setMethod("t", signature(x = "ntpMatrix"),
          function(x) as(callGeneric(as(x, "ntrMatrix")), "ntpMatrix"))
setMethod("t", signature(x = "nspMatrix"),
          function(x) as(callGeneric(as(x, "nsyMatrix")), "nspMatrix"))

setMethod("!", "ntrMatrix",
	  function(e1) {
	      e1@x <- !e1@x
	      ## And now we must fill one triangle with '!FALSE' results :

	      ## TODO: the following should be .Call using
	      ##	a variation of make_array_triangular:
	      r <- as(e1, "ngeMatrix")
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

setMethod("!", "ntpMatrix", function(e1) !as(e1, "ntrMatrix"))

## for the other ldense* ones
setMethod("!", "ngeMatrix",
          function(e1) { e1@x <- !e1@x ; e1 })
## FIXME : this loses symmetry "nsy" and "nsp":
setMethod("!", "ndenseMatrix",
          function(e1) !as(e1, "ngeMatrix"))


setMethod("|", signature(e1="ngeMatrix", e2="ngeMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      e1@x <- e1@x | e2@x
	      e1
	  })
setMethod("&", signature(e1="ngeMatrix", e2="ngeMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      e1@x <- e1@x & e2@x
	      e1
	  })

setMethod("|", signature(e1="ndenseMatrix", e2="ndenseMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      as(e1, "ngeMatrix") | as(e2, "ngeMatrix")
	  })

setMethod("&", signature(e1="ndenseMatrix", e2="ndenseMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      as(e1, "ngeMatrix") & as(e2, "ngeMatrix")
	  })

setMethod("as.vector", signature(x = "ndenseMatrix", mode = "missing"),
          function(x) as(x, "ngeMatrix")@x)
