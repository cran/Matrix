#### --- All method definitions for  "!" (not) ---

### R >= 2.6.0 needs (x)
### R <= 2.5.x had   (e1)

if(getRversion() >= "2.6.0") {

## Divert everything to  "lMatrix" and its subclasses :
setMethod("!", "Matrix", function(x) !as(x, "lMatrix"))

## -- diag ---

setMethod("!", "ldiMatrix", function(x) {
    if(x@diag == "N")
	x@x <- !x@x
    else { ## "U"
	x@diag <- "N"
	x@x <- rep.int(FALSE, x@Dim[1])
    }
    x
})

## -- lsparse --

setMethod("!", "lsparseMatrix",
          ## turns FALSE to TRUE --> dense matrix
          function(x) !as(x, "denseMatrix"))# was "lgeMatrix"

## Use "Matrix" method !as(. , "lMatrix")
## setMethod("!", "nsparseMatrix",
##           ## turns FALSE to TRUE --> dense matrix
##           function(x) !as(x, "ngeMatrix"))


## -- ldense ---

setMethod("!", "ltrMatrix",
	  function(x) {
	      x@x <- !x@x
	      ## And now we must fill one triangle with '!FALSE' results :

	      ## TODO: the following should be .Call using
	      ##	a variation of make_array_triangular:
	      r <- as(x, "lgeMatrix")
	      n <- x@Dim[1]
	      coli <- rep(1:n, each=n)
	      rowi <- rep(1:n, n)
	      Udiag <- x@diag == "U"
	      log.i <-
		  if(x@uplo == "U") {
		      if(Udiag) rowi >= coli else rowi > coli
		  } else {
		      if(Udiag) rowi <= coli else rowi < coli
		  }
	      r@x[log.i] <- TRUE
	      r
	  })

setMethod("!", "ltpMatrix", function(x) !as(x, "ltrMatrix"))

## for the other ldense* ones
setMethod("!", "lgeMatrix",
	  function(x) { x@x <- !x@x ; x })
## FIXME : this loses symmetry "lsy" and "lsp":
setMethod("!", "ldenseMatrix",
	  function(x) !as(x, "lgeMatrix"))

## -- ndense ---

setMethod("!", "ntrMatrix",
	  function(x) {
	      x@x <- !x@x
	      ## And now we must fill one triangle with '!FALSE' results :

	      ## TODO: the following should be .Call using
	      ##	a variation of make_array_triangular:
	      r <- as(x, "ngeMatrix")
	      n <- x@Dim[1]
	      coli <- rep(1:n, each=n)
	      rowi <- rep(1:n, n)
	      Udiag <- x@diag == "U"
	      log.i <-
		  if(x@uplo == "U") {
		      if(Udiag) rowi >= coli else rowi > coli
		  } else {
		      if(Udiag) rowi <= coli else rowi < coli
		  }
	      r@x[log.i] <- TRUE
	      r
	  })

setMethod("!", "ntpMatrix", function(x) !as(x, "ntrMatrix"))

## for the other ldense* ones
setMethod("!", "ngeMatrix",
          function(x) { x@x <- !x@x ; x })
## FIXME : this loses symmetry "nsy" and "nsp":
setMethod("!", "ndenseMatrix",
          function(x) !as(x, "ngeMatrix"))

##-----------------------------------------------
} else { ## R <= 2.5.x ---- the same thing, with 'e1' instead of 'x'
##-----------------------------------------------

## Divert everything to  "lMatrix" and its subclasses :
setMethod("!", "Matrix", function(e1) !as(e1, "lMatrix"))

## -- diag ---

setMethod("!", "ldiMatrix", function(e1) {
    if(e1@diag == "N")
	e1@x <- !e1@x
    else { ## "U"
	e1@diag <- "N"
	e1@x <- rep.int(FALSE, e1@Dim[1])
    }
    e1
})

## -- lsparse --

setMethod("!", "lsparseMatrix",
          ## turns FALSE to TRUE --> dense matrix
          function(e1) !as(e1, "denseMatrix"))# was "lgeMatrix"

## Use "Matrix" method !as(. , "lMatrix")
## setMethod("!", "nsparseMatrix",
##           ## turns FALSE to TRUE --> dense matrix
##           function(e1) !as(e1, "ngeMatrix"))


## -- ldense ---

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

## -- ndense ---

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

}
