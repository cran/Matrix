#### ------- All "Ops"  group methods for all Matrix classes ------------
###               ===
### Note that the "Ops" group consists of
### sub-groups   "Arith", "Compare", and "Logic"
###               -----    -------        -----
### where 'Arith'   :=  '"+"', '"-"', '"*"', '"^"', '"%%"', '"%/%"', '"/"'
###       'Compare' := '"=="', '">"', '"<"', '"!="', '"<="', '">="'
###       'Logic'   :=  '"&"', '"|"'  (( but *not* '"!"' since that has
###			                 only one argument ))

## in shell, find them with
##    egrep 'Method\("(Ops|Compare|Arith|Logic)"' R/*R

### "Ops" ---- remember Ops = {Arith, Compare, Logic}  [Logic: since R 2.4.1]
### -----


### Note: diagonalMatrix are handled by special methods -> ./diagMatrix.R
###       --------------                                   ~~~~~~~~~~~~~~


### --  0 -- (not dense *or* sparse) -----------------------------------

##-------- originally from ./Matrix.R --------------------

## Some ``Univariate'' "Arith":
setMethod("+", signature(e1 = "Matrix", e2 = "missing"), function(e1) e1)
## "fallback":
setMethod("-", signature(e1 = "Matrix", e2 = "missing"),
	  function(e1) {
	      warning("inefficient method used for \"- e1\"")
	      0-e1
	  })

## old-style matrices are made into new ones
setMethod("Ops", signature(e1 = "Matrix", e2 = "matrix"),
	  function(e1, e2) callGeneric(e1, Matrix(e2)))
##	    callGeneric(e1, Matrix(e2, sparse=is(e1,"sparseMatrix"))))
setMethod("Ops", signature(e1 = "matrix", e2 = "Matrix"),
	  function(e1, e2) callGeneric(Matrix(e1), e2))

## bail-outs -- on highest possible level, hence "Ops", not "Compare"/"Arith" :
setMethod("Ops", signature(e1 = "Matrix", e2 = "Matrix"),
          function(e1, e2) {
              d <- dimCheck(e1,e2)
              .bail.out.2(.Generic, class(e1), class(e2))
          })
setMethod("Ops", signature(e1 = "Matrix", e2 = "ANY"),
          function(e1, e2) .bail.out.2(.Generic, class(e1), class(e2)))
setMethod("Ops", signature(e1 = "ANY", e2 = "Matrix"),
          function(e1, e2) .bail.out.2(.Generic, class(e1), class(e2)))

## "General principle"
##  - - - - - - - - -
## For "Arith" it is sufficient (though not optimal, once we have "iMatrix"!)
## to define "dMatrix" methods and coerce all other "[nli]Matrix" to "dMatrix"
setMethod("Arith", signature(e1 = "Matrix", e2 = "Matrix"),
          function(e1, e2) callGeneric(as(e1, "dMatrix"), as(e2, "dMatrix")))

## For "Compare", this would be feasible too, but is clearly suboptimal,
## particularly for "==" and "!="
## and for "lMatrix" and "nMatrix"  should not coerce at all
## setMethod("Compare", signature(e1 = "Matrix", e2 = "Matrix"),
##           function(e1, e2) {
##               if(is.na(match(.Generic, c("==", "!="))))
##                   callGeneric(as(e1, "dMatrix"), as(e2, "dMatrix"))
##               else { ## no coercion needed for "==" or "!="
##               }
##           })


## Working entirely on "matching" x slot:
## can be done for matching-dim "*geMatrix", and also
## matching-{dim + uplo} for *packed* (only!) symmetric+triangular
.Ops.via.x <- function(e1,e2) {
    d <- dimCheck(e1, e2)
    e1@x <- callGeneric(e1@x, e2@x)
    e1
}


##-------- originally from ./dMatrix.R --------------------

## Note that there extra methods for <sparse> o <sparse> !
if(FALSE) ## use the (more sophisticated) "Ops" method below
setMethod("Arith", signature(e1 = "dMatrix", e2 = "dMatrix"),
          ##  "+", "-", "*", "^", "%%", "%/%", "/"
          ## Going -> dense* (= ddense*) -> dgeMatrix
	  function(e1, e2) {
	      d <- dimCheck(e1,e2)
	      callGeneric(as(e1, "denseMatrix"),
			  as(e2, "denseMatrix"))
	  })

## "Compare" -> returning  logical Matrices
setMethod("Compare", signature(e1 = "numeric", e2 = "dMatrix"),
	  function(e1,e2) {
	      ## "swap RHS and LHS" and use the method below:
	      switch(.Generic,
		     "==" =, "!=" = callGeneric(e2, e1),
		     "<"  = e2 >  e1,
		     "<=" = e2 >= e1,
		     ">"  = e2 <  e1,
		     ">=" = e2 <= e1)
	  })

setMethod("Compare", signature(e1 = "dMatrix", e2 = "numeric"),
	  function(e1, e2) { ## result will inherit from "lMatrix"
	      lClass <- class2(cl <- class(e1), "l")
              cl1 <- getClassDef(cl)
	      fullCl <-
                  if(extends(cl1, "symmetricMatrix")) "lsyMatrix" else "lgeMatrix"
	      r	 <- callGeneric(e1@x, e2)
	      r0 <- callGeneric(0, e2)
              d <- e1@Dim
	      ## trivial case first (beware of NA)
	      if(isTRUE(r0 && all(r))) {
		  r <- new(fullCl)
		  r@Dim <- d
		  r@Dimnames <- e1@Dimnames
		  r@x <- rep.int(TRUE, prod(d))
	      }
	      else if(extends(cl1, "denseMatrix")) {
		  full <- !isPacked(e1) # << both "dtr" and "dsy" are 'full'
		  if(full || identical(r0, FALSE) || extends(cl1,"symmetricMatrix"))
		      r <- new(lClass, x = r, Dim = d, Dimnames = dimnames(e1))
		  else { ## packed matrix with structural 0 and r0 is not FALSE:
		      ##--> result cannot be packed anymore
                      ## [dense & packed & not symmetric ] ==> must be "dtp*" :
                      if(!extends(cl1, "dtpMatrix"))
                          stop("internal bug in \"Compare\" method for \"dMatrix\"; please report")
                      rx <- rep.int(r0, d[1]*d[2])
                      rx[indTri(d[1], upper = (e1@uplo == "U"))] <- r
                      r <- new(fullCl, x = rx, Dim = d, Dimnames = dimnames(e1))
		  }

	      }
	      else { ## dsparseMatrix => lClass is "lsparse*"
		  if(extends(cl1, "triangularMatrix") && e1@diag == "U") {
		      e1 <- .diagU2N(e1, cl = cl1)
		      if(extends(cl1, "CsparseMatrix")) {
			  r <- callGeneric(e1@x, e2)
		      }
		      else {
			  ## CARE: assume that diagU2N appends the diagonal entries at end !
			  r <- c(r, rep.int(callGeneric(1, e2), d[1]))
		      }
		  }
		  if(identical(r0, FALSE)) { ## things remain sparse
		      if(!any(is.na(r)) && ((Ar <- all(r)) || !any(r))) {
			  r <- new(lClass)
			  r@Dim <- d
			  r@Dimnames <- dimnames(e1)
			  if(Ar) { # 'TRUE' instead of 'x': same sparsity:
			      r@x <- rep.int(TRUE, length(e1@x))
			      for(n in intersect(c("i","j","p"), slotNames(r)))
				  slot(r, n) <- slot(e1, n)
			  }
			  else { ## all FALSE: keep empty 'r' matrix
			      ## but may need a valid 'pointer' slot:
			      if(extends(lClass, "CsparseMatrix"))
				  r@p <- rep.int(0L, 1+ncol(r))
			      else if(extends(lClass, "RsparseMatrix"))
				  r@p <- rep.int(0L, 1+nrow(r))
			  }
		      } else { # some TRUE, FALSE, NA : go via unique 'Tsparse'
			  M <- asTuniq(e1)
			  nCl <- class2(class(M), 'l') # logical Tsparse
			  r <- new(nCl)
			  r@x <- callGeneric(M@x, e2)
			  ## copy "the other slots" (important for "tr"/"sym"):
			  ## "%w/o%" <- function(x,y) x[is.na(match(x, y))]
			  sN <- slotNames(nCl)
			  for(n in sN[is.na(match(sN, "x"))])
			      slot(r, n) <- slot(M, n)
			  if(extends(cl1, "CsparseMatrix"))
			      r <- as(r, "CsparseMatrix")
			  else if(extends(cl1, "RsparseMatrix"))
			      r <- as(r, "RsparseMatrix")
		      }
		  } else {
		      ## non sparse result
		      if(getOption("verbose"))
			  message(sprintf("sparse to dense (%s) coercion in '%s'",
					  lClass, .Generic))
		      rx <- rep.int(r0, d[1]*d[2])
		      rx[1L + encodeInd(non0ind(e1, cl1), nr = d[1])] <- r
		      r <- new(fullCl, x = rx, Dim = d, Dimnames = dimnames(e1))
		  }
	      }
	      r
	  })

## "dMatrix <-> work with 'x' slot {was originally just for "Compare":
##  ------- this is *not* used for 'Arith" which
setMethod("Ops", signature(e1 = "dMatrix", e2 = "dMatrix"),
	  function(e1, e2) {
	      d <- dimCheck(e1,e2)
	      if((dens1 <- is(e1, "denseMatrix"))) gen1 <- is(e1, "generalMatrix")
	      if((dens2 <- is(e2, "denseMatrix"))) gen2 <- is(e2, "generalMatrix")

	      if(dens1 && dens2) { ## both inherit from ddense*

		  if(!gen1) e1 <- as(e1, "dgeMatrix")
		  if(!gen2) e2 <- as(e2, "dgeMatrix")
		  ## now, both are dge {ddense* & general*}

		  r <- callGeneric(e1@x, e2@x)
		  r <- new(paste(.M.kind(r), "geMatrix", sep=''),
			   x = r, Dim = d, Dimnames = dimnames(e1))
	      }
	      else {
		  if(!dens1 && !dens2) {
		      ## both e1 _and_ e2 are sparse
		      ## should not happen since we have <sparse> o <sparse> methods
		      stop("Mistaken intended method dispatch -- please report to ",
			   packageDescription("Matrix")$Author)
		  }
		  ## else
		  if(dens1 && !dens2) ## go to dense
		      r <- callGeneric(e1, as(e2, "denseMatrix"))
		  else ## if(!dens1 && dens2)
		      r <- callGeneric(as(e1, "denseMatrix"), e2)

		  ## criterion "2 * nnz(.) < ." as in sparseDefault() in Matrix()  [./Matrix.R] :
		  if(2 * nnzero(r, na.counted = TRUE) < prod(d))
		      r <- as(r, "sparseMatrix")
	      }
	      r
	  })

### --  I -- dense -----------------------------------------------------------

##-------- originally from ./dgeMatrix.R --------------------

## ----- only work with NAMESPACE importFrom(methods, ..)

setMethod("Arith", signature(e1 = "dgeMatrix", e2 = "dgeMatrix"),
          ##  "+", "-", "*", "^", "%%", "%/%", "/"
	  function(e1, e2) {
	      ## NB:  triangular, symmetric, etc may need own method
	      d1 <- e1@Dim
	      d2 <- e2@Dim
	      eqD <- d1 == d2
	      if (!eqD[1])
		  stop("Matrices must have same number of rows for arithmetic")
	      same.dim <- eqD[2]
	      if (same.dim) {
		  d <- d1
		  dn <- dimNamesCheck(e1, e2)
	      }
	      else { # nrows differ ----> maybe recycling
		  if(d2[2] %% d1[2] == 0) { # nrow(e2) is a multiple
		      e1@x <- rep.int(e1@x, d2[2] %/% d1[2])
		      d <- d2
		      dn <- e2@Dimnames
		  } else if(d1[2] %% d2[2] == 0) { # nrow(e1) is a multiple
		      e2@x <- rep.int(e2@x, d1[2] %/% d2[2])
		      d <- d1
		      dn <- e1@Dimnames
		  } else
		      stop("number of rows are not compatible for ", .Generic)
	      }

	      ## be smart and preserve, e.g., triangular, or symmetric
	      ## but this sucks: For these,
	      ## 'uplo' and 'diag' also must coincide or be dealt with properly

	      ## ==> triangular, symmetric, etc may need own method
	      ##     also since their @x is `non-typical'

##		 if(same.dim) {
##		     if(extends(class(e1), class(e2))) {
##			 e2@x <- callGeneric(e1@x, e2@x)
##			 e2@Dimnames <- dn
##			 e2
##		     }
##		     else if(extends(class(e2), class(e1))) {
##			 e1@x <- callGeneric(e1@x, e2@x)
##			 e1@Dimnames <- dn
##			 e1
##		     }
##		 }
##		 else
		  new("dgeMatrix", Dim = d, Dimnames = dn,
		      x = callGeneric(e1@x, e2@x))
	  })

setMethod("Arith", signature(e1 = "dgeMatrix", e2 = "numeric"),
	  function(e1, e2) {
	      d <- e1@Dim
	      le <- length(e2)
	      if(le == 1 || le == d[1] || prod(d) == le) { # matching dim
		  e1@x <- callGeneric(e1@x, as.vector(e2))
		  e1
	      } else stop ("length of 2nd arg does not match dimension of first")
	  })

setMethod("Arith", signature(e1 = "numeric", e2 = "dgeMatrix"),
	  function(e1, e2) {
	      d <- e2@Dim
	      le <- length(e1)
	      if(le == 1 || le == d[1] || prod(d) == le) { # matching dim
		  e2@x <- callGeneric(as.vector(e1), e2@x)
		  e2
	      } else stop ("length of 1st arg does not match dimension of 2nd")
	  })

##-------- originally from ./ddenseMatrix.R --------------------

## Cheap version: work via "dgeMatrix" and use the group methods there:
## FIXME(?): try to preserve "symmetric", "triangular", ...
setMethod("Arith", signature(e1 = "ddenseMatrix", e2 = "ddenseMatrix"),
          function(e1, e2) callGeneric(as(e1, "dgeMatrix"),
                                       as(e2, "dgeMatrix")))
setMethod("Arith", signature(e1 = "ddenseMatrix", e2 = "numeric"),
          function(e1, e2) callGeneric(as(e1, "dgeMatrix"), e2))
setMethod("Arith", signature(e1 = "numeric", e2 = "ddenseMatrix"),
          function(e1, e2) callGeneric(e1, as(e2, "dgeMatrix")))

## "Logic"
## -------

##-------- originally from ./ldenseMatrix.R --------------------

## These all had "Logic", but now are for all "Ops":

setMethod("Ops", signature(e1="lgeMatrix", e2="lgeMatrix"), .Ops.via.x)

setMethod("Ops", signature(e1="ngeMatrix", e2="ngeMatrix"), .Ops.via.x)

setMethod("Ops", signature(e1="ldenseMatrix", e2="ldenseMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      callGeneric(as(e1, "lgeMatrix"), as(e2, "lgeMatrix"))
	  })

setMethod("Ops", signature(e1="ndenseMatrix", e2="ndenseMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      callGeneric(as(e1, "ngeMatrix"), as(e2, "ngeMatrix"))
	  })


setMethod("Logic", signature(e1 = "logical", e2 = "Matrix"),
          function(e1, e2) callGeneric(e1, as(e2, "lMatrix")))
setMethod("Logic", signature(e1 = "Matrix", e2 = "logical"),
          function(e1, e2) callGeneric(as(e1, "lMatrix"), e2))

### "ANY" here means "any non-Matrix" (since "Ops"(ANY) has already bailout above):
setMethod("Logic", signature(e1 = "ANY", e2 = "Matrix"),
          function(e1, e2) callGeneric(as.logical(e1), as(e2, "lMatrix")))
setMethod("Logic", signature(e1 = "Matrix", e2 = "ANY"),
          function(e1, e2) callGeneric(as(e1, "lMatrix"), as.logical(e2)))

setMethod("Logic", signature(e1 = "logical", e2 = "lMatrix"),
          ## "swap RHS and LHS" and use the method below:
	  function(e1,e2) callGeneric(e2, e1))

## This is modelled after the "Compare" ("dMatrix", "numeric") above
setMethod("Logic", signature(e1 = "lMatrix", e2 = "logical"),
	  function(e1, e2) { ## result will typically be "like" e1:
	      r	 <- callGeneric(e1@x, e2)
	      r0 <- callGeneric(FALSE, e2)
	      d <- e1@Dim
	      ## trivial case first (beware of NA)
	      if(isTRUE(r0 && all(r))) {
		  r <- new(if(d[1] == d[2]) "lsyMatrix" else "lgeMatrix")
		  r@Dim <- d
		  r@Dimnames <- e1@Dimnames
		  r@x <- rep.int(TRUE, prod(d))
	      }
	      else if(extends(cl1 <- getClassDef(cl <- class(e1)),
			      "denseMatrix")) {
		  full <- !isPacked(e1) # << both "dtr" and "dsy" are 'full'
		  if(full || identical(r0, FALSE) || extends(cl1,"symmetricMatrix"))
		      r <- new(cl, x = r, Dim = d, Dimnames = dimnames(e1))
		  else { ## packed matrix with structural 0 and r0 is not FALSE:
		      ##--> result cannot be packed anymore
		      ## [dense & packed & not symmetric ] ==> must be "ltp*" :
		      if(!extends(cl1, "ltpMatrix"))
			  stop("internal bug in \"Logic\" method for \"lMatrix\"; please report")
		      rx <- rep.int(r0, d[1]*d[2])
		      rx[indTri(d[1], upper = (e1@uplo == "U"))] <- r
		      r <- new("lgeMatrix", x = rx,
			       Dim = d, Dimnames = dimnames(e1))
		  }

	      }
	      else { ## lsparseMatrix

		  if(identical(r0, FALSE)) { ## things remain sparse
		      if(!any(is.na(r)) && ((Ar <- all(r)) || !any(r))) {
			  r <- new(cl)
			  r@Dim <- d
			  r@Dimnames <- dimnames(e1)
			  if(Ar) { # 'TRUE' instead of 'x': same sparsity:
			      r@x <- rep.int(TRUE, length(e1@x))
			      for(n in intersect(c("i","j","p"), slotNames(r)))
				  slot(r, n) <- slot(e1, n)
			  }
			  else { ## all FALSE: keep empty 'r' matrix
			      ## but may need a valid 'pointer' slot:
			      if(extends(cl, "CsparseMatrix"))
				  r@p <- rep.int(0L, 1+ncol(r))
			      else if(extends(cl, "RsparseMatrix"))
				  r@p <- rep.int(0L, 1+nrow(r))
			  }
		      } else { # some TRUE, FALSE, NA : go via unique 'Tsparse'
			  M <- asTuniq(e1)
			  nCl <- class(M) # logical Tsparse
			  r <- new(nCl)
			  r@x <- callGeneric(M@x, e2)
			  ## copy "the other slots" (important for "tr"/"sym"):
			  ## "%w/o%" <- function(x,y) x[is.na(match(x, y))]
			  sN <- slotNames(nCl)
			  for(n in sN[is.na(match(sN, "x"))])
			      slot(r, n) <- slot(M, n)
			  if(extends(cl1, "CsparseMatrix"))
			      r <- as(r, "CsparseMatrix")
			  else if(extends(cl1, "RsparseMatrix"))
			      r <- as(r, "RsparseMatrix")
		      }
		  } else {
		      ## non sparse result
		      if(getOption("verbose"))
			  message(sprintf("sparse to dense (%s) coercion in '%s'",
					  cl, .Generic))
		      rx <- rep.int(r0, d[1]*d[2])
		      if(extends(cl1, "triangularMatrix") && e1@diag == "U")
			  r <- c(r, rep.int(callGeneric(1, e2),d[1]))
		      rx[1L + encodeInd(non0ind(e1, cl1), nr = d[1])] <- r
		      fullCl <-
			  if(extends(cl1, "symmetricMatrix")) "lsyMatrix" else "lgeMatrix"
		      r <- new(fullCl, x = rx, Dim = d, Dimnames = dimnames(e1))
		  }
	      }
	      r
	  })


### -- II -- sparse ----------------------------------------------------------

## Have lgC o lgC  and then lgT o lgT  Logic - quite similarly -

## Here's the common functionality
.do.Logic.lsparse <- function(e1,e2, d, dn, isOR, ij1, ij2) {
    ii <- WhichintersectInd(ij1, ij2, nrow=d[1])
    I1 <- ii[[1]] ; has1 <- length(I1) > 0
    I2 <- ii[[2]] ; has2 <- length(I2) > 0

    ## 1) common indices
    i <- ij1[I1, 1]
    j <- ij1[I1, 2]

    if(isOR) { ## i.e. .Generic == "|" i.e. not "&"
	x <- e1@x[I1] | e2@x[I2]

	## 2) "e1 o  FALSE":
	x2 <- if(has1) e1@x[- I1] else e1@x # == callGeneric(e1@x[- I1], FALSE)
	## 3) "0  o e1":
	x3 <- if(has2) e2@x[- I2] else e2@x # == callGeneric(FALSE, e2@x[- I2])
	i <- c(i,
	       if(has1) ij1[-I1, 1] else ij1[, 1],
	       if(has2) ij2[-I2, 1] else ij2[, 1])
	j <- c(j,
	       if(has1) ij1[-I1, 2] else ij1[, 2],
	       if(has2) ij2[-I2, 2] else ij2[, 2])
	x <- c(x, x2, x3)
    } else { ## AND
	x <- e1@x[I1] & e2@x[I2]
    }

    if(any(!x)) { ## drop 'FALSE's
	i <- i[x]
	j <- j[x]
	x <- x[x]
    }
    new("lgTMatrix", Dim = d, Dimnames = dn, i = i, j = j, x = x)
}

setMethod("Logic", signature(e1="lgCMatrix", e2="lgCMatrix"),
	  function(e1, e2) {
	      d <- dimCheck(e1, e2)
	      dn <- dimNamesCheck(e1, e2)
	      ## Very easy case first :
	      if(identical(e1@i, e2@i) && identical(e1@p, e2@p)) {
		  e1@x <- callGeneric(e1@x, e2@x)
		  return(e1)
	      }
	      ## else :

	      .Call(Tsparse_to_Csparse,
		    .do.Logic.lsparse(e1, e2, d = d, dn = dn,
				      isOR = .Generic == "|",
				      ij1 = .Call(compressed_non_0_ij, e1, TRUE),
				      ij2 = .Call(compressed_non_0_ij, e2, TRUE)),
		    FALSE)
	  })

setMethod("Logic", signature(e1="lgTMatrix", e2="lgTMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      dn <- dimNamesCheck(e1, e2)
	      ## Very easy case first :
	      if(identical(e1@i, e2@i) && identical(e1@j, e2@j)) {
		  e1@x <- callGeneric(e1@x, e2@x)
		  return(e1)
	      }
	      ## else :
	      cld <- getClassDef("lgTMatrix")
	      .do.Logic.lsparse(e1, e2, d = d, dn = dn,
				isOR = .Generic == "|",
				ij1 = non0ind(e1, cld),
				ij2 = non0ind(e2, cld))
          })

## Now the other "Ops" for the "lgT" and "lgC" cases:
setMethod("Arith", signature(e1="lgCMatrix", e2="lgCMatrix"),
	  function(e1, e2) callGeneric(as(e1, "dgCMatrix"), as(e2, "dgCMatrix")))
setMethod("Arith", signature(e1="lgTMatrix", e2="lgTMatrix"),
	  function(e1, e2) callGeneric(as(e1, "dgTMatrix"), as(e2, "dgTMatrix")))

## FIXME: These are really too cheap: currently almost all go via dgC*() :
## setMethod("Compare", signature(e1="lgCMatrix", e2="lgCMatrix"),
## setMethod("Compare", signature(e1="lgTMatrix", e2="lgTMatrix"),
## setMethod("Compare", signature(e1="lsparseMatrix", e2="lsparseMatrix"),
## 	  function(e1, e2) callGeneric(as(e1, "dgCMatrix"), as(e2, "dgCMatrix")))
setMethod("Compare", signature(e1="lsparseMatrix", e2="lsparseMatrix"),
	  function(e1, e2) callGeneric(as(e1, "CsparseMatrix"),
				       as(e2, "CsparseMatrix")))

## setMethod("Compare", signature(e1="lgTMatrix", e2="lgTMatrix"), ## coerce to Csparse

## 	  function(e1, e2) callGeneric(as(e1, "dgCMatrix"), as(e2, "dgCMatrix")))


##--- Sparse ... ----------

.Ops.lsparse <- function(e1, e2) {
    ## of course, this is still suboptimal, e.g. if both are "lsT" or "lsC"
    if(is(e1, "TsparseMatrix"))
        callGeneric(as(e1, "lgTMatrix"), as(e2, "lgTMatrix"))
    else
	callGeneric(as(e1, "lgCMatrix"), as(e2, "lgCMatrix"))
}

setMethod("Ops", signature(e1="lsparseMatrix", e2="lsparseMatrix"),
	  .Ops.lsparse)

setMethod("Logic", signature(e1="lsparseMatrix", e2="ldenseMatrix"),
	  .Ops.lsparse)

setMethod("Logic", signature(e1="ldenseMatrix", e2="lsparseMatrix"),
	  function(e1, e2) {
	      if(is(e2, "TsparseMatrix"))
		  callGeneric(as(e1, "lgTMatrix"), as(e2, "lgTMatrix"))
	      else
		  callGeneric(as(e1, "lgCMatrix"), as(e2, "lgCMatrix"))
	  })



## FIXME: also want (symmetric o symmetric) , (triangular o triangular)
## -----
setMethod("Arith", signature(e1 = "dsCMatrix", e2 = "dsCMatrix"),
	  function(e1, e2) {
	      if(getOption("verbose"))
		  message("suboptimal implementation of sparse 'symm. o symm.'")
	      forceSymmetric(callGeneric(as(e1, "dgCMatrix"), as(e2, "dgCMatrix")))
	  })


##-------- originally from ./dgCMatrix.R --------------------

.Arith.Csparse <- function(e1, e2, Generic, triangular = FALSE)
{
    ## Generic is one of  "+", "-", "*", "^", "%%", "%/%", "/"

    ## triangular:  TRUE  iff e1,e2 are triangular  _and_  e1@uplo == e2@uplo
    d <- dimCheck(e1, e2)
    dn <- dimNamesCheck(e1, e2)
    if(triangular) {
	## need these for the 'x' slots in any case
	if (e1@diag == "U") e1 <- .Call(Csparse_diagU2N, e1)
	if (e2@diag == "U") e2 <- .Call(Csparse_diagU2N, e2)
        ## slightly more efficient than non0.i() or non0ind():
	ij1 <- .Call(compressed_non_0_ij, e1, isC=TRUE)
	ij2 <- .Call(compressed_non_0_ij, e2, isC=TRUE)

	newTMat <- function(i,j,x)
	    new("dtTMatrix", Dim = d, Dimnames = dn, i = i, j = j, x = x,
		uplo = e1@uplo)
	dmat <- "dtrMatrix"
    } else {
	ij1 <- non0ind(e1)
	ij2 <- non0ind(e2)

	newTMat <- function(i,j,x)
	    new("dgTMatrix", Dim = d, Dimnames = dn, i = i, j = j, x = x)
	dmat <- "dgeMatrix"
    }

    switch(Generic,
	   "+" = , "-" = {
	       ## special "T" convention: repeated entries are *summed*
	       .Call(Tsparse_to_Csparse,
		     newTMat(i = c(ij1[,1], ij2[,1]),
			     j = c(ij1[,2], ij2[,2]),
			     x = if(Generic == "+") c(e1@x, e2@x) else c(e1@x, - e2@x)),
		     triangular)
	   },

	   "*" =
       { ##  X * 0 == 0 * X == 0 --> keep common non-0
	   ii <- WhichintersectInd(ij1, ij2, nrow=d[1])
	   ij <- ij1[ii[[1]], , drop = FALSE]
	   .Call(Tsparse_to_Csparse,
		 newTMat(i = ij[,1],
			 j = ij[,2],
			 x = e1@x[ii[[1]]] * e2@x[ii[[2]]]),
		 triangular)
       },

	   "^" =
       {
	   ii <- WhichintersectInd(ij1, ij2, nrow=d[1])
	   ## 3 cases:
	   ## 1) X^0 := 1  (even for X=0) ==> dense
	   ## 2) 0^Y := 0  for Y != 0	      =====
	   ## 3) x^y :

	   ## FIXME:	dgeM[cbind(i,j)] <- V  is not yet possible
	   ##	    nor dgeM[ i_vect   ] <- V
	   ## r <- as(e2, "dgeMatrix")
	   ## ...
	   r <- as(e2, "matrix")
	   Yis0 <- is0(r)
	   r[complementInd(ij1, dim=d)] <- 0 ## 2)
	   r[1L + ij2[ii[[2]], , drop=FALSE]] <-
	       e1@x[ii[[1]]] ^ e2@x[ii[[2]]] ## 3)
	   r[Yis0] <- 1			     ## 1)
	   as(r, dmat)
       },

	   "%%" = , "%/%" = , "/" = ## 0 op 0	 |-> NaN => dense
	   get(Generic)(as(e1, dmat), e2)

	   )# end{switch(..)}
}

setMethod("Arith", signature(e1 = "dgCMatrix", e2 = "dgCMatrix"),
	  function(e1,e2) .Arith.Csparse(e1,e2, .Generic))

setMethod("Arith", signature(e1 = "dtCMatrix", e2 = "dtCMatrix"),
	  function(e1, e2) {
	      U1 <- e1@uplo
	      isTri <- U1 == e2@uplo # will the result definitely be triangular?
	      if(isTri) {
		  .Arith.Csparse(e1,e2, .Generic, triangular=TRUE)
	      }
	      else { ## lowerTri  o  upperTri: |--> "all 0" {often} -- FIXME?
		  callGeneric(as(e1, "dgCMatrix"), as(e2, "dgCMatrix"))
	      }
	  })



setMethod("Arith", signature(e1 = "dgCMatrix", e2 = "numeric"),
	  function(e1, e2) {
	      if((l2 <- length(e2)) == 0)
		  stop("<Matrix>",.Generic,"numeric(0) is undefined")
	      f0 <- callGeneric(0, e2)
	      if(all(is0(f0))) { # remain sparse
		  if(l2 > 1) { # length(e2) > 1 : "recycle" e2 "carefully"
		      n <- prod(d <- dim(e1))
		      if(n < l2)
			  stop("<Matrix>",.Generic,"numeric(<too-long>)")
		      if(n %% l2 != 0) ## identical warning as in main/arithmetic.c
			  warning(gettextf("longer object length\n\tis not a multiple of shorter object length"))
		      ## TODO: construction of [1L + in0 %%l2] via one .Call()
		      ## 0-based indices:
		      in0 <- encodeInd(.Call(compressed_non_0_ij, e1, TRUE), d[1])
		      e2 <- e2[1L + in0 %% l2]
		  }
		  e1@x <- callGeneric(e1@x, e2)
		  e1
	      }
	      else { ## non-sparse, since '0 o e2' is not (all) 0
		  r <- as(e1, "matrix")
		  if(length(e2) == 1) {
		      r[] <- f0
		      r[non0ind(e1, getClassDef("dgCMatrix")) + 1L] <-
                          callGeneric(e1@x, e2)
		      as(r, "dgeMatrix")
		  } else {
		      as(callGeneric(r, e2), "dgeMatrix")
		  }
	      }})


setMethod("Arith", signature(e1 = "numeric", e2 = "dgCMatrix"),
	  function(e1, e2) {
	      if((l1 <- length(e1)) == 0)
		  stop("numeric(0)",.Generic,"<Matrix> is undefined")
	      f0 <- callGeneric(e1, 0)
	      if(all(is0(f0))) { # remain sparse
		  if(l1 > 1) { # length(e1) > 1 : "recycle" e1 "carefully"
		      n <- prod(d <- dim(e2))
		      if(n < l1)
			  stop("numeric(<too-long>)",.Generic,"<Matrix>")
		      if(n %% l1 != 0) ## identical warning as in main/arithmetic.c
			  warning(gettextf("longer object length\n\tis not a multiple of shorter object length"))
		      ## TODO: construction of [1L + in0 %% l1] via one .Call()
		      ## 0-based indices:
		      in0 <- encodeInd(.Call(compressed_non_0_ij, e2, TRUE), d[1])
		      e1 <- e1[1L + in0 %% l1]
		  }
		  e2@x <- callGeneric(e1, e2@x)
		  e2
	      }
	      else { ## non-sparse, since '0 o e2' is not (all) 0
		  r <- as(e2, "matrix")
		  if(length(e1) == 1) {
		      r[] <- f0
		      r[non0ind(e2, getClassDef("dgCMatrix")) + 1L] <-
                        callGeneric(e1, e2@x)
		      as(r, "dgeMatrix")
		  } else {
		      as(callGeneric(e1, r), "dgeMatrix")
		  }
	      }})

##-------- originally from ./Csparse.R --------------------
## TODO : Consider going a level up, and do this for all "Ops"

setMethod("Arith", signature(e1 = "CsparseMatrix", e2 = "CsparseMatrix"),
	  function(e1, e2) {
	      ## go via "symmetric" if both are symmetric, etc...
	      s1 <- .M.shape(e1, c1 <- getClassDef(class(e1)))
	      s2 <- .M.shape(e2, c2 <- getClassDef(class(e2)))
	      viaCl <- paste("d", if(s1 == s2) s1 else "g", "CMatrix", sep='')
	      callGeneric(as(as(e1, "dMatrix"), viaCl),
			  as(as(e2, "dMatrix"), viaCl))
          })

setMethod("Arith", signature(e1 = "CsparseMatrix", e2 = "numeric"),
	  function(e1, e2) {
	      if(length(e2) == 1) { ## e.g.,  Mat ^ a
		  f0 <- callGeneric(0, e2)
		  if(is0(f0)) { # remain sparse, symm., tri.,...
                      e1 <- diagU2N(as(e1, "dMatrix"))
		      e1@x <- callGeneric(e1@x, e2)
		      return(e1)
		  }
	      }
	      ## all other (potentially non-sparse) cases: give up symm, tri,..
	      callGeneric(as(e1, paste(.M.kind(e1), "gCMatrix", sep='')), e2)
	  })

## The same,  e1 <-> e2 :
setMethod("Arith", signature(e1 = "numeric", e2 = "CsparseMatrix"),
	  function(e1, e2) {
	      if(length(e1) == 1) {
		  f0 <- callGeneric(e1, 0)
		  if(is0(f0)) {
                      e2 <- diagU2N(as(e2, "dMatrix"))
		      e2@x <- callGeneric(e1, e2@x)
		      return(e2)
		  }
	      }
	      callGeneric(e1, as(e2, paste(.M.kind(e2), "gCMatrix", sep='')))
	  })

setMethod("Compare", signature(e1 = "CsparseMatrix", e2 = "CsparseMatrix"),
	  function(e1, e2) {
	      d <- dimCheck(e1,e2)

	      ## How do the "0" or "FALSE" entries compare?
	      ## Depends if we have an "EQuality RELation" or not:
	      EQrel <- switch(.Generic,
			      "==" =, "<=" =, ">=" = TRUE,
			      "!=" =, "<"  =, ">"  = FALSE)
	      if(EQrel) {
		  ## The (0 op 0) or  (FALSE op FALSE) comparison gives TRUE
		  ## -> result becomes *dense*; the following may be suboptimal
		  return( callGeneric(as(e1, "denseMatrix"),
				      as(e2, "denseMatrix")))
	      }

	      ## else: INequality:   0 op 0 gives FALSE ---> remain sparse!

              cD1 <- getClassDef(class(e1))
              cD2 <- getClassDef(class(e2))

	      ## NB non-diagonalMatrix := Union{ general, symmetric, triangular}
	      gen1 <- extends(cD1, "generalMatrix")
	      gen2 <- extends(cD2, "generalMatrix")
	      sym1 <- !gen1 && extends(cD1, "symmetricMatrix")
	      sym2 <- !gen2 && extends(cD2, "symmetricMatrix")
	      tri1 <- !gen1 && !sym1
	      tri2 <- !gen2 && !sym2
	      G <- gen1 && gen2
	      S <- sym1 && sym2 && e1@uplo == e2@uplo
	      T <- tri1 && tri2 && e1@uplo == e2@uplo

	      if(T && e1@diag != e2@diag) {
		  ## one is "U" the other "N"
		  if(e1@diag == "U")
		      e1 <- diagU2N(e1)
		  else ## (e2@diag == "U"
		      e2 <- diagU2N(e2)
	      }
	      else if(!G && !S && !T) {
                  ## e.g. one symmetric, one general
                  ## coerce to generalMatrix and go :
		  if(!gen1) e1 <- as(e1, "generalMatrix", strict = FALSE)
		  if(!gen2) e2 <- as(e2, "generalMatrix", strict = FALSE)
	      }

	      dn <- dimNamesCheck(e1, e2)
              ## the result object:
	      newC <- sub("^.", "l", class(e1))
	      r <- new(newC)

              ## Very easy case:
              if(identical(e1@i, e2@i) && identical(e1@p, e2@p)) {
		  if(extends(cD2, "nMatrix")) {
		      ## non-equality of identical pattern matrices: all FALSE
		  } else { # have 'x' slot
		      r@x <- callGeneric(e1@x, e2@x)
		      ## and all others are  '0 op 0' which give FALSE
		      r@i <- e1@i
		      r@p <- e1@p
		  }
		  r@Dim <- d
		  r@Dimnames <- dn
		  r
              }
              else {
                  ## now the 'x' slots ``match'' insofar as they are for the
                  ## same "space" (triangle for tri* and symm*; else rectangle)

                  ## not non0ind() which gives more;
                  ## want only those which correspond to 'x' slot
                  ij1 <- .Call(compressed_non_0_ij, e1, TRUE)
                  ij2 <- .Call(compressed_non_0_ij, e2, TRUE)
                  ii <- WhichintersectInd(ij1, ij2, nrow=d[1])
                  I1 <- ii[[1]]
                  I2 <- ii[[2]]

                  ## 1) common
                  x <- callGeneric(e1@x[I1],
                                   e2@x[I2])
                  ## 2) "e1 o  0":
                  x2 <- callGeneric(e1@x[- I1], 0)
                  ## 3) "0  o e1":
                  x3 <- callGeneric(0, e2@x[- I2])

		  i <- c(ij1[I1, 1], ij1[-I1, 1], ij2[-I2, 1])
		  j <- c(ij1[I1, 2], ij1[-I1, 2], ij2[-I2, 2])
		  x <- c( x,	    x2,		 x3)
		  if(any(!x)) { # drop 'FALSE's
		      i <- i[x]
		      j <- j[x]
		      x <- x[x]
		  }
		  as(new("lgTMatrix", Dim = d, Dimnames = dn,
			 i = i, j = j, x = x),
		     "lgCMatrix")
              }
	  })


##-------- originally from ./sparseMatrix.R --------------------

## "Arith" short cuts / exceptions
setMethod("-", signature(e1 = "sparseMatrix", e2 = "missing"),
	  function(e1, e2) { e1 <- diagU2N(e1); e1@x <- -e1@x; e1 })
## with the following exceptions:
setMethod("-", signature(e1 = "nsparseMatrix", e2 = "missing"),
          function(e1,e2) callGeneric(as(as(e1, "dMatrix"), "dgCMatrix")))
setMethod("-", signature(e1 = "pMatrix", e2 = "missing"),
          function(e1,e2) callGeneric(as(e1, "ngTMatrix")))

## Group method  "Arith"

## have CsparseMatrix methods above
## which may preserve "symmetric", "triangular" -- simply defer to those:

## these where 'Arith', now generalized:
setMethod("Ops", signature(e1 = "sparseMatrix", e2 = "sparseMatrix"),
	  function(e1, e2) callGeneric(as(e1, "CsparseMatrix"),
				       as(e2, "CsparseMatrix")))
setMethod("Ops", signature(e1 = "sparseMatrix", e2 = "numeric"),
	  function(e1, e2) callGeneric(as(e1, "CsparseMatrix"), e2))
setMethod("Ops", signature(e1 = "numeric", e2 = "sparseMatrix"),
	  function(e1, e2) callGeneric(e1, as(e2, "CsparseMatrix")))

## setMethod("Compare", signature(e1 = "sparseMatrix", e2 = "sparseMatrix"),
## 	  function(e1, e2) callGeneric(as(e1, "CsparseMatrix"),
## 				       as(e2, "CsparseMatrix")))

###-------- sparseVector -------------
###-------- ============ -------------

## Catch all ( ==> better error message than default):
setMethod("Ops", signature(e1 = "sparseVector", e2 = "ANY"),
          function(e1, e2) .bail.out.2(.Generic, class(e1), class(e2)))
setMethod("Ops", signature(e1 = "ANY", e2 = "sparseVector"),
          function(e1, e2) .bail.out.2(.Generic, class(e1), class(e2)))

## 1)  spVec  o  (sp)Vec : -------------

## FIXME:
##   1. <spVec>  o  <scalar>   should  *NOT* go via  <spV> o <spV>
##   2. <spVec>  o  <non-NA numeric>  should also happen directly and
##                                    |-> sparse for o = {'*', "/", '&&'

setMethod("Ops", signature(e1 = "sparseVector", e2 = "atomicVector"),
	  function(e1, e2) callGeneric(e1, as(e2, "sparseVector")))
setMethod("Ops", signature(e1 = "atomicVector", e2 = "sparseVector"),
	  function(e1, e2) callGeneric(as(e1, "sparseVector"), e2))

setMethod("Arith", signature(e1 = "sparseVector", e2 = "sparseVector"),
          function(e1, e2) callGeneric(as(e1, "dsparseVector"),
                                       as(e2, "dsparseVector")))
setMethod("Logic", signature(e1 = "sparseVector", e2 = "sparseVector"),
          function(e1, e2) callGeneric(as(e1, "lsparseVector"),
                                       as(e2, "lsparseVector")))

setMethod("Arith", signature(e1 = "dsparseVector", e2 = "dsparseVector"),

          ##  "+", "-", "*", "^", "%%", "%/%", "/"

          function(e1, e2) {
              n1 <- e1@length
              n2 <- e2@length
              if(n1 != n2) {
                  if(n1 < n2) {
                      n <- n1 ; N <- n2
                  } else {
                      n <- n2 ; N <- n1
                  }
                  if(N %% n != 0) ## require this here, for conveniense
                      ## for regular vectors, this is only a warning:
                      stop("longer object length\n\t",
                           "is not a multiple of shorter object length")
                  if(n == 1) { # simple case, do not really recycle
                      if(n1 < n2) return(callGeneric(sp2vec(e1, "double"), e2))
                      else        return(callGeneric(e1, sp2vec(e2, "double")))
                  }
                  ## else : 2 <= n < N --- recycle the shorter one
                  q <- N %/% n
                  if(n1 < n2) {
                      e1@i <- rep.int(e1@i, q)
                      e1@x <- rep.int(e1@x, q)
                  } else {
                      e2@i <- rep.int(e2@i, q)
                      e2@x <- rep.int(e2@x, q)
                  }
              } else { ## n1 == n2
                  N <- n1
              }
              r <- new("dsparseVector", length = N)
	      switch(.Generic,
		     "+" = , "-" =  ## X +- 0 == 0 +- X == X

		     ii <- union(e1@i, e2@i),

		     "*" =
		 { ##  X * 0 == 0 * X == 0 --> keep common non-0
		     ii <- intersect(e1@i, e2@i)

		 },

		     "^" =
		 {

		     ii <- intersect(e1@i, e2@i)
		     ## 3 cases:
		     ## 1) X^0 := 1  (even for X=0) ==> dense
		     ## 2) 0^Y := 0  for Y != 0		=====
		     ## 3) x^y :

		     ## FIXME:	dgeM[cbind(i,j)] <- V  is not yet possible
		     ##	    nor dgeM[ i_vect   ] <- V
		     ## r <- as(e2, "dgeMatrix")
		     ## ...
		 },

		     "%%" = , "%/%" = , "/" = ## 0 op 0	 |-> NaN => dense

		     )

              .bail.out.2(.Generic, class(e1), class(e2))
          })

## "Arith"  exception (shortcut)
setMethod("-", signature(e1 = "dsparseVector", e2 = "missing"),
          function(e1) { e1@x <- -e1@x ; e1 })


setMethod("Logic", signature(e1 = "lsparseVector", e2 = "lsparseVector"),
          function(e1, e2) {
              .bail.out.2(.Generic, class(e1), class(e2))
          })

## 2)  spVec  o  [Mm]atrix : -------------
