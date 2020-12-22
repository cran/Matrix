####--- All "Math" and "Math2" group methods for all Matrix classes (incl sparseVector) ------
####	     ====	=====

## "Design-bug":  log(x, base)  has *two* arguments // ditto for  "trunc()" !!
## ---> need "log" methods "everywhere to catch 2-arg case !


###--------- Csparse

Math.vecGenerics <- grep("^cum", getGroupMembers("Math"), value=TRUE)
## "cummax" .. "cumsum" : work on full *vector* and return vector also for matrix input
setMethod("Math", "CsparseMatrix", function(x)
{
    if(.Generic %nin% Math.vecGenerics && is0(callGeneric(0.))) {
        ## sparseness, symm., triang.,... preserved
        cl <- class(x)
        has.x <- !extends(cl, "nsparseMatrix")
        ## has.x  <==> *not* nonzero-pattern == "nMatrix"
        if(has.x) {
            type <- storage.mode(x@x)
            r <- callGeneric(x@x)
        } else { ## nsparseMatrix
            type <- ""
            r <- rep.int(as.double(callGeneric(TRUE)),
                         switch(.sp.class(cl),
                                CsparseMatrix = length(x@i),
                                TsparseMatrix = length(x@i),
                                RsparseMatrix = length(x@j)))
        }
        if(type == storage.mode(r)) {
            x@x <- r
            x
        } else { ## e.g. abs( <lgC> ) --> integer Csparse
            ## FIXME: when we have 'i*' classes, use them here:
            rx <- new(sub("^.", "d", MatrixClass(cl)))
            rx@x <- as.double(r)
            ## result is "same"
            sNams <- slotNames(cl)
            for(nm in sNams[sNams != "x"])
                slot(rx, nm) <- slot(x, nm)
            rx
        }
    } else { ## no sparseness (or no matrix!); C2dense() returns *numeric*
        callGeneric(C2dense(x))
    }
}) ## {Math}

setMethod("log", "CsparseMatrix", function(x, base = exp(1)) log(C2dense(x), base))

###--------- ddenseMatrix

##' Used for  dt[rp]Matrix, ds[yp]Matrix (and subclasses, e.g. dpo*(), cor*() !):
##' as dgeMatrix has direct method:
setMethod("Math", "ddenseMatrix", function(x)
    {
	if(.Generic %in% Math.vecGenerics) # vector result
	    callGeneric(as(x,"dgeMatrix")@x)
	else if(is(x, "symmetricMatrix")) { ## -> result symmetric: keeps class
            cld <- getClassDef(cl <- class(x))
	    if((po <- extends(cld, "dpoMatrix")) || extends(cld, "dppMatrix")) { # result is *not* pos.def!
		x <- as(x, if(po) "dsyMatrix" else "dspMatrix")
            }
            ## "symmetricMatrix" has 'factors' slot:
            if(!is.null(x@factors)) x@factors <- list()
            x@x <- callGeneric(x@x)
            x
	}
	else { ## triangularMatrix (no need for testing), includes, e.g. "corMatrix"!
	    ## if(is0(f0 <- callGeneric(0.))) { ## -> result remains triangular
	    if(is0(callGeneric(0.))) { ## -> result remains triangular
		cld <- getClassDef(cl <- class(x))
                if(extends(cld, "triangularMatrix")) {
                    if((isF <- extends(cld, "MatrixFactorization")) || extends(cld, "corMatrix")) {
                        x <- as(x, if(isF && .isPacked(x)) "dtpMatrix" else "dtrMatrix")
                    }
                } else
                    if(inherits(x, "compMatrix")) # has 'factors' slot
                        if(!is.null(x@factors)) x@factors <- list()
		x@x <- callGeneric(x@x)
		x
	    }
	    else {
                if(inherits(x, "compMatrix")) # has 'factors' slot
                    if(!is.null(x@factors)) x@factors <- list()
                ## result is general: *could* use f0 <- callGeneric(0.) for the whole 0-triangle,
		## but this is much easier:
		callGeneric(as(x,"dgeMatrix"))
	    }
	}
    })

## "log" with *two* arguments
setMethod("log", "ddenseMatrix", function(x, base = exp(1))
{
    if(is(x, "symmetricMatrix")) { ## -> result symmetric: keeps class
        cld <- getClassDef(class(x))
        if((po <- extends(cld, "dpoMatrix")) || extends(cld, "dppMatrix")) { # result is *not* pos.def!
            x <- as(x, if(po) "dsyMatrix" else "dspMatrix")
        }
        ## "symmetricMatrix" has 'factors' slot:
        if(!is.null(x@factors)) x@factors <- list()
        x@x <- log(x@x, base)
        x
    }
    else { ## triangularMatrix or generalMatrix, includes, e.g. "corMatrix"!
        if(inherits(x, "compMatrix")) # has 'factors' slot
            if(!is.null(x@factors)) x@factors <- list()
        ## result is general: *could* use  -Inf   for the whole 0-triangle,
        ## but this is much easier:
        log(as(x,"dgeMatrix"), base)
    }
})


###--------- denseMatrix

## FIXME: Once we have integer (idense..),  sign(), abs(.) may need different:
setMethod("Math", signature(x = "denseMatrix"),
	  function(x) callGeneric(as(x, "dMatrix"))) # -> "ddenseMatrix" above
setMethod("log", "denseMatrix", function(x, base = exp(1)) log(as(x, "dMatrix"), base))

###--------- dgeMatrix

setMethod("Math", signature(x = "dgeMatrix"),
	  function(x) {
	      if(.Generic %in% Math.vecGenerics)
		  callGeneric(x@x)
	      else {
		  x@x <- callGeneric(x@x)
		  x
	      }
	  })
setMethod("log", "dgeMatrix", function(x, base = exp(1)) {
    x@x <- log(x@x, base)
    x
})

###--------- diagMatrix

## Till 2014-08-04, went via "dtC" (triangular)
setMethod("Math", signature(x = "diagonalMatrix"),
	  function(x) {
	      if(.Generic %in% Math.vecGenerics) # vector result
		  callGeneric(.diag2mat(x))
	      ## else if(is0(f0 <- callGeneric(0.))) { ## result remains diagonal
	      else if(is0(callGeneric(0.))) { ## result remains diagonal
		  cl <- class(x)
		  if(!extends(cl, "ddiMatrix"))
		      cl <- class(x <- as(x, "dMatrix"))
		  ##d type <- storage.mode(x@x)
                  if(x@diag == "U") {
		      ##d if((f1 <- callGeneric(as1(mod=type))) == 1 && type == "double")
		      if((f1 <- callGeneric(1.)) == 1)
			  return(x) # [ddi] as f(0) = 0, f(1) = 1
		      else {
			  n <- x@Dim[1]
			  return( Diagonal(n=n, x = rep.int(f1, n)) )
		      }
                  }
                  r <- callGeneric(x@x)
		  ##d if(type == storage.mode(r)) {
		      x@x <- r
		      x
		  ##d } else { ## e.g. abs( <lgC> ) --> integer Csparse
		  ##d     ## FIXME: when we have 'i*' classes, use them here:
		  ##d     rx <- new(sub("^.", "d", cl))
		  ##d     rx@x <- as.double(r)
		  ##d     ## result is "same"
		  ##d     sNams <- slotNames(cl)
		  ##d     for(nm in sNams[sNams != "x"])
		  ##d         slot(rx, nm) <- slot(x, nm)
		  ##d     rx
		  ##d }
	      } else { ## no sparseness, i.e., no diagonal, but still symmetric:
		  ## FIXME: gain efficiency by reusing f0  for *all* off-diagonal entries!
		  callGeneric(as(as(as(.diag2sT(x), "dMatrix"), "denseMatrix"), "dspMatrix"))
	      }
	  }) ## {Math}

setMethod("log", "diagonalMatrix", function(x, base = exp(1)) {
    ## no sparseness, i.e., no diagonal, but still symmetric:
    r <- as(as(as(.diag2sT(x), "dMatrix"), "denseMatrix"), "dspMatrix")
    diag(r) <- if(x@diag == "U") 0 else log(x@x, base)
    ## Assign log(0, <base>) == -Inf  to all off-diagonal elements;
    ## indices depend crucially on  uplo = "U" / "L" :
    n <- x@Dim[[1L]]
    if(n >= 1L) {
        k <- seq_len(n)
	i <- k*(k+1)/2 # as r@uplo == "U"
        ## } else { # uplo == "L"
        ##     cumsum(c(1, if(n>1) n:2))
        ## }
	r@x[-i] <- -Inf # = log(0, <base>)
    }
    r
})


## NB: "Math2" (round, signif) for diagMatrix is perfectly via "dMatrix"


###--------- dMatrix

## Use these as "catch-all" -- more specific methods are for sub-classes (sparse)

setMethod("Math2", signature(x = "dMatrix"),
          ## Assume that  Generic(u, k) |--> u for u in {0,1}
          ## which is true for round(), signif() ==> all structure maintained
	  function(x, digits) {
              x@x <- callGeneric(x@x, digits = digits)
              x
          })
## the same, first coercing to "dMatrix":
setMethod("Math2", signature(x = "Matrix"),
	  function(x, digits) {
	      x <- as(x, "dMatrix")
	      x@x <- callGeneric(x@x, digits = digits)
	      x
	  })


###--------- sparseMatrix

setMethod("Math", signature(x = "sparseMatrix"),
	  function(x) callGeneric(as(x, "CsparseMatrix")))

setMethod("log", "sparseMatrix", function(x, base = exp(1)) log(as(x, "CsparseMatrix"), base))


###--------- sparseVector

setMethod("Math", signature(x = "sparseVector"),
	  function(x) {
	      if(.Generic %nin% Math.vecGenerics && is0(callGeneric(0.))) {
		  ## sparseness preserved
		  cld <- getClassDef(class(x))
		  kind <- .M.kindC(cld)# "d", "n", "l", "i", "z", ...
		  has.x <- kind != "n"
		  if(has.x) {
		      rx <- callGeneric(x@x)
		      if(kind == "d") {
			  x@x <- rx
			  x
		      }
		      else {
			  new("dsparseVector", x = rx, i = x@i, length = x@length)
		      }
		  } else { ## kind == "n"
		      new("dsparseVector", x = rep.int(callGeneric(1), length(x@i)),
			  i = x@i, length = x@length)
		  }
	      } else { ## dense
		  callGeneric(sp2vec(x))
	      }
	  })

setMethod("log", "sparseVector", function(x, base = exp(1)) log(sp2vec(x), base))

setMethod("Math2", signature(x = "dsparseVector"),
          ## Assume that  Generic(u, k) |--> u for u in {0,1}
          ## which is true for round(), signif() ==> all structure maintained
	  function(x, digits) {
              x@x <- callGeneric(x@x, digits = digits)
              x
          })
## the same, first coercing to "dsparseVector":
setMethod("Math2", signature(x = "sparseVector"),
	  function(x, digits) {
	      x <- as(x, "dsparseVector")
	      x@x <- callGeneric(x@x, digits = digits)
	      x
	  })

