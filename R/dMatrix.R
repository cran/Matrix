### Define Methods that can be inherited for all subclasses

##-> "dMatrix" <--> "lMatrix"   ---> ./lMatrix.R

## these two are parallel to "n <-> l" in the above :
setAs("nMatrix", "dMatrix",
      function(from) {
	  cld <- getClassDef(cl <- class(from))
	  isSp <- extends(cld, "sparseMatrix")
	  ## faster(not "nicer"): any(substr(cl,3,3) == c("C","T","R"))
	  sNams <- slotNames(cld)
	  r <- copyClass(from, sub("^n", "d", cl),
			 if(isSp) sNams else sNams[sNams != "x"])
	  r@x <- if(isSp) rep.int(1., nnzSparse(from)) else as.double(from@x)
	  r
      })

## NOTE: This is *VERY* parallel to  ("lMatrix" -> "nMatrix") in ./lMatrix.R :
setAs("dMatrix", "nMatrix",
      function(from) {
	  if(any(is.na(from@x)))
	      stop("\"dMatrix\" object with NAs cannot be coerced to \"nMatrix\"")
	  ## i.e. from@x are only TRUE (or FALSE in dense case)
	  cld <- getClassDef(cl <- class(from))
	  if(extends(cld, "diagonalMatrix")) { # have no "ndi*" etc class
	      cl <- class(from <- as(from, "sparseMatrix"))
	      isSp <- TRUE
	  } else {
	      isSp <- extends(cld, "sparseMatrix")
	      if(isSp && any(from@x == 0)) {
		  from <- drop0(from, cld)
		  if(cl != (c. <- class(from)))
		      cld <- getClassDef(cl <- c.)
	      }
	  }
	  sNams <- slotNames(cld)
	  copyClass(from, sub("^d", "n", cl),
		    if(isSp) sNams[sNams != "x"] else sNams)
      })


## Methods for operations where one argument is integer
## No longer made use of (and confusing hence) since R version 2.1.0
## where "integer" goes as part of "numeric"

## Note: Use as.matrix() {not directly array()} :
##  1) to ensure consistency with "numeric" (non-matrix)
##  2) names -> dimnames {potentially}
## setMethod("%*%", signature(x = "dMatrix", y = "integer"),
##           function(x, y) callGeneric(x, as.numeric(y)))

## setMethod("%*%", signature(x = "integer", y = "dMatrix"),
##           function(x, y) callGeneric(as.numeric(x), y))

## setMethod("crossprod", signature(x = "dMatrix", y = "integer"),
##           function(x, y = NULL) callGeneric(x, as.numeric(y)))

## setMethod("crossprod", signature(x = "integer", y = "dMatrix"),
##           function(x, y = NULL) callGeneric(as.numeric(x), y))

## setMethod("solve", signature(a = "dMatrix", b = "integer"),
##           function(a, b, ...) callGeneric(a, as.numeric(b)))


## Group Methods, see ?Arith (e.g.)
## -----
## >>> More specific methods for sub-classes (sparse), use these as "catch-all":

## the non-Ops ones :
setMethod("Math2",
          ## Assume that  Generic(u, k) |--> u for u in {0,1}
          ## which is true for round(), signif() ==> all structure maintained
          signature(x = "dMatrix"),
	  function(x, digits) {
              x@x <- callGeneric(x@x, digits = digits)
              x
          })

## at installation time:
## "max" "min" "range"  "prod" "sum"   "any" "all" :
summGenerics <- getGroupMembers("Summary")
## w/o "prod" & "sum":
summGener1 <- summGenerics[match(summGenerics, c("prod","sum"), 0) == 0]

## [also needs extra work in ./AllGeneric.R ] :
setMethod("Summary", signature(x = "ddenseMatrix", na.rm = "ANY"),
	  function(x, ..., na.rm) {
	      d <- x@Dim
	      if(any(d == 0)) return(callGeneric(numeric(0), ..., na.rm=na.rm))
	      clx <- getClassDef(class(x))
	      if(extends(clx, "generalMatrix"))
		  callGeneric(x@x, ..., na.rm = na.rm)
	      else if(extends(clx, "symmetricMatrix")) { # incl packed, pos.def.
		  if(.Generic %in% summGener1) {
		      callGeneric(if (length(x@x) < prod(d)) x@x
				  else x@x[indTri(d[1], upper= x@uplo == "U",
						  diag= TRUE)],
				  ..., na.rm = na.rm)
		  } else callGeneric(as(x, "dgeMatrix")@x, ..., na.rm = na.rm)
	      }
	      else { ## triangular , packed
		  if(.Generic %in% summGener1)
		      callGeneric(x@x, 0, if(x@diag == "U") 1, ..., na.rm = na.rm)
		  else callGeneric(as(x, "dgeMatrix")@x, ..., na.rm = na.rm)
	      }
	  })

setMethod("Summary", signature(x = "dsparseMatrix", na.rm = "ANY"),
	  function(x, ..., na.rm)
      {
	  ne <- prod(d <- dim(x))
	  if(ne == 0) return(callGeneric(numeric(0), ..., na.rm=na.rm))
	  l.x <- length(x@x)
	  if(l.x < ne) {
	      clx <- getClassDef(class(x))
	      if(extends(clx, "symmetricMatrix") && l.x == choose(d[1]+1,2)) {
		  ## fully non-zero - very rare!
		  callGeneric((if(.Generic %in% summGener1) x else
			       as(x, "generalMatrix"))@x,
			      ..., na.rm = na.rm)
	      }
	      else { ## has at least one structural 0 (e.g. triangular)	 --normal case--
		  if(.Generic == "prod") 0 else
		  callGeneric((if(.Generic %in% summGener1) diagU2N(x) else
			       as(x, "generalMatrix"))@x,
			      0, ..., na.rm = na.rm)
	      }
	  }
	  else { ## fully non-zero - very rare
	      callGeneric(x@x, ..., na.rm = na.rm)
	  }
      })


## "Ops" ("Arith", "Compare", "Logic") --> ./Ops.R

## -- end{group generics} -----------------------




## Methods for single-argument transformations

setMethod("zapsmall", signature = list(x = "dMatrix"),
          function(x, digits = getOption("digits")) {
              x@x <- zapsmall(x@x, digits)
              x
          })

## -- end(single-argument transformations) ------
