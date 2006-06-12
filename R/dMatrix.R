### Define Methods that can be inherited for all subclasses


setAs("dMatrix", "matrix",
      function(from) as(as(from, "dgeMatrix"), "matrix"))

### Idea: Coercion between *VIRTUAL* classes -- as() chooses "closest" classes
### ----  should also work e.g. for  triangular --> triangular !
##
## setAs("dMatrix", "lMatrix",
##       function(from) {
##       })

setAs("dMatrix", "sparseMatrix", function(from) as(from, "dgCMatrix"))

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

setMethod("expm", signature(x = "dMatrix"),
          function(x) callGeneric(as(x, "dgeMatrix")))


## Group Methods, see ?Arith (e.g.)
## -----

## Cheap version: work via "dgeMatrix" and use the group methods there:
## FIXME: To improve by doing some of these for dsparse* !
setMethod("Arith", ##  "+", "-", "*", "^", "%%", "%/%", "/"
          signature(e1 = "dMatrix", e2 = "dMatrix"),
          function(e1, e2) callGeneric(as(e1, "dgeMatrix"),
                                       as(e2, "dgeMatrix")))
setMethod("Arith",
          signature(e1 = "dMatrix", e2 = "numeric"),
          function(e1, e2) callGeneric(as(e1, "dgeMatrix"), e2))
setMethod("Arith",
          signature(e1 = "numeric", e2 = "dMatrix"),
          function(e1, e2) callGeneric(e1, as(e2, "dgeMatrix")))

setMethod("Math",
          signature(x = "dMatrix"),
          function(x) callGeneric(as(x, "dgeMatrix")))

setMethod("Math2",
          ## Assume that  Generic(u, k) |--> u for u in {0,1}
          ## which is true for round(), signif() ==> all structure maintained
          signature(x = "dMatrix", digits = "numeric"),
	  function(x, digits) {
              x@x <- callGeneric(x@x, digits = digits)
              x
          })

## This needs extra work in ./AllGeneric.R :
setMethod("Summary", signature(x = "dMatrix", na.rm = "ANY"),
          function(x, ..., na.rm) callGeneric(x@x, ..., na.rm = na.rm))

## "Compare" -> returning  logical Matrices
setMethod("Compare", signature(e1 = "numeric", e2 = "dMatrix"),
          function(e1,e2) {
              switch(.Generic,
                     "==" =, "!=" = callGeneric(e2, e1),
                     "<" =, ">" =, "<=" =, ">=" = !callGeneric(e2, e1))
          })

setMethod("Compare", signature(e1 = "dMatrix", e2 = "numeric"),
	  function(e1, e2) {
              lClass <- dClass2(class(e1), "l")
              fullCl <- if(isSymmetric(e1)) "lsyMatrix" else "lgeMatrix"
	      ## Dbg cat("Compare", class(e1), "|-> ",lClass, "\n")
	      r  <- callGeneric(e1@x, e2)
              r0 <- callGeneric(0, e2)
	      if(is(e1, "denseMatrix")) {
                  full <- !isPacked(e1) # << both "dtr" and "dsy" are 'full'
                  if(full || identical(r0, FALSE) || is(e1, "symmetricMatrix"))
                      r <- new(lClass, x = r,
                               Dim = dim(e1), Dimnames = dimnames(e1))
                  else { ## packed matrix with structural 0 and r0 is not FALSE:
                      ##--> result cannot be packed anymore
                      dr <- as(r, fullCl)
                      ## FIXME: implement this:
                      dr[ind.0(e1)] <- r0
                  }
	      }
              else { ## dsparseMatrix
                  r <- new(lClass, Dim = dim(e1), Dimnames = dimnames(e1),
                           x = r)
                  for(n in setdiff(c("i","j","p"), slotNames(r)))
                      slot(r, n) <- slot(e1, n)

		  if(!identical(r0, FALSE)) {
                      warning("sparse to dense coercion in ",.Generic)

                      .bail.out.2(.Generic, class(e1), class(e2))
                      ## NOT YET:

                      ## non sparse result
                      dr <- as(r,
                               if(isSymmetric(r))"lsyMatrix" else "lgeMatrix")
                      ## FIXME: implement this:
                      dr[ind.0(e1)] <- r0
		  }
	      }
              r
	  })

## "dMatrix <-> work with 'x' slot
setMethod("Compare", signature(e1 = "dMatrix", e2 = "dMatrix"),
          function(e1, e2) {
              d <- dimCheck(e1,e2)
              lClass <- dClass2(class(e1), "l")

              ## FIXME: if (the 'x' are slots compatible)
	      r <- callGeneric(e1@x, e2@x)
	      if(is(e1, "denseMatrix")) {
		  r <- new(lClass, x = r,
                           Dim = dim(e1), Dimnames = dimnames(e1))
	      }
              else { ## dsparseMatrix

		  if(identical(FALSE, r0 <- callGeneric(0, e2))) {
		      ## return (potentially even more) sparse logical Matrix
		      r <- new(lClass, x = r,
                               Dim = dim(e1), Dimnames = dimnames(e1))

		  } else { ## non sparse result

		  stop("'Compare' for sparse dMatrix not yet implemented for all cases")
### FIXME
		  }
	      }
              r
	  })

## -- end{group generics} -----------------------




## Methods for single-argument transformations

setMethod("zapsmall", signature = list(x = "dMatrix"),
          function(x, digits = getOption("digits")) {
              x@x <- zapsmall(x@x, digits)
              x
          })

## -- end(single-argument transformations) ------
