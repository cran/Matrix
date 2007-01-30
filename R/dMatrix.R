### Define Methods that can be inherited for all subclasses


setAs("dMatrix", "matrix",
      function(from) as(as(from, "dgeMatrix"), "matrix"))

##-> "dMatrix" <--> "lMatrix"   ---> ./lMatrix.R

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
## >>> More specific methods for sub-classes (sparse), use these as "catch-all":

## the non-Ops ones :
setMethod("Math2",
          ## Assume that  Generic(u, k) |--> u for u in {0,1}
          ## which is true for round(), signif() ==> all structure maintained
          signature(x = "dMatrix", digits = "numeric"),
	  function(x, digits) {
              x@x <- callGeneric(x@x, digits = digits)
              x
          })

## round(x) == round(x, 0)  etc
setMethod("Math2",
	  signature(x = "dMatrix", digits = "missing"),
	  function(x, digits) callGeneric(x, digits = 0))

## This needs extra work in ./AllGeneric.R :
setMethod("Summary", signature(x = "dMatrix", na.rm = "ANY"),
          function(x, ..., na.rm) callGeneric(x@x, ..., na.rm = na.rm))

## "Ops" ("Arith", "Compare", "Logic") --> ./Ops.R

## -- end{group generics} -----------------------




## Methods for single-argument transformations

setMethod("zapsmall", signature = list(x = "dMatrix"),
          function(x, digits = getOption("digits")) {
              x@x <- zapsmall(x@x, digits)
              x
          })

## -- end(single-argument transformations) ------
