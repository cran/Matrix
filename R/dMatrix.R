### Define Methods that can be inherited for all subclasses

if(paste(R.version$major, R.version$minor, sep=".") < "2.0.1") {
    ## for R 2.0.0 (and earlier)
    setAs("dMatrix", "dgeMatrix",
	  function(from) new("dgeMatrix", x = from@x,
			     Dim = from@Dim, Dimnames = from@Dimnames))
}

setAs("dMatrix", "matrix",
      function(from) as(as(from, "dgeMatrix"), "matrix"))

## Methods for operations where one argument is integer

setMethod("%*%", signature(x = "dMatrix", y = "integer"),
          function(x, y)
          callGeneric(x, array(as.numeric(y), c(length(y), 1))),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "integer", y = "dMatrix"),
          function(x, y)
          callGeneric(array(as.numeric(x), c(1, length(x))), y),
          valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dMatrix", y = "integer"),
          function(x, y = NULL)
          callGeneric(x, array(as.numeric(y), c(length(y), 1))),
          valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "integer", y = "dMatrix"),
          function(x, y = NULL)
          callGeneric(array(as.numeric(x), c(1, length(x))), y),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dMatrix", b = "integer"),
          function(a, b, ...)
          callGeneric(a, array(as.numeric(b), c(length(b), 1))),
          valueClass = "dgeMatrix")


## Group Methods, see ?Arith (e.g.)
## -----

## Cheap version: work via "dgeMatrix" and use the group methods there:

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
          signature(x = "dMatrix", digits = "numeric"),
          function(x, digits) callGeneric(as(x, "dgeMatrix"), digits = digits))

## TODO :  "Compare" -> returning  logical Matrices

## -- end{group generics} -----------------------
