### Define Methods that can be inherited for all subclasses

if(paste(R.version$major, R.version$minor, sep=".") < "2.0.1") {
    ## for R 2.0.0 (and earlier)
    setAs("dMatrix", "dgeMatrix",
	  function(from) new("dgeMatrix", x = from@x,
			     Dim = from@Dim, Dimnames = from@Dimnames))
}

setAs("dMatrix", "matrix",
      function(from) as(as(from, "dgeMatrix"), "matrix"))

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
