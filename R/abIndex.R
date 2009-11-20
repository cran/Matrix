#### Methods for the "abIndex" := ``abstract Index'' class

### Note: this partly builds on ideas and code from  Jens Oehlschlaegel,
### ----  as implemented (in the GPL'ed part of) package 'ff'.

## Basic idea:  a vector  x  of integer indices often has long stretches
##              i, i+1, i+2, ...  such that diff(x) has stretches of '1'.
## 		Now keep x[1] =: first and diff(x) =: d,
## 		and use rle() to encode d.  Here, use a C version for rle()
rleMaybe <- function(i) .Call(Matrix_int_rle, i)

setAs("numeric", "abIndex",
      function(from) {
	  ans <- new("abIndex")
	  r <- ## NB: diff(.) is coerced to integer on C level
	      rleMaybe(diff(from))
	  if(is.null(r)) { # no "compression"
	      ans@kind <- if(is.integer(from)) "int32" else "double"
	      ans@x <- from
	  } else {
	      ans@kind <- "rleDiff"
	      ## ans@x <- integer(0) # <- prototype does that
	      ans@rleD <- new("rleDiff", first = from[1], rle = r)
	  }
	  ans
      })

setMethod("length", "abIndex", function(x)
	  if(identical(x@kind, "rleDiff"))
	  sum(x@rleD@rle$lengths)+ 1L else length(x@x))

abI2num <- function(from) {
    switch(from@kind,
	   "rleDiff" = {
	       x <- from@rleD
	       ## as  inverse.rle():
	       cumsum(c(x@first, rep.int(x@rle$values, x@rle$lengths)))
	   },
	   "int32" =, "double" = from@x)
}
setAs("abIndex", "numeric", abI2num)
setAs("abIndex", "integer", function(from) as.integer(abI2num(from)))
## for faint hearted useRs:
setMethod(as.numeric, "abIndex", function(x) abI2num(x))

## Need   max(<i>), min(<i>),   all(<i> == <j>)   any(<i> == <j>)

## --->  Groups  "Summary"  and "Compare"  (maybe all "Ops")

## For that, we really  need  "[" and/or  "rep"() methods -- TODO --
##
setMethod("[", signature(x = "abIndex", i = "index"),
	  function (x, i, j, ..., drop)
{
    switch(x@kind,
           "rleDiff" = {
               ## intIv() in ./sparseVector.R -- not memory-efficient (??)
	       n <- length(x)
               ii <- intIv(i, n)## ii : 1-based integer indices
               d <- x@rleD
               ## Now work with the equivalent of
               ##   cumsum(c(d@first, rep.int(d@rle$values, d@rle$lengths)))

               stop("<abIndex>[i]  is not yet implemented")
           },
	   "int32" =, "double" =
	   ## as it's not rle-packed, can remain simple:
	   x@x[i])
})

## Endpoints of all linear stretches -- auxiliary for range(.)
##' @param x an "rleDiff" object

##' @return numeric vector of end points of all linear stretches of x.
ends.rleD <- function(x)
{
    rl <- x@rle
    stopifnot(length(lens <- rl$lengths) == (m <- length(vals <- rl$values)))
    cumsum(c(x@first, lens*vals))
}

## Summary: { max, min, range, prod, sum, any, all } :
## have  'summGener1' := those without prod, sum

setMethod("Summary", signature(x = "abIndex", na.rm = "ANY"),
          function(x, ..., na.rm)
{
    switch(x@kind,
           "rleDiff" =
       {
           d <- x@rleD
           if(.Generic %in% c("range","min","max")) {
               callGeneric(ends.rleD(d), ..., na.rm=na.rm)
           } else { ## "sum", "prod" :
               switch(.Generic,
                      "all" = {
                          ## these often, but *not* always come in pairs
                          en <- ends.rleD(d)
                          ## so maybe it does not really help!

                          stop("all(<abIndex>) is not yet implemented")

                          ## all(c(d@first, d@rle$values), ..., na.rm=na.rm)
                      },
                      "any" = any(c(d@first, d@rle$values), ..., na.rm=na.rm),
                      "sum" = {
                          stop("sum(<abIndex>) is not yet implemented")
                      },
                      "prod"= {
                          stop("prod(<abIndex>) is not yet implemented")
                      })
           }
       },
           "int32" =, "double" = callGeneric(x@x, ..., na.rm = na.rm)
           )
})

setMethod("Ops", signature(e1 = "abIndex", e2 = "abIndex"),
	  function(e1, e2)
{
    l1 <- length(e1)
    l2 <- length(e2)
    mM <- range(l1,l2)

 stop("not yet implemented")

    if(mM[1] != mM[2]) { ## lengths differ
        if(mM[1] %% mM[2] != 0) ## identical warning as in main/arithmetic.c
            warning("longer object length\n\tis not a multiple of shorter object length")
        if(l1 < l2) {

        } else { ## l1 > l2

        }
    }
    switch(e1@kind,
	   "rleDiff" = {

	   },
	   "int32" =, "double" = {

           })

})

## FIXME:  1 + <abIndex>  is trivial: just modify 'first' !
setMethod("Ops", signature(e1 = "abIndex", e2 = "numeric"), ## "numeric" or "ANY"
	  function(e1, e2) callGeneric(e1, as(e2, "abIndex")))

setMethod("Ops", signature(e1 = "numeric", e2 = "abIndex"),
	  function(e1, e2) callGeneric(as(e1, "abIndex"), e2))





## Then I want something like  get.ind.sel(.)  [ ./Tsparse.R ] working,
## i.e. possibly   match(i, <abI>, nomatch = 0)
