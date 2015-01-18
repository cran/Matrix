### Define Methods that can be inherited for all subclasses

##-> "dMatrix" <--> "lMatrix"   ---> ./lMatrix.R

## these two are parallel to "n <-> l" in the above :
setAs("nMatrix", "dMatrix",
      function(from) {
	  cld <- getClassDef(cl <- MatrixClass(class(from)))
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
	  cld <- getClassDef(cl <- MatrixClass(class(from)))
	  if(extends(cld, "diagonalMatrix")) { # have no "ndi*" etc class
	      cl <- class(from <- as(from, "sparseMatrix"))
	      isSp <- TRUE
	  } else {
	      isSp <- extends(cld, "sparseMatrix")
	      if(isSp && any(from@x == 0)) {
		  from <- drop0(from) # was drop0(from, cld)
		  if(cl != (c. <- class(from)))
		      cld <- getClassDef(cl <- c.)
	      }
	  }
	  sNams <- slotNames(cld)
	  r <- copyClass(from, sub("^d", "n", cl), sNams[sNams != "x"])
	  if(!isSp) #  'x' slot |--> logical
	      r@x <- as.logical(from@x)
	  r
      })


## Group Methods:
## -----
## "Math", "Math2" in			--> ./Math.R
## "Summary"				--> ./Summary.R
## "Ops" ("Arith", "Compare", "Logic")	--> ./Ops.R



## Methods for single-argument transformations

setMethod("zapsmall", signature(x = "dMatrix"),
          function(x, digits = getOption("digits")) {
              x@x <- zapsmall(x@x, digits)
              x
          })

## -- end(single-argument transformations) ------


