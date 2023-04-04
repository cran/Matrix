## Methods for virtual class "symmetricMatrix" of symmetric matrices

## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("Matrix", "symmetricMatrix", ..M2sym)
setAs("matrix", "symmetricMatrix", ..M2sym)

## MJ: prefer more general methods above going via forceSymmetric(),
##     for which specialized methods should be defined
## MJ: methods above _do_ tolerate numerical fuzz; those below _did not_ ...
if(FALSE) {
setAs("denseMatrix", "symmetricMatrix", # checking symmetry vvvv
      function(from) .Call(dense_to_symmetric, from, "U",   TRUE))
setAs("matrix", "symmetricMatrix",      # checking symmetry vvvv
      function(from) .Call(dense_to_symmetric, from, "U",   TRUE))

## autogenerate coercions
##  as(*,  "symmetricMatrix")
##  ~~~~~~~~~~~~~~~~~~~~~~~~~
## Basic problem:
## This should work at package install time when package:Matrix does not exist!
if(FALSE)
local({
    allCl <- getClasses("package:Matrix") ## << fails at install time!!!!
    clss <- allCl[sapply(allCl, extends, class2 = "Matrix")]
    virt <- sapply(clss, isVirtualClass)
    ## Now ensure coercions for all  non-virtual "Matrix" inheriting classes:
    for(cl in clss[!virt]) {
        cld <- getClassDef(cl)
        if(extends(cld, "symmetricMatrix"))
            cat("\tsymmetric:\t", cl,"\n")
        else if(extends(cld, "triangularMatrix"))
            cat("\ttriangular:\t", cl,"\n")
        else if(extends(cld, "diagonalMatrix"))
            cat("\tdiagonal:\t", cl,"\n")
        else {
            cat("do ",cl,"\n")
##             setAs(cl, "symmetricMatrix",
##                   function(from) as(from, ".s.Matrix"))
        }
    }## for
})
} ## MJ


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dimnames", signature(x = "symmetricMatrix"),
          function(x) symmDN(x@Dimnames))

setMethod("isSymmetric", signature(object = "symmetricMatrix"),
          function(object, ...) TRUE)

setMethod("isTriangular", signature(object = "symmetricMatrix"),
          function(object, upper = NA, ...) {
              if(!isDiagonal(object))
                  FALSE
              else if(is.na(upper))
                  `attr<-`(TRUE, "kind", "U")
              else TRUE
          })

## MJ: no longer needed ...
##     methods for [CRT]sparseMatrix, (un)?packedMatrix do the same in C
##     while handling Hermitian zMatrix properly
if(FALSE) {
setMethod("symmpart", signature(x = "symmetricMatrix"),
          function(x) x)
setMethod("skewpart", signature(x = "symmetricMatrix"),
          function(x) .setZero(x))
} ## MJ

## MJ: no longer needed ...
##     methods for [CRT]sparseMatrix, (un)?packedMatrix do the same in C
##     and it is best to avoid hacks like these where possible
if(FALSE) {
.sM.force1 <- function(x, uplo) x
.sM.force2 <- function(x, uplo) if(uplo == x@uplo) x else t(x)
.sM.subclasses <- names(getClassDef("symmetricMatrix")@subclasses)
for (.cl in setdiff(.sM.subclasses, c("dpoMatrix", "dppMatrix", "corMatrix"))) {
    setMethod("forceSymmetric", signature(x = .cl, uplo = "missing"),
	      .sM.force1)
    setMethod("forceSymmetric", signature(x = .cl, uplo = "character"),
	      .sM.force2)
}
rm(.sM.force1, .sM.force2, .cl)
} ## MJ

## MJ: no longer needed ... replacement in ./(un)?packedMatrix.R
if(FALSE) {
## forceSymmetric() coerces to "symmetricMatrix"  withOUT  testing
## ---------------- contrary to  as(M, <symmetric>)  which should only
## work when 'M' is a symmetric matrix __ in the sense of isSymmetric() __
## i.e., with allowing a little bit of asymmetric numeric fuzz:
setMethod("forceSymmetric", signature(x = "denseMatrix", uplo = "missing"),
	  function(x, uplo) {
	      uplo <- if(is(x, "triangularMatrix")) x@uplo else "U"
	      .Call(dense_to_symmetric, x, uplo, FALSE) ## FIXME? diagU2N()?
          })
setMethod("forceSymmetric", signature(x = "denseMatrix", uplo = "character"),
	  function(x, uplo) .Call(dense_to_symmetric, x, uplo, FALSE))
setMethod("forceSymmetric", signature(x = "matrix", uplo = "missing"),
	  function(x, uplo) .Call(dense_to_symmetric, x,  "U", FALSE))
setMethod("forceSymmetric", signature(x = "matrix", uplo = "character"),
	  function(x, uplo) .Call(dense_to_symmetric, x, uplo, FALSE))

###------- pack() and unpack() --- for *dense*  symmetric & triangular matrices:
packM <- function(x, Mtype, kind, unpack=FALSE) {
    cd <- getClassDef(cx <- class(x))
    if(extends(cd, "sparseMatrix"))
	stop(sprintf("(un)packing only applies to dense matrices, class(x)='%s'",
		     cx))
    if(!missing(kind) && kind == "symmetric") { ## use 'unpack' but not 'Mtype'
	## special treatment for positive definite ones:
	as(x, if(unpack) {
	    if(extends(cd, "dppMatrix")) "dpoMatrix"
	    else paste0(.M.kindC(cd), "syMatrix")
	} else { ## !unpack : "pack" :
	    if(extends(cd, "dpoMatrix")) "dppMatrix"
	    else paste0(.M.kindC(cd), "spMatrix")
	})
    } else as(x, paste0(.M.kindC(cd), Mtype))
}
setMethod("unpack", "symmetricMatrix",
          function(x, ...) packM(x, kind="symmetric", unpack=TRUE))
setMethod("pack",   "symmetricMatrix", function(x, ...) packM(x, kind="symmetric"))
setMethod("unpack", "triangularMatrix",
	  function(x, ...) packM(x, "trMatrix", unpack=TRUE))
setMethod("pack",   "triangularMatrix", function(x, ...) packM(x, "tpMatrix"))
## to produce a nicer error message:
pckErr <- function(x, ...)
    stop(sprintf("(un)packing only applies to dense matrices, class(x)='%s'",
		 class(x)))
setMethod("unpack", "sparseMatrix", pckErr)
setMethod("pack",   "sparseMatrix", pckErr)
rm(pckErr)

##' pack (<matrix>)  -- smartly:
setMethod("pack", signature(x = "matrix"),
	  function(x, symmetric=NA, upperTri = NA, ...) {
	      if(is.na(symmetric)) ## must guess symmetric / triangular
		  symmetric <- isSymmetric.matrix(x)
	      if(symmetric) {
		  pack(.Call(dense_to_symmetric, x, "U", TRUE), ...)
	      } else { # triangular
		  ## isTriMat(..) : should still check fully (speed up ?) ..
		  if(isTr <- isTriMat(x, upper=upperTri)) {
		      uplo <- attr(isTr, "kind")
		      pack(new(paste0(.M.kind(x),"tpMatrix"),
			       x = x[indTri(nrow(x), upper=(uplo == "U"), diag=TRUE)],
			       Dim = dim(x), Dimnames = .M.DN(x), uplo = uplo), ...)
		  } else
		      stop("'x' is not symmetric nor triangular")
	      }
	  })

## {"traditional"} specific methods
setMethod("unpack", "dspMatrix",
	  function(x, ...) dsp2dsy(x), valueClass = "dsyMatrix")
setMethod("unpack", "dtpMatrix",
	  function(x, ...) dtp2dtr(x), valueClass = "dtrMatrix")
} ## MJ
