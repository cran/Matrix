#### All  %*%, crossprod() and tcrossprod() methods of the Matrix package
#### ^^^  ----------------------------------------------------------
###  with EXCEPTIONS:	./diagMatrix.R
###			./indMatrix.R ./pMatrix.R

### NOTA BENE:   vector %*% Matrix  _and_  Matrix %*% vector
### ---------   The k-vector is treated as  (1,k)-matrix *or* (k,1)-matrix
### on both sides when ever it "helps fit" the matrix dimensions:
##--- ./products.Rout
##    ~~~~~~~~~~~~~~~
## ========> in a M.v or v.M operation ,
##           you *must* look at dim(M) to see how to treat  v  !!!!!!!!!!!!!!!!

## For %*% (M = Matrix; v = vector (double, integer,.. or "sparsevector"):
## Drawback / bug: for (dense)vectors, the *names* are lost [sparsevectors have no names!]
.M.v <- function(x, y) { #
    dim(y) <- if(ncol(x) == (n <- length(y)))
        c(n, 1L) else c(1L, n) ## which works when m == 1, otherwise errors
    x %*% y
}

## For %*% :
.v.M <- function(x, y) {
    dim(x) <- if(nrow(y) == (n <- length(x))) c(1L, n) else c(n, 1L)
    x %*% y
}

## For tcrossprod() :
.v.Mt <- function(x, y) {
    dim(x) <- if(ncol(y) == (n <- length(x))) c(1L, n) else c(n, 1L)
    tcrossprod(x, y)
}

###-- I --- %*% ------------------------------------------------------

## General method for dense matrix multiplication in case specific methods
## have not been defined.
setMethod("%*%", signature(x = "ddenseMatrix", y = "ddenseMatrix"),
	  function(x, y) .Call(dgeMatrix_matrix_mm,
			       .Call(dup_mMatrix_as_dgeMatrix, x), y, FALSE),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y) .Call(dgeMatrix_matrix_mm, x, y, FALSE),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y) .Call(dgeMatrix_matrix_mm, x, y, FALSE),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y) .Call(dgeMatrix_matrix_mm, y, x, TRUE),
	  valueClass = "dgeMatrix")

.dsy_m_mm <- function(x, y) .Call(dsyMatrix_matrix_mm, x, y, FALSE)
setMethod("%*%", signature(x = "dsyMatrix", y = "matrix"),  .dsy_m_mm)
setMethod("%*%", signature(x = "dsyMatrix", y = "ddenseMatrix"),  .dsy_m_mm)
## for disambiguity :
setMethod("%*%", signature(x = "dsyMatrix", y = "dsyMatrix"),  .dsy_m_mm)
## or even
## for(yCl in .directSubClasses(getClass("ddenseMatrix")))
##     setMethod("%*%", signature(x = "dsyMatrix", y = yCl), .dsy_m_mm)

setMethod("%*%", signature(x = "ddenseMatrix", y = "dsyMatrix"),
          function(x, y) .Call(dsyMatrix_matrix_mm, y, x, TRUE))
setMethod("%*%", signature(x = "matrix", y = "dsyMatrix"),
          function(x, y) .Call(dsyMatrix_matrix_mm, y, x, TRUE))

setMethod("%*%", signature(x = "dspMatrix", y = "ddenseMatrix"),
          function(x, y) .Call(dspMatrix_matrix_mm, x, y),
          valueClass = "dgeMatrix")
setMethod("%*%", signature(x = "dspMatrix", y = "matrix"),
          function(x, y) .Call(dspMatrix_matrix_mm, x, y),
          valueClass = "dgeMatrix")


## Not needed because of c("numeric", "Matrix") method
##setMethod("%*%", signature(x = "numeric", y = "CsparseMatrix"),
##	    function(x, y) t(.Call(Csparse_dense_crossprod, y, t(x))),
##	    valueClass = "dgeMatrix")

## FIXME -- do the "same" for "dtpMatrix" {also, with [t]crossprod()}
## all just like these "%*%" :
setMethod("%*%", signature(x = "dtrMatrix", y = "dtrMatrix"),
	  function(x, y) .Call(dtrMatrix_dtrMatrix_mm, x, y, FALSE, FALSE))

setMethod("%*%", signature(x = "dtrMatrix", y = "ddenseMatrix"),
	  function(x, y) .Call(dtrMatrix_matrix_mm, x, y, FALSE, FALSE),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dtrMatrix", y = "matrix"),
	  function(x, y) .Call(dtrMatrix_matrix_mm, x, y, FALSE, FALSE),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "ddenseMatrix", y = "dtrMatrix"),
	  function(x, y) .Call(dtrMatrix_matrix_mm, y, x, TRUE, FALSE),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "matrix", y = "dtrMatrix"),
	  function(x, y) .Call(dtrMatrix_matrix_mm, y, x, TRUE, FALSE),
	  valueClass = "dgeMatrix")



setMethod("%*%", signature(x = "dtpMatrix", y = "ddenseMatrix"),
	  function(x, y) .Call(dtpMatrix_matrix_mm, x, y, FALSE, FALSE))
setMethod("%*%", signature(x = "dgeMatrix", y = "dtpMatrix"),
	  function(x, y) .Call(dgeMatrix_dtpMatrix_mm, x, y))

## dtpMatrix <-> matrix : will be used by the "numeric" one
setMethod("%*%", signature(x = "dtpMatrix", y = "matrix"),
          function(x, y) .Call(dtpMatrix_matrix_mm, x, y, FALSE, FALSE))
setMethod("%*%", signature(x = "matrix", y = "dtpMatrix"),
          function(x, y) ..2dge(x) %*% y)

## dtpMatrix <-> numeric : the auxiliary functions are R version specific!
##setMethod("%*%", signature(x = "dtpMatrix", y = "numeric"), .M.v)
##setMethod("%*%", signature(x = "numeric", y = "dtpMatrix"), .v.M)


## For multiplication operations, sparseMatrix overrides other method
## selections.	Coerce a ddensematrix argument to a lsparseMatrix.
setMethod("%*%", signature(x = "lsparseMatrix", y = "ldenseMatrix"),
	  function(x, y) x %*% as(y, "sparseMatrix"))

setMethod("%*%", signature(x = "ldenseMatrix", y = "lsparseMatrix"),
	  function(x, y) as(x, "sparseMatrix") %*% y)

## and coerce lsparse* to lgC*
setMethod("%*%", signature(x = "lsparseMatrix", y = "lsparseMatrix"),
	  function(x, y) as(x, "lgCMatrix") %*% as(y, "lgCMatrix"))


for(c.x in c("lMatrix", "nMatrix")) {
    setMethod("%*%", signature(x = c.x, y = "dMatrix"),
	      function(x, y) as(x, "dMatrix") %*% y)
    setMethod("%*%", signature(x = "dMatrix", y = c.x),
	      function(x, y) x %*% as(y, "dMatrix"))
    for(c.y in c("lMatrix", "nMatrix"))
    setMethod("%*%", signature(x = c.x, y = c.y),
	      function(x, y) as(x, "dMatrix") %*% as(y, "dMatrix"))
}; rm(c.x, c.y)

setMethod("%*%", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
	  function(x, y) .Call(Csparse_Csparse_prod, x, y))

setMethod("%*%", signature(x = "CsparseMatrix", y = "ddenseMatrix"),
	  function(x, y) .Call(Csparse_dense_prod, x, y))
setMethod("%*%", signature(x = "CsparseMatrix", y = "matrix"),
	  function(x, y) .Call(Csparse_dense_prod, x, y)) # was  x %*% Matrix(y)
setMethod("%*%", signature(x = "CsparseMatrix", y = "numeric"),
	  function(x, y) .Call(Csparse_dense_prod, x, y))


## Not yet.  Don't have methods for y = "CsparseMatrix" and general x
#setMethod("%*%", signature(x = "ANY", y = "TsparseMatrix"),
#	   function(x, y) callGeneric(x, as(y, "CsparseMatrix")))

setMethod("%*%", signature(x = "TsparseMatrix", y = "ANY"),
	  function(x, y) .T.2.C(x) %*% y)
setMethod("%*%", signature(x = "ANY", y = "TsparseMatrix"),
	  function(x, y) x %*% .T.2.C(y))
setMethod("%*%", signature(x = "TsparseMatrix", y = "Matrix"),
	  function(x, y) .T.2.C(x) %*% y)
setMethod("%*%", signature(x = "Matrix", y = "TsparseMatrix"),
	  function(x, y) x %*% .T.2.C(y))
setMethod("%*%", signature(x = "TsparseMatrix", y = "TsparseMatrix"),
	  function(x, y) .T.2.C(x) %*% .T.2.C(y))



##-------- Work via  as(*, lgC) : ------------

## For multiplication operations, sparseMatrix overrides other method
## selections.	Coerce a ddensematrix argument to a nsparseMatrix.
setMethod("%*%", signature(x = "nsparseMatrix", y = "ndenseMatrix"),
	  function(x, y) x %*% as(y, "nsparseMatrix"))

setMethod("%*%", signature(x = "ndenseMatrix", y = "nsparseMatrix"),
	  function(x, y) as(x, "nsparseMatrix") %*% y)
## and coerce nsparse* to lgC*
setMethod("%*%", signature(x = "nsparseMatrix", y = "nsparseMatrix"),
	  function(x, y) as(x, "ngCMatrix") %*% as(y, "ngCMatrix"))


## FIXME(2): These three are sub-optimal : has  2 x  t(<dense>)  :
##   *faster*: provide   dense_Csparse_prod()
## x %*% y =  t(crossprod(y, t(x)))  unless when x is vector
setMethod("%*%", signature(x = "ddenseMatrix", y = "CsparseMatrix"),
	  function(x, y) t(.Call(Csparse_dense_crossprod, y, t(x))),
	  valueClass = "dgeMatrix")
setMethod("%*%", signature(x = "matrix", y = "CsparseMatrix"),
	  function(x, y) t(.Call(Csparse_dense_crossprod, y, t(x))),
	  valueClass = "dgeMatrix")
setMethod("%*%", signature(x = "numLike", y = "CsparseMatrix"),
	  function(x, y) t(.Call(Csparse_dense_crossprod, y, x)),
	  valueClass = "dgeMatrix")

## "Matrix"
## Methods for operations where one argument is numeric
setMethod("%*%", signature(x = "Matrix", y = "numLike"), .M.v)
setMethod("%*%", signature(x = "numLike", y = "Matrix"), .v.M)

setMethod("%*%", signature(x = "Matrix", y = "matrix"),
	  function(x, y) x %*% Matrix(y))
setMethod("%*%", signature(x = "matrix", y = "Matrix"),
	  function(x, y) Matrix(x) %*% y)

## bail-out methods in order to get better error messages
.local.bail.out <- function (x, y)
    stop(gettextf('not-yet-implemented method for <%s> %%*%% <%s>',
		  class(x), class(y)), domain=NA)
setMethod("%*%", signature(x = "ANY", y = "Matrix"), .local.bail.out)
setMethod("%*%", signature(x = "Matrix", y = "ANY"), .local.bail.out)


### sparseVector
sp.x.sp <- function(x, y) Matrix(sum(x * y), 1L, 1L, sparse=FALSE)
    ## inner product -- no sense to return sparse!
setMethod("%*%", signature(x = "Matrix", y = "sparseVector"), .M.v)
setMethod("%*%", signature(x = "sparseVector", y = "Matrix"), .v.M)
setMethod("%*%", signature(x = "sparseVector", y = "sparseVector"), sp.x.sp)
## setMethod("%*%", signature(x = "sparseMatrix", y = "sparseVector"),
##           function(x, y) x %*% .sparseV2Mat(y))

###--- II --- crossprod -----------------------------------------------------

setMethod("crossprod", signature(x = "dgeMatrix", y = "missing"),
	  function(x, y = NULL) .Call(dgeMatrix_crossprod, x, FALSE),
	  valueClass = "dpoMatrix")

## crossprod (x,y)
setMethod("crossprod", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y = NULL) .Call(dgeMatrix_dgeMatrix_crossprod, x, y, FALSE),
	  valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y = NULL) .Call(dgeMatrix_matrix_crossprod, x, y, FALSE),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "dgeMatrix", y = "numeric"),
	  function(x, y = NULL) .Call(dgeMatrix_matrix_crossprod, x, y, FALSE),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y = NULL) crossprod(..2dge(x), y),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "numLike", y = "dgeMatrix"),
	  function(x, y = NULL) crossprod(as.matrix(as.double(x)), y),
	  valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "ddenseMatrix", y = "missing"),
	  function(x, y = NULL) crossprod(as(x, "dgeMatrix")))

setMethod("crossprod", signature(x = "dtrMatrix", y = "missing"),
	  function(x, y = NULL) crossprod(as(x, "dgeMatrix")),
	  valueClass = "dpoMatrix")

## "dtrMatrix" - remaining (uni)triangular if possible
setMethod("crossprod", signature(x = "dtrMatrix", y = "dtrMatrix"),
	  function(x, y) .Call(dtrMatrix_dtrMatrix_mm, x, y, FALSE, TRUE))

setMethod("crossprod", signature(x = "dtrMatrix", y = "ddenseMatrix"),
	  function(x, y) .Call(dtrMatrix_matrix_mm, x, y, FALSE, TRUE),
	  valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dtrMatrix", y = "matrix"),
	  function(x, y) .Call(dtrMatrix_matrix_mm, x, y, FALSE, TRUE),
	  valueClass = "dgeMatrix")
## "dtpMatrix"
if(FALSE) ## not yet in C
setMethod("crossprod", signature(x = "dtpMatrix", y = "dtpMatrix"),
	  function(x, y) .Call(dtpMatrix_dtpMatrix_mm, x, y, FALSE, TRUE))

setMethod("crossprod", signature(x = "dtpMatrix", y = "ddenseMatrix"),
	  function(x, y) .Call(dtpMatrix_matrix_mm, x, y, FALSE, TRUE),
	  valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dtpMatrix", y = "matrix"),
	  function(x, y) .Call(dtpMatrix_matrix_mm, x, y, FALSE, TRUE),
	  valueClass = "dgeMatrix")



## "crossprod" methods too ...
## setMethod("crossprod", signature(x = "dgTMatrix", y = "missing"),
##	     function(x, y = NULL)
##	     .Call(csc_crossprod, as(x, "dgCMatrix")))

## setMethod("crossprod", signature(x = "dgTMatrix", y = "matrix"),
##	     function(x, y = NULL)
##	     .Call(csc_matrix_crossprod, as(x, "dgCMatrix"), y))

##setMethod("crossprod", signature(x = "dgTMatrix", y = "numeric"),
##	    function(x, y = NULL)
##	    .Call(csc_matrix_crossprod, as(x, "dgCMatrix"), as.matrix(y)))

## setMethod("tcrossprod", signature(x = "dgTMatrix", y = "missing"),
##	     function(x, y = NULL)
##	     .Call(csc_tcrossprod, as(x, "dgCMatrix")))

setMethod("crossprod", signature(x = "CsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
	      if (is(x, "symmetricMatrix"))
		  ## crossprod() should give "symmetric*":
		  forceSymmetric(x %*% x, uplo = x@uplo)
	      else
		  .Call(Csparse_crossprod, x, trans = FALSE, triplet = FALSE)
	  })

setMethod("crossprod", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
	  function(x, y = NULL)
	  .Call(Csparse_Csparse_crossprod, x, y, trans = FALSE))

## FIXME: Generalize the class of y.  This specific method is to replace one
##	  in dgCMatrix.R
setMethod("crossprod", signature(x = "CsparseMatrix", y = "ddenseMatrix"),
	  function(x, y = NULL) .Call(Csparse_dense_crossprod, x, y))
setMethod("crossprod", signature(x = "CsparseMatrix", y = "matrix"),
	  function(x, y = NULL) .Call(Csparse_dense_crossprod, x, y))
setMethod("crossprod", signature(x = "CsparseMatrix", y = "numeric"),
	  function(x, y = NULL) .Call(Csparse_dense_crossprod, x, y))


setMethod("crossprod", signature(x = "TsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
	      if (is(x, "symmetricMatrix"))
		  ## crossprod() should give "symmetric*":
		  forceSymmetric(x %*% x, uplo = x@uplo)
	      else
		  .Call(Csparse_crossprod, x, trans = FALSE, triplet = TRUE)
	  })

setMethod("crossprod", signature(x = "TsparseMatrix", y = "ANY"),
	  function(x, y = NULL) crossprod(.T.2.C(x), y))
setMethod("crossprod", signature(x = "ANY", y = "TsparseMatrix"),
	  function(x, y = NULL) crossprod(x, .T.2.C(y)))
setMethod("crossprod", signature(x = "TsparseMatrix", y = "Matrix"),
	  function(x, y = NULL) crossprod(.T.2.C(x), y))
setMethod("crossprod", signature(x = "Matrix", y = "TsparseMatrix"),
	  function(x, y = NULL) crossprod(x, .T.2.C(y)))
setMethod("crossprod", signature(x = "TsparseMatrix", y = "TsparseMatrix"),
	  function(x, y = NULL) crossprod(.T.2.C(x), .T.2.C(y)))


setMethod("crossprod", signature(x = "dsparseMatrix", y = "ddenseMatrix"),
	  function(x, y = NULL)
	  .Call(Csparse_dense_crossprod, as(x, "CsparseMatrix"), y))

setMethod("crossprod", signature(x = "ddenseMatrix", y = "dgCMatrix"),
	  function(x, y = NULL) t(.Call(Csparse_dense_crossprod, y, x)))
setMethod("crossprod", signature(x = "ddenseMatrix", y = "dsparseMatrix"),
	  function(x, y = NULL)
	  t(.Call(Csparse_dense_crossprod, as(y, "CsparseMatrix"), x)))

setMethod("crossprod", signature(x = "dgCMatrix", y = "dgeMatrix"),
	  function(x, y = NULL) .Call(Csparse_dense_crossprod, x, y))
setMethod("crossprod", signature(x = "dsparseMatrix", y = "dgeMatrix"),
	  function(x, y = NULL)
	  .Call(Csparse_dense_crossprod, as(x, "CsparseMatrix"), y))

## NB: there's already
##     ("CsparseMatrix", "missing") and ("TsparseMatrix", "missing") methods

## infinite recursion:
## setMethod("crossprod", signature(x = "dgeMatrix", y = "dsparseMatrix"),
##	  function(x, y = NULL) crossprod(x, as(y, "dgCMatrix")))


setMethod("crossprod", signature(x = "lsparseMatrix", y = "ldenseMatrix"),
	  function(x, y = NULL) crossprod(x, as(y, "sparseMatrix")))

setMethod("crossprod", signature(x = "ldenseMatrix", y = "lsparseMatrix"),
	  function(x, y = NULL) crossprod(as(x, "sparseMatrix"), y))

setMethod("crossprod", signature(x = "lsparseMatrix", y = "lsparseMatrix"),
	  function(x, y = NULL)
	  crossprod(as(x, "lgCMatrix"), as(y, "lgCMatrix")))

setMethod("crossprod", signature(x = "nsparseMatrix", y = "ndenseMatrix"),
	  function(x, y = NULL) crossprod(x, as(y, "sparseMatrix")))

setMethod("crossprod", signature(x = "ndenseMatrix", y = "nsparseMatrix"),
	  function(x, y = NULL) crossprod(as(x, "sparseMatrix"), y))

setMethod("crossprod", signature(x = "nsparseMatrix", y = "nsparseMatrix"),
	  function(x, y = NULL)
	  crossprod(as(x, "ngCMatrix"), as(y, "ngCMatrix")))


## FIXME(3): slightly sub-optimal : t(<dense>)	:
setMethod("crossprod", signature(x = "ddenseMatrix", y = "CsparseMatrix"),
	  function(x, y) t(.Call(Csparse_dense_crossprod, y, x)))
setMethod("crossprod", signature(x = "matrix",	     y = "CsparseMatrix"),
	  function(x, y) t(.Call(Csparse_dense_crossprod, y, x)))
setMethod("crossprod", signature(x = "numeric",	     y = "CsparseMatrix"),
	  function(x, y) t(.Call(Csparse_dense_crossprod, y, x)))


## "Matrix" : cbind(), rbind() do  names -> dimnames
setMethod("crossprod", signature(x = "Matrix", y = "numLike"),
	  function(x, y) crossprod(x, cbind(y,deparse.level=0)))
setMethod("crossprod", signature(x = "numLike", y = "Matrix"), 
	  function(x, y) crossprod(rbind(x,deparse.level=0), y))

setMethod("crossprod", signature(x = "Matrix", y = "matrix"),
	  function(x, y) crossprod(x, Matrix(y)))
setMethod("crossprod", signature(x = "matrix", y = "Matrix"),
	  function(x, y) crossprod(Matrix(x), y))

## sparseVector
setMethod("crossprod", signature(x = "Matrix", y = "sparseVector"),
	  function(x, y) crossprod(x, .sparseV2Mat(y)))
setMethod("crossprod", signature(x = "sparseVector", y = "Matrix"), 
	  function(x, y)
	  crossprod(spV2M(x, nrow = length(x), ncol = 1L, check = FALSE), y))

setMethod("crossprod", signature(x = "sparseVector", y = "sparseVector"), sp.x.sp)
setMethod("crossprod", signature(x = "sparseVector", y = "missing"),
	  function(x, y=NULL) sp.x.sp(x,x))

## cheap fallbacks
setMethod("crossprod", signature(x = "Matrix", y = "Matrix"),
	  function(x, y) t(x) %*% y)
setMethod("crossprod", signature(x = "Matrix", y = "missing"),
	  function(x, y) t(x) %*% x)
setMethod("crossprod", signature(x = "Matrix", y = "ANY"),
	  function(x, y) t(x) %*% y)
setMethod("crossprod", signature(x = "ANY", y = "Matrix"),
	  function(x, y) t(x) %*% y)

###--- III --- tcrossprod ---------------------------------------------------

setMethod("tcrossprod", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y = NULL) .Call(dgeMatrix_dgeMatrix_crossprod, x, y, TRUE),
	  valueClass = "dgeMatrix")

setMethod("tcrossprod", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y = NULL) .Call(dgeMatrix_matrix_crossprod, x, y, TRUE),
	  valueClass = "dgeMatrix")
setMethod("tcrossprod", signature(x = "dgeMatrix", y = "numLike"),
	  function(x, y = NULL) .Call(dgeMatrix_matrix_crossprod, x, y, TRUE),
	  valueClass = "dgeMatrix")
setMethod("tcrossprod", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y = NULL) tcrossprod(..2dge(x), y),
	  valueClass = "dgeMatrix")
setMethod("tcrossprod", signature(x = "numLike", y = "dgeMatrix"), .v.Mt,
	  valueClass = "dgeMatrix")


setMethod("tcrossprod", signature(x = "dgeMatrix", y = "missing"),
	  function(x, y = NULL) .Call(dgeMatrix_crossprod, x, TRUE),
	  valueClass = "dpoMatrix")

if(FALSE) { ## this would mask 'base::tcrossprod'
setMethod("tcrossprod", signature(x = "matrix", y = "missing"),
	  function(x, y = NULL)
	  .Call(dgeMatrix_crossprod, ..2dge(x), TRUE),
	  valueClass = "dpoMatrix")

setMethod("tcrossprod", signature(x = "numLike", y = "missing"),
	  function(x, y = NULL) tcrossprod(as.matrix(as.double(x))))
}

setMethod("tcrossprod", signature(x = "ddenseMatrix", y = "missing"),
	  function(x, y = NULL) tcrossprod(as(x, "dgeMatrix")))


setMethod("tcrossprod", signature(x = "dtrMatrix", y = "dtrMatrix"),
	  function(x, y) .Call(dtrMatrix_dtrMatrix_mm, y, x, TRUE, TRUE))


## Must  have 1st arg. = "dtrMatrix" in  dtrMatrix_matrix_mm ():
## would need another way, to define  tcrossprod()  --- TODO? ---
##
## setMethod("tcrossprod", signature(x = "dtrMatrix", y = "ddenseMatrix"),
## 	  function(x, y) .Call(dtrMatrix_matrix_mm, y, x, TRUE, TRUE))

## setMethod("tcrossprod", signature(x = "dtrMatrix", y = "matrix"),
## 	  function(x, y) .Call(dtrMatrix_matrix_mm, y, x, TRUE, TRUE))

setMethod("tcrossprod", signature(x = "ddenseMatrix", y = "dtrMatrix"),
 	  function(x, y) .Call(dtrMatrix_matrix_mm, y, x, TRUE, TRUE))

setMethod("tcrossprod", signature(x = "matrix", y = "dtrMatrix"),
 	  function(x, y) .Call(dtrMatrix_matrix_mm, y, x, TRUE, TRUE))

if(FALSE) { ## TODO in C
setMethod("tcrossprod", signature(x = "ddenseMatrix", y = "dtpMatrix"),
 	  function(x, y) .Call(dtpMatrix_matrix_mm, y, x, TRUE, TRUE))

setMethod("tcrossprod", signature(x = "matrix", y = "dtpMatrix"),
 	  function(x, y) .Call(dtpMatrix_matrix_mm, y, x, TRUE, TRUE))
}




setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
	  function(x, y = NULL)
	  .Call(Csparse_Csparse_crossprod, x, y, trans = TRUE))

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
	      if (is(x, "symmetricMatrix"))
		  ## tcrossprod() should give "symmetric*":
		  forceSymmetric(x %*% x, uplo = x@uplo)
	      else
		  .Call(Csparse_crossprod, x, trans = TRUE, triplet = FALSE)
	  })

### FIXME (suboptimal):  one t(<dense>):
setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "ddenseMatrix"),
	  function(x, y) .Call(Csparse_dense_prod, x, t(y)))
setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "matrix"),
	  function(x, y) .Call(Csparse_dense_prod, x, t(y)))
setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "numLike"),
	  function(x, y) .Call(Csparse_dense_prod, x, rbind(y, deparse.level=0)))


### FIXME (speed):  Csparse_dense_crossprod  should also get a 'trans = TRUE'
## so we could have one less t(<dense>) in this:
## -- xy' = (yx')'
setMethod("tcrossprod", signature(x = "ddenseMatrix", y = "CsparseMatrix"),
	  function(x, y) t(.Call(Csparse_dense_prod, y, t(x))))
setMethod("tcrossprod", signature(x = "matrix",	      y = "CsparseMatrix"),
	  function(x, y) t(.Call(Csparse_dense_prod, y, t(x))))
setMethod("tcrossprod", signature(x = "numLike",      y = "CsparseMatrix"),
                                        # x or t(x) depending on dimension of y !
          .v.Mt)#<- FIXME more efficient



setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
	      if (is(x, "symmetricMatrix"))
		  ## tcrossprod() should give "symmetric*":
		  forceSymmetric(x %*% x, uplo = x@uplo)
	      else
		  .Call(Csparse_crossprod, x, trans = TRUE, triplet = TRUE)
	  })

setMethod("tcrossprod", signature(x = "ANY", y = "TsparseMatrix"),
	  function(x, y = NULL) tcrossprod(x, .T.2.C(y)))
setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "ANY"),
	  function(x, y = NULL) tcrossprod(.T.2.C(x), y))
setMethod("tcrossprod", signature(x = "Matrix", y = "TsparseMatrix"),
	  function(x, y = NULL) tcrossprod(x, .T.2.C(y)))
setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "Matrix"),
	  function(x, y = NULL) tcrossprod(.T.2.C(x), y))
setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "TsparseMatrix"),
	  function(x, y = NULL) tcrossprod(.T.2.C(x), .T.2.C(y)))


## "Matrix"
setMethod("tcrossprod", signature(x = "Matrix", y = "numLike"), 
	  function(x, y) x %*% rbind(y,deparse.level=0))
setMethod("tcrossprod", signature(x = "numLike", y = "Matrix"), .v.Mt)
setMethod("tcrossprod", signature(x = "Matrix", y = "matrix"),
	  function(x, y = NULL) tcrossprod(x, Matrix(y)))
setMethod("tcrossprod", signature(x = "matrix", y = "Matrix"),
	  function(x, y = NULL) tcrossprod(Matrix(x), y))

## sparseVector
setMethod("tcrossprod", signature(x = "Matrix", y = "sparseVector"), 
	  function(x, y) tcrossprod(x, .sparseV2Mat(y)))
setMethod("tcrossprod", signature(x = "sparseVector", y = "Matrix"), .v.Mt)
setMethod("tcrossprod", signature(x = "sparseMatrix", y = "sparseVector"),
	  function(x, y) tcrossprod(x, .sparseV2Mat(y)))
setMethod("tcrossprod", signature(x = "sparseVector", y = "sparseMatrix"), .v.Mt)
setMethod("tcrossprod", signature(x = "sparseVector", y = "sparseVector"),
	  function(x, y) .sparseV2Mat(x) %*%
          spV2M(y, nrow=1L, ncol=length(y), check=FALSE))
setMethod("tcrossprod", signature(x = "sparseVector", y = "missing"),
	  ## could be speeded: spV2M(x, *) called twice with different ncol/nrow
	  function(x, y=NULL) .sparseV2Mat(x) %*%
	  spV2M(x, nrow=1L, ncol=length(x), check=FALSE))


## cheap fallbacks
setMethod("tcrossprod", signature(x = "Matrix", y = "Matrix"),
	  function(x, y = NULL) x %*% t(y))
setMethod("tcrossprod", signature(x = "Matrix", y = "missing"),
	  function(x, y = NULL) x %*% t(x))
setMethod("tcrossprod", signature(x = "Matrix", y = "ANY"),
	  function(x, y = NULL) x %*% t(y))
setMethod("tcrossprod", signature(x = "ANY", y = "Matrix"),
	  function(x, y = NULL) x %*% t(y))

## Local variables:
## mode: R
## page-delimiter: "^###---"
## End:
