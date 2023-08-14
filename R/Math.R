####--- All "Math" and "Math2" group methods for all Matrix classes (incl sparseVector) ------
####	     ====	=====

## "Design-bug":  log(x, base)  has *two* arguments // ditto for  "trunc()" !!
## ---> need "log" methods "everywhere to catch 2-arg case !


### ~~~~ Math, log ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## FIXME:
## Once we start having non-virtual [iz]Matrix,
## many of these will need adjustment ...

## cum(min|max|sum|prod) return vector also for matrix arguments
Math.vecGenerics <- grep("^cum", getGroupMembers("Math"), value = TRUE)

###--------- dgeMatrix

setMethod("Math", signature(x = "dgeMatrix"), function(x)
{
    if(.Generic %in% Math.vecGenerics)
        callGeneric(x@x)
    else {
        x@x <- callGeneric(x@x)
        x
    }
})

setMethod("log", "dgeMatrix", function(x, base = exp(1))
{
    x@x <- log(x@x, base)
    x
})

###--------- ddenseMatrix

## Used for dt[rp]Matrix, ds[yp]Matrix (and subclasses, e.g., dpo*, cor*)
##' _only_, as dgeMatrix has its own method above

setMethod("Math", signature(x = "ddenseMatrix"), function(x)
{
    if(.Generic %in% Math.vecGenerics)
        ## Result is a vector
        return(callGeneric(.M2gen(x, ".")@x))
    cld <- getClassDef(class(x))
    if(extends(cld, "symmetricMatrix")) {
        ## Argument and result are symmetricMatrix
        if((po <- extends(cld, "dpoMatrix")) || extends(cld, "dppMatrix"))
            ## But result is _not_ positive definite!
            x <- as(x, if(po) "dsyMatrix" else "dspMatrix")
        x@x <- callGeneric(x@x)
        x@factors <- list()
        x
    } else if(is0(callGeneric(0))) {
        ## Argument and result are triangularMatrix
        if(extends(cld, "MatrixFactorization"))
            ## But result is _not_ a factor or correlation
            x <- as(x, if(.isPacked(x)) "dtpMatrix" else "dtrMatrix")
        x@x <- callGeneric(x@x)
        if(x@diag != "N" && isN1(f1 <- callGeneric(1)))
            diag(x) <- f1
        x
    } else {
        ## Argument is triangularMatrix, result is generalMatrix
        callGeneric(.M2gen(x, "."))
    }
})

## "log" with *two* arguments
setMethod("log", signature(x = "ddenseMatrix"), function(x, base = exp(1))
{
    cld <- getClassDef(class(x))
    if(extends(cld, "symmetricMatrix")) {
        ## Argument and result are symmetricMatrix
        if((po <- extends(cld, "dpoMatrix")) || extends(cld, "dppMatrix"))
            ## But result is _not_ positive definite
            x <- as(x, if(po) "dsyMatrix" else "dspMatrix")
        x@x <- log(x@x, base)
        x@factors <- list()
        x
    } else {
        ## Argument is triangularMatrix, result is generalMatrix
        log(.M2gen(x, "."), base)
    }
})

###--------- denseMatrix

setMethod("Math", signature(x = "denseMatrix"),
	  function(x) callGeneric(.M2kind(x, "d")))

setMethod("log", signature(x = "denseMatrix"),
          function(x, base = exp(1)) log(.M2kind(x, "d"), base))

###--------- CsparseMatrix

setMethod("Math", signature(x = "CsparseMatrix"), function(x)
{
    if(.Generic %in% Math.vecGenerics)
        ## Result is a vector
        return(callGeneric(.M2m(x)))
    if(isN0(callGeneric(0)))
        ## Result is a denseMatrix
        return(callGeneric(.sparse2dense(x)))
    ## Result preserves sparseness and structure (symmetric, triangular)
    cld <- getClassDef(cl <- class(x))
    if(isN1(callGeneric(1)))
        x <- .Call(R_sparse_diag_U2N, x)
    if(extends(cld, "nsparseMatrix")) {
        ## No 'x' slot
        r <- rep.int(callGeneric(1), length(x@i))
    } else {
        r <- callGeneric(x@x)
        if(typeof(r) == typeof(x@x)) {
            x@x <- r
            return(x)
        }
    }
    ## e.g., abs( <lgC> ) -> dgC
    y <- new(`substr<-`(MatrixClass(cl, cld), 1L, 1L, "d"))
    y@x <- as.double(r)
    nms <- slotNames(cld)
    for(nm in nms[nms != "x"])
        slot(y, nm) <- slot(x, nm)
    y
}) ## {Math}

setMethod("log", signature(x = "CsparseMatrix"),
          function(x, base = exp(1)) log(.sparse2dense(x), base))

###--------- diagonalMatrix

setMethod("Math", signature(x = "diagonalMatrix"), function(x)
{
    if(.Generic %in% Math.vecGenerics)
        ## Result is a vector
        return(callGeneric(.M2m(x)))
    unit <- x@diag != "N"
    r <- callGeneric(if(unit) 1 else x@x)
    if(isN0(f0 <- callGeneric(0))) {
        ## Result is dense, symmetric
        ## MJ: hmm ... what if the 'Dimnames' are asymmetric?
        y <- new("dspMatrix")
        n <- (y@Dim <- x@Dim)[1L]
        y@Dimnames <- symmDN(x@Dimnames)
        y@x <- rep.int(f0, 0.5 * n * (n + 1))
        if(n > 0L)
            diag(y) <- r
        y
    } else if(typeof(r) == typeof(x@x)) {
        ## Result is diagonal ... modify 'x'
        if(!unit) {
            x@x <- r
        } else if(isN1(r)) {
            x@x <- rep.int(r, x@Dim[1L])
            x@diag <- "N"
        }
        x
    } else {
        ## Result is diagonal ... modify new()
        y <- new("ddiMatrix")
        y@Dim <- x@Dim
        y@Dimnames <- x@Dimnames
        if(!unit)
            y@x <- as.double(r)
        else if(isN1(r))
            y@x <- rep.int(as.double(r), x@Dim[1L])
        else
            y@diag <- "U"
        y
    }
}) ## {Math}

setMethod("log", "diagonalMatrix", function(x, base = exp(1))
{
    ## Result is dense, symmetric
    ## MJ: hmm ... what if the 'Dimnames' are asymmetric?
    y <- new("dspMatrix")
    n <- (y@Dim <- x@Dim)[1L]
    y@Dimnames <- symmDN(x@Dimnames)
    y@x <- rep.int(-Inf, 0.5 * n * (n + 1))
    if(n > 0L)
        diag(y) <- if(x@diag == "N") log(x@x, base) else 0
    y
})

###--------- sparseMatrix

setMethod("Math", signature(x = "sparseMatrix"),
	  function(x) callGeneric(as(x, "CsparseMatrix")))

setMethod("log", signature(x = "sparseMatrix"),
          function(x, base = exp(1)) log(as(x, "CsparseMatrix"), base))

###--------- sparseVector

setMethod("Math", signature(x = "sparseVector"), function(x)
{
    if(.Generic %in% Math.vecGenerics || isN0(callGeneric(0)))
        ## Result is a (traditional) vector
        return(callGeneric(sp2vec(x)))
    ## Result is a sparseVector
    cld <- getClassDef(class(x))
    if(extends(cld, "dsparseVector")) {
        x@x <- callGeneric(x@x)
        x
    } else {
        y <- new("dsparseVector")
        y@x <-
            if(extends(cld, "nsparseVector"))
                rep.int(callGeneric(1), length(x@i))
            else callGeneric(x@x)
        y@i <- x@i
        y@length <- x@length
        y
    }
})

setMethod("log", "sparseVector", function(x, base = exp(1))
{
    lx <- rep.int(-Inf, x@length)
    if(length(x@i) > 0L)
        lx[x@i] <- if(is(x, "nsparseVector")) 0 else log(x@x, base)
    lx
})


### ~~~~ Math2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## NB: For round(), signif(), we have that Generic(u, k) |-> u
##     for all u in {0,1}, for all k, implying that "structure"
##     is invariant ... hence minimal "cases" are needed here

setMethod("Math2", signature(x = "dMatrix"),
          function(x, digits) {
              x@x <- callGeneric(x@x, digits = digits)
              x
          })

## As above, but first coercing to dMatrix:
setMethod("Math2", signature(x = "Matrix"),
	  function(x, digits) {
	      x <- as(x, "dMatrix")
	      x@x <- callGeneric(x@x, digits = digits)
	      x
	  })

setMethod("Math2", signature(x = "dsparseVector"),
          function(x, digits) {
              x@x <- callGeneric(x@x, digits = digits)
              x
          })

## As above, but first coercing to dsparseVector:
setMethod("Math2", signature(x = "sparseVector"),
	  function(x, digits) {
	      x <- as(x, "dsparseVector")
	      x@x <- callGeneric(x@x, digits = digits)
	      x
	  })


## ~~~~ Not group generic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("zapsmall", signature(x = "dMatrix"),
          function(x, digits = getOption("digits")) {
              x@x <- zapsmall(x@x, digits)
              x
          })
