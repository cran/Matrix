## METHODS FOR CLASS: Matrix (virtual)
## mother class containing all matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ~~~~ COERCIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs(   "ANY", "Matrix", function(from) Matrix(as.matrix(from)))
setAs("matrix", "Matrix", function(from) Matrix(          from ))

## Need 'base' functions calling as.*() to dispatch to our S4 methods:
as.vector.Matrix <- function(x, mode = "any") as.vector(as(x, "matrix"), mode)
as.matrix.Matrix <- function(x, ...)                    as(x, "matrix")
 as.array.Matrix <- function(x, ...)                    as(x, "matrix")

## FIXME: compare with methods for subclasses; avoid duplication

setMethod("as.vector", signature(x = "Matrix"), as.vector.Matrix)
setMethod("as.matrix", signature(x = "Matrix"), as.matrix.Matrix)
setMethod( "as.array", signature(x = "Matrix"),  as.array.Matrix)

setMethod("as.logical", signature(x = "Matrix"),
          function(x, ...) as.logical(as.vector(x)))
setMethod("as.numeric", signature(x = "Matrix"),
          function(x, ...) as.numeric(as.vector(x)))

setAs("Matrix",  "vector", function(from)  as.vector(as(from, "matrix")))
setAs("Matrix", "logical", function(from) as.logical(as(from, "matrix")))
setAs("Matrix", "integer", function(from) as.integer(as(from, "matrix")))
setAs("Matrix", "numeric", function(from) as.numeric(as(from, "matrix")))
setAs("Matrix", "complex", function(from) as.complex(as(from, "matrix")))


## ~~~~ CONSTRUCTORS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Matrix <- function(data = NA, nrow = 1, ncol = 1, byrow = FALSE,
                   dimnames = NULL, sparse = NULL,
                   doDiag = TRUE, forceCheck = FALSE)
{
    i.M <- i.sM <- i.dM <- i.sV <- i.m <- FALSE
    mnrow <- missing(nrow)
    mncol <- missing(ncol)
    if(isS4(data)) {
        cld <- getClassDef(class(data))
        i.M <- extends(cld, "Matrix")
        if(i.M) {
            i.sM <- extends(cld, "sparseMatrix")
            i.dM <- i.sM && extends(cld, "diagonalMatrix")
        } else if(extends(cld, "sparseVector")) {
            ## need to transmit missingness to 'spV2M'
            call. <- quote(spV2M(x = data, nrow =, ncol =, byrow = byrow))
            if(!mnrow)
                call.[[3L]] <- quote(nrow)
            if(!mncol)
                call.[[4L]] <- quote(ncol)
            data <- eval(call.)
            i.M <- i.sM <- i.sV <- forceCheck <- TRUE
        }
    } else {
        i.m <- is.matrix(data)
    }
    if(!i.M) {
        ## validate non-Matrix 'data', throwing type errors _early_
        if(is.object(data)) {
            if(i.m)
                class(data) <- NULL # retaining 'dim'
            else
                data <- as.vector(data)
        }
        mode. <- mode(data)
        kind <- switch(mode., numeric = "d", logical = "l",
                       stop("invalid 'data'"))
    }
    if(i.M || i.m) {
        ## 'data' is a Matrix or a numeric or logical matrix
        ## without a 'class' attribute
        if(!i.sV && !(mnrow && mncol && missing(byrow)))
            warning("'nrow', 'ncol', 'byrow' disregarded for [mM]atrix 'data'")
        if(!is.null(dimnames))
            dimnames(data) <- dimnames
        if(is.null(sparse))
            sparse <- sparseDefault(data)
        if(i.M) {
            ## return early in these cases:
            if(i.dM)
                ## !doDiag has been documented to result in a coercion to
                ## symmetricMatrix; we must use diag2*() below because the
                ## "usual" as(<diagonalMatrix>, "(Csparse|unpacked)Matrix")
                ## inherits from triangularMatrix, _not_ symmetricMatrix
                return(if(doDiag)
                           data
                       else if(sparse)
                           .diag2sparse(data, ".sC", "U", FALSE)
                       else .diag2dense(data, ".sy", "U"))
            if(!forceCheck)
                return(if(i.sM == sparse)
                           data
                       else if(sparse)
                           as(data, "CsparseMatrix")
                       else as(data, "unpackedMatrix"))
        }
    } else {
        ## 'data' is a numeric or logical vector or non-matrix array
        ## without a 'class' attribute
        if(length(data) == 1L && !is.na(data) && data == 0 &&
           (is.null(sparse) || sparse)) {
            ## Matrix(0, ...): sparseMatrix unless sparse=FALSE
            ## MJ: we should _try_ to behave as R's do_matrix()
            ##     in the edge cases ... integer overflow is "OK"
            ##     since anyNA(Dim) is caught by validity methods
            if(mnrow == mncol) {
                nrow <- as.integer(nrow)
                ncol <- as.integer(ncol)
            } else if(mnrow) {
                ncol <- as.integer(ncol)
                if(ncol == 0L)
                    stop("data is too long")
                nrow <- as.integer(ceiling(1 / ncol))
            } else {
                nrow <- as.integer(nrow)
                if(nrow == 0L)
                    stop("data is too long")
                ncol <- as.integer(ceiling(1 / nrow))
            }
            square <- nrow == ncol
            if(is.null(dimnames))
                dimnames <- list(NULL, NULL)
            if(square && doDiag)
                return(new(paste0(kind, "diMatrix"),
                           Dim = c(nrow, ncol),
                           Dimnames = dimnames,
                           x = vector(mode., nrow)))
            data <- new(paste0(kind, if(square) "s" else "g", "CMatrix"),
                        Dim = c(nrow, ncol),
                        Dimnames = dimnames,
                        p = integer(ncol + 1))
            i.M <- i.sM <- sparse <- TRUE
        } else {
            ## usual case: vector|array->matrix
            data <- .External(Mmatrix,
                              data, nrow, ncol, byrow, dimnames, mnrow, mncol)
            if(is.null(sparse))
                sparse <- sparseDefault(data)
            i.m <- TRUE
        }
    }

    ## 'data' is a Matrix (but _not_ a diagonalMatrix) or a
    ## numeric or logical matrix without a 'class' attribute
    if(doDiag && isDiagonal(data))
        ## as(<[mM]atrix>, "diagonalMatrix") uses check = TRUE (a waste)
        return(forceDiagonal(data))
    if(i.m || i.sM != sparse) {
        data <- as(data, if(sparse) "CsparseMatrix" else "unpackedMatrix")
        if(i.m)
            ## as(<matrix>, "CsparseMatrix"), as(<matrix>, "unpackedMatrix")
            ## already check for symmetric, triangular structure
            return(data)
    }
    if(!is(data, "generalMatrix"))
        data
    else if(isSymmetric(data))
        forceSymmetric(data)
    else if(!(it <- isTriangular(data)))
        data
    else if(attr(it, "kind") == "U")
        triu(data)
    else tril(data)
}


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dim", signature(x = "Matrix"),
          function(x) x@Dim)

setMethod("length", "Matrix",
          function(x) prod(x@Dim))

setMethod("dimnames", signature(x = "Matrix"),
          function(x) x@Dimnames)

setMethod("dimnames<-", signature(x = "Matrix", value = "list"),
          function(x, value) {
              x@Dimnames <- fixupDN.if.valid(value, x@Dim)
              x
          })

setMethod("dimnames<-", signature(x = "Matrix", value = "NULL"),
          function(x, value) {
              x@Dimnames <- list(NULL, NULL)
              x
          })

setMethod("dimnames<-", signature(x = "compMatrix", value = "list"),
          function(x, value) {
              if(length(x@factors))
                  x@factors <- list()
              x@Dimnames <- fixupDN.if.valid(value, x@Dim)
              x
          })

setMethod("dimnames<-", signature(x = "compMatrix", value = "NULL"),
          function(x, value) {
              if(length(x@factors))
                  x@factors <- list()
              x@Dimnames <- list(NULL, NULL)
              x
          })

setMethod("unname", signature(obj = "Matrix"),
          function(obj, force = FALSE) {
              obj@Dimnames <- list(NULL, NULL)
              obj
          })

setMethod("drop", signature(x = "Matrix"),
          function(x) if(any(x@Dim == 1L)) drop(as(x, "matrix")) else x)

## These work nicely as long as methods are defined for '[' :
setMethod("head", signature(x = "Matrix"),
          head.matrix)
setMethod("tail", signature(x = "Matrix"),
          tail.matrix)
setMethod("diff", signature(x = "Matrix"),
          ## Mostly cut and paste of 'base::diff.default' :
          function(x, lag = 1L, differences = 1L, ...) {
              if(length(lag) != 1L || length(differences) > 1L ||
                  lag < 1L || differences < 1L)
                  stop("'lag' and 'differences' must be integers >= 1")
              if(lag * differences >= x@Dim[1L])
                  return(x[0L])
              i1 <- -seq_len(lag)
              for(i in seq_len(differences)) {
                  m <- x@Dim[1L]
                  x <- x[i1, , drop = FALSE] -
                      x[-m:-(m - lag + 1L), , drop = FALSE]
              }
              x
          })

if(FALSE) { ## still does not work for c(1, Matrix(2))
## For the same reason (and just in case) also do both S3 and S4 here:
c.Matrix <- function(...) unlist(lapply(list(...), as.vector))
## NB: Must use   signature  '(x, ..., recursive = FALSE)' :
setMethod("c", "Matrix", function(x, ..., recursive) c.Matrix(x, ...))
## The above is not sufficient for  c(NA, 3:2, <Matrix>, <matrix>)
setMethod("c", "numMatrixLike", function(x, ..., recursive) c.Matrix(x, ...))
}# not yet

## We want to use all.equal.numeric() *and* make sure that uses
## not just base::as.vector but the generic with our methods:
all.equal_num <- base::all.equal.numeric
##               ^^^^^^<R>/src/library/base/R/all.equal.R
environment(all.equal_num) <- environment() # our namespace

all.equal_Mat <- function(target, current, check.attributes = TRUE,
                          factorsCheck = FALSE, ...)
{
    msg <- attr.all_Mat(target, current, check.attributes=check.attributes,
                        factorsCheck=factorsCheck, ...)
    if(is.list(msg)) msg[[1]]
    else .a.e.comb(msg,
                   all.equal_num(as.vector(target), as.vector(current),
                                 check.attributes=check.attributes, ...))
}

## The all.equal() methods for dense matrices (and fallback):
setMethod("all.equal", c(target = "Matrix", current = "Matrix"),
          all.equal_Mat)
setMethod("all.equal", c(target = "Matrix", current = "ANY"),
          all.equal_Mat)
setMethod("all.equal", c(target = "ANY", current = "Matrix"),
          all.equal_Mat)
rm(all.equal_Mat)
## -> ./sparseMatrix.R, ./sparseVector.R  have specific methods

### "[<-" : -----------------

## A[ ij ] <- value,  where ij is (i,j) 2-column matrix :
## ----------------
## The cheap general method, now only used for "pMatrix","indMatrix"
## sparse all use  .TM.repl.i.mat()
## NOTE:  need '...' below such that setMethod() does
##	  not use .local() such that nargs() will work correctly:
.M.repl.i.2col <- function (x, i, j, ..., value) {
    nA <- nargs()
    if(nA == 3) { ##  M [ cbind(ii,jj) ] <- value  or M [ Lmat ] <- value
        if(!is.integer(nc <- ncol(i)))
            stop(".M.repl.i.2col(): 'i' has no integer column number;\n should never happen; please report")
        else if(!is.numeric(i) || nc != 2)
            stop("such indexing must be by logical or 2-column numeric matrix")
        if(is.logical(i)) {
            message(".M.repl.i.2col(): drop 'matrix' case ...")
            ## c(i) : drop "matrix" to logical vector
            return( callGeneric(x, i=c(i), value=value) )
        }
        if(!is.integer(i)) storage.mode(i) <- "integer"
        if(any(i < 0))
            stop("negative values are not allowed in a matrix subscript")
        if(anyNA(i))
            stop("NAs are not allowed in subscripted assignments")
        if(any(i0 <- (i == 0))) # remove them
            i <- i[ - which(i0, arr.ind = TRUE)[,"row"], ]
        ## now have integer i >= 1
        m <- nrow(i)
        ## mod.x <- .type.kind[.M.kind(x)]
        if(length(value) > 0 && m %% length(value) != 0)
            warning("number of items to replace is not a multiple of replacement length")
        ## recycle:
        value <- rep_len(value, m)
        i1 <- i[,1]
        i2 <- i[,2]
        if(m > 2)
            message("m[ <ij-matrix> ] <- v: inefficiently treating single elements")
        ## inefficient -- FIXME -- (also loses "symmetry" unnecessarily)
        for(k in seq_len(m))
            x[i1[k], i2[k]] <- value[k]
        x
    } else
        stop(gettextf("nargs() = %d.  Extraneous illegal arguments inside '[ .. ]' ?",
                      nA),
             domain = NA)
}

setReplaceMethod("[",
                 signature(x = "Matrix", i = "matrix", j = "missing",
                           value = "replValue"),
                 .M.repl.i.2col)

## Three catch-all methods ... would be very inefficient for sparse*
## --> extra methods in ./sparseMatrix.R
setReplaceMethod("[",
                 signature(x = "Matrix", i = "missing", j = "ANY",
                           value = "Matrix"),
                 function(x, i, j, ..., value)
                     callGeneric(x=x, , j=j, value = as.vector(value)))

setReplaceMethod("[",
                 signature(x = "Matrix", i = "ANY", j = "missing",
                           value = "Matrix"),
                 function(x, i, j, ..., value)
                     if(nargs() == 3)
                         callGeneric(x=x, i=i, value = as.vector(value))
                     else
                         callGeneric(x=x, i=i, , value = as.vector(value)))

setReplaceMethod("[",
                 signature(x = "Matrix", i = "ANY", j = "ANY",
                           value = "Matrix"),
                 function(x, i, j, ..., value)
                     callGeneric(x=x, i=i, j=j, value = as.vector(value)))


setReplaceMethod("[",
                 signature(x = "Matrix", i = "missing", j = "ANY",
                           value = "matrix"),
                 function(x, i, j, ..., value)
                     callGeneric(x=x, , j=j, value = c(value)))

setReplaceMethod("[",
                 signature(x = "Matrix", i = "ANY", j = "missing",
                           value = "matrix"),
                 function(x, i, j, ..., value)
                     if(nargs() == 3)
                         callGeneric(x=x, i=i, value = c(value))
                     else
                         callGeneric(x=x, i=i, , value = c(value)))

setReplaceMethod("[",
                 signature(x = "Matrix", i = "ANY", j = "ANY",
                           value = "matrix"),
                 function(x, i, j, value)
                     callGeneric(x=x, i=i, j=j, value = c(value)))


##  M [ <lMatrix> ] <- value; used notably for x = "CsparseMatrix"
.repl.i.lDMat <- function (x, i, j, ..., value)
    `[<-`(x, i=which(as.vector(i)), value=value)
setReplaceMethod("[",
                 signature(x = "Matrix", i = "ldenseMatrix", j = "missing",
                           value = "replValue"),
                 .repl.i.lDMat)
setReplaceMethod("[",
                 signature(x = "Matrix", i = "ndenseMatrix", j = "missing",
                           value = "replValue"),
                 .repl.i.lDMat)
rm(.repl.i.lDMat)

.repl.i.lSMat <- function (x, i, j, ..., value)
    `[<-`(x, i=which(as(i, "sparseVector")), value=value)
setReplaceMethod("[",
                 signature(x = "Matrix", i = "lsparseMatrix", j = "missing",
                           value = "replValue"),
                 .repl.i.lSMat)
setReplaceMethod("[",
                 signature(x = "Matrix", i = "nsparseMatrix", j = "missing",
                           value = "replValue"),
                 .repl.i.lSMat)
rm(.repl.i.lSMat)

## (ANY,ANY,ANY) is used when no `real method' is implemented :
setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY",
                                value = "ANY"),
                 function (x, i, j, value) {
                     if(!is.atomic(value))
                         stop(gettextf("RHS 'value' (class %s) matches 'ANY', but must match matrix class %s",
                                       class(value), class(x)),
                              domain = NA)
                     else stop("not-yet-implemented 'Matrix[<-' method")
                 })
