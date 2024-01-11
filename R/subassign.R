## METHODS FOR GENERIC: [<-
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## GOAL: automate method definitions and eventually replace the ones
##       collected below ...
##
##       need to write C-level functions
##
##             *_subassign_1ary    (x, i,    value)
##             *_subassign_1ary_mat(x, i,    value)
##             *_subassign_2ary    (x, i, j, value)
##
##       for * = unpackedMatrix,packedMatrix,
##               CsparseMatrix,RsparseMatrix,TsparseMatrix,
##               diagonalMatrix,indMatrix

if(FALSE) { # TODO
.subassign.invalid <- function(value) {
    if(is.object(i))
        gettextf("invalid subassignment value class \"%s\"", class(i)[1L])
    else gettextf("invalid subassignment value type \"%s\"", typeof(i))
}

.subassign.1ary <- function(x, i, value) {

}

..subassign.1ary <- function(x, i, value) {

}

.subassign.1ary.mat <- function(x, i, value) {

}

..subassign.1ary.mat <- function(x, i, value) {

}

.subassign.2ary <- function(x, i, j, value) {

}

..subassign.2ary <- function(x, i, j, value) {

}

setMethod("[<-", signature(x = "Matrix", i = "missing", j = "missing",
                           value = "ANY"),
          function(x, i, j, ..., value) {
              if(missing(value))
                  stop("missing subassignment value")
              na <- nargs()
              if(na <= 4L)
                  ## x[], x[, ] <- value
                  .subassign.1ary(x, , value)
              else
                  ## x[, , ], etc. <- value
                  stop("incorrect number of dimensions")
          })

setMethod("[<-", signature(x = "Matrix", i = "index", j = "missing",
                           value = "ANY"),
          function(x, i, j, ..., value) {
              if(missing(value))
                  stop("missing subassignment value")
              na <- nargs()
              if(na == 3L)
                  ## x[i=] <- value
                  .subassign.1ary(x, i, value)
              else if(na == 4L)
                  ## x[i=, ], x[, i=] <- value
                  .subassign.2ary(x, i, , value)
              else
                  ## x[i=, , ], etc. <- value
                  stop("incorrect number of dimensions")
          })

setMethod("[<-", signature(x = "Matrix", i = "missing", j = "index",
                           value = "ANY"),
          function(x, i, j, ..., value) {
              if(missing(value))
                  stop("missing subassignment value")
              na <- nargs()
              if(na == 3L)
                  ## x[j=] <- value
                  .subassign.1ary(x, j, value)
              else if(na == 4L)
                  ## x[j=, ], x[, j=] <- value
                  .subassign.2ary(x, , j, value)
              else
                  ## x[, j=, ], etc. <- value
                  stop("incorrect number of dimensions")
          })

setMethod("[<-", signature(x = "Matrix", i = "index", j = "index",
                           value = "ANY"),
          function(x, i, j, ..., value) {
              if(missing(value))
                  stop("missing subassignment value")
              na <- nargs()
              if(na == 4L)
                  ## x[i=, j=], x[j=, i=] <- value
                  .subassign.2ary(x, i, j, value)
              else
                  ## x[i=, j=, ], etc. <- value
                  stop("incorrect number of dimensions")
          })

for(.cl in c("matrix", "nMatrix", "lMatrix"))
setMethod("[<-", signature(x = "Matrix", i = .cl, j = "missing",
                           value = "ANY"),
          function(x, i, j, ..., value) {
              if(missing(value))
                  stop("missing subassignment value")
              na <- nargs()
              if(na == 3L)
                  ## x[i=] <- value
                  .subassign.1ary.mat(x, i, value)
              else if(na == 4L)
                  ## x[i=, ], x[, i=] <- value
                  .subassign.2ary(x, i, , value)
              else
                  ## x[i=, , ], etc. <- value
                  stop("incorrect number of dimensions")
          })
rm(.cl)

setMethod("[<-", signature(x = "Matrix", i = "NULL", j = "ANY",
                           value = "ANY"),
          function(x, i, j, ..., value) {
              i <- integer(0L)
              callGeneric()
          })

setMethod("[<-", signature(x = "Matrix", i = "ANY", j = "NULL",
                           value = "ANY"),
          function(x, i, j, ..., value) {
              j <- integer(0L)
              callGeneric()
          })

setMethod("[<-", signature(x = "Matrix", i = "NULL", j = "NULL",
                           value = "ANY"),
          function(x, i, j, ..., value) {
              i <- integer(0L)
              j <- integer(0L)
              callGeneric()
          })

setMethod("[<-", signature(x = "sparseVector", i = "missing", j = "missing",
                           value = "ANY"),
          function(x, i, j, ..., value) {
              if(missing(value))
                  stop("missing subassignment value")
              if(nargs() > 3L)
                  stop("incorrect number of dimensions")
              if(isS4(value)) {
                  if(!.isVector(value))
                      stop(.subassign.invalid(value), domain = NA)
              } else
                  value <- switch(typeof(value),
                                  "logical" =,
                                  "integer" =,
                                  "double" =,
                                  "complex" = .m2V(value),
                                  stop(.subassign.invalid(value), domain = NA))
              n.x <- length(x)
              n.value <- length(value)
              if(n.x > 0L && n.value == 0L)
                  stop("replacement has length zero")
              k.x <- .M2kind(x)
              k.value <- .M2kind(value)
              if(k.x != k.value) {
                  map <- `names<-`(1:5, c("n", "l", "i", "d", "z"))
                  if(map[[k.value]] < map[[k.x]])
                      value <- .V2kind(value, k.x)
              }
              if(n.value == 0L)
                  return(value)
              if(n.x %% n.value != 0L)
                  warning("number of items to replace is not a multiple of replacement length")
              .V.rep.len(value, n.x)
          })

setMethod("[<-", signature(x = "sparseVector", i = "index", j = "missing",
                           value = "ANY"),
          function(x, i, j, ..., value) {
              if(missing(value))
                  stop("missing subassignment value")
              if(nargs() > 3L)
                  stop("incorrect number of dimensions")
              if(isS4(value)) {
                  if(!.isVector(value))
                      stop(.subassign.invalid(value), domain = NA)
              } else
                  value <- switch(typeof(value),
                                  "logical" =,
                                  "integer" =,
                                  "double" =,
                                  "complex" = .m2V(value),
                                  stop(.subassign.invalid(value), domain = NA))
              switch(typeof(i),
                     "logical" = {},
                     "integer" = {},
                     "double" = {},
                     stop(.subscript.invalid(value), domain = NA))
              k.x <- .M2kind(x)
              k.value <- .M2kind(value)
              if(k.x != k.value) {
                  map <- `names<-`(1:5, c("n", "l", "i", "d", "z"))
                  if(map[[k.x]] < map[[k.value]])
                      x <- .V2kind(x, k.value)
              }
              ## TODO
          })

setMethod("[<-", signature(x = "sparseVector", i = "nsparseVector", j = "missing",
                           value = "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              if(missing(value))
                  stop("missing subassignment value")
              if(nargs() > 3L)
                  stop("incorrect number of dimensions")
              x[.subscript.recycle(i, length(x), TRUE)] <- value
              x
          })

setMethod("[<-", signature(x = "sparseVector", i = "lsparseVector", j = "missing",
                           value = "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              if(missing(value))
                  stop("missing subassignment value")
              if(nargs() > 3L)
                  stop("incorrect number of dimensions")
              x[.subscript.recycle(i, length(x), FALSE)] <- value
              x
          })

setMethod("[<-", signature(x = "sparseVector", i = "NULL", j = "ANY",
                           value = "ANY"),
          function(x, i, j, ..., value) {
              i <- integer(0L)
              callGeneric()
          })
} # TODO

## ==== Matrix =========================================================

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


## ==== denseMatrix ====================================================

## x[] <- value :
setReplaceMethod("[", signature(x = "denseMatrix", i = "missing", j = "missing",
                                value = "ANY"),## double/logical/...
                 function (x, value) {
                     x <- .M2gen(x)
                     x@x[] <- value
                     validObject(x)# check if type and lengths above match
                     x
                 })

## FIXME: 1) These are far from efficient
## -----
setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "missing",
                                value = "replValue"),
                 function (x, i, j, ..., value) {
                     r <- as(x, "matrix")
                     ## message("`[<-` with nargs()= ",nargs())
                     if((na <- nargs()) == 3)
                         r[i] <- value
                     else if(na == 4)
                         r[i, ] <- value
                     else stop(gettextf("invalid nargs()= %d", na), domain=NA)
                     .m2dense(r, paste0(.M.kind(x), "ge"))
                 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "missing", j = "index",
                                value = "replValue"),
                 function (x, i, j, ..., value) {
                     r <- as(x, "matrix")
                     r[, j] <- value
                     .m2dense(r, paste0(.M.kind(x), "ge"))
                 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "index",
                                value = "replValue"),
                 function (x, i, j, ..., value) {
                     r <- as(x, "matrix")
                     r[i, j] <- value
                     as_denseClass(r, class(x)) ## was as(r, class(x))
                 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "matrix",  # 2-col.matrix
                                j = "missing", value = "replValue"),
                 function(x, i, j, ..., value) {
                     r <- as(x, "matrix")
                     r[ i ] <- value
                     .m2dense(r, paste0(.M.kind(x), "ge"))
                 })


## ==== sparseMatrix ===================================================

## x[] <- value :
setReplaceMethod("[", signature(x = "sparseMatrix", i = "missing", j = "missing",
                                value = "ANY"),## double/logical/...
                 function (x, i, j,..., value) {
                     if(all0(value)) { # be faster
                         cld <- getClassDef(class(x))
                         x <- diagU2N(x, cl = cld)
                         for(nm in intersect(nsl <- names(cld@slots),
                                             c("x", "i","j", "factors")))
                             length(slot(x, nm)) <- 0L
                         if("p" %in% nsl)
                             x@p <- rep.int(0L, ncol(x)+1L)
                     } else {
                         ## typically non-sense: assigning to full sparseMatrix
                         x[TRUE] <- value
                     }
                     x
                 })

## Do not use as.vector() (see ./Matrix.R ) for sparse matrices :
setReplaceMethod("[", signature(x = "sparseMatrix", i = "missing", j = "ANY",
                                value = "sparseMatrix"),
                 function (x, i, j, ..., value)
                     callGeneric(x=x, , j=j, value=as(value, "sparseVector")))

setReplaceMethod("[", signature(x = "sparseMatrix", i = "ANY", j = "missing",
                                value = "sparseMatrix"),
                 function (x, i, j, ..., value)
                     if(nargs() == 3)
                         callGeneric(x=x, i=i, value=as(value, "sparseVector"))
                     else
                         callGeneric(x=x, i=i, , value=as(value, "sparseVector")))

setReplaceMethod("[", signature(x = "sparseMatrix", i = "ANY", j = "ANY",
                                value = "sparseMatrix"),
                 function (x, i, j, ..., value)
                     callGeneric(x=x, i=i, j=j, value=as(value, "sparseVector")))


## ==== CsparseMatrix ==================================================

## workhorse for "[<-" -- for d*, l*, and n..C-sparse matrices :
## ---------     -----
replCmat <- function (x, i, j, ..., value) {
    di <- dim(x)
    dn <- dimnames(x)
    iMi <- missing(i)
    jMi <- missing(j)
    na <- nargs()
    Matrix.message("replCmat[x,i,j,..,val] : nargs()=", na, "; ",
	       if(iMi || jMi) sprintf("missing (i,j) = (%d,%d)", iMi, jMi),
	       .M.level = 2)
    if(na == 3L) { ## vector (or 2-col) indexing M[i] <- v : includes M[TRUE] <- v or M[] <- v !
	x <- .M2T(x)
	x[i] <- value # may change class, e.g., from dtT* to dgT*
	cl.C <- sub(".Matrix$", "CMatrix", class(x))
	if(.hasSlot(x, "x") && any0(x@x))
	    ## drop all values that "happen to be 0"
	    drop0(x, is.Csparse = FALSE)
        else as_CspClass(x, cl.C)
    } else ## nargs() == 4 :
	replCmat4(x,
		  i1 = if(iMi)
                           seq.int(from = 0L, length.out = di[1L])
                       else .ind.prep2(i, 1L, di, dn),
		  i2 = if(jMi)
                           seq.int(from = 0L, length.out = di[2L])
                       else .ind.prep2(j, 2L, di, dn),
                  iMi = iMi, jMi = jMi, value = value)
} ## replCmat

replCmat4 <- function(x, i1, i2, iMi, jMi, value,
                      spV = is(value, "sparseVector")) {
    dind <- c(length(i1), length(i2)) # dimension of replacement region
    lenRepl <- prod(dind)
    lenV <- length(value)
    if(lenV == 0) {
	if(lenRepl != 0L)
	    stop("nothing to replace with")
	return(x)
    }
    ## else: lenV := length(value)	 is > 0
    if(lenRepl %% lenV != 0L)
	stop("number of items to replace is not a multiple of replacement length")
    if(lenV > lenRepl)
	stop("too many replacement values")

    clx <- class(x)
    clDx <- getClassDef(clx) # extends() , is() etc all use the class definition

    ## keep "symmetry" if changed here:
    x.sym <- extends(clDx, "symmetricMatrix")
    if(x.sym) { ## only half the indices are there..
	## using array() for large dind is a disaster...
	mkArray <- if(spV) # TODO: room for improvement
	    function(v, dim) spV2M(v, dim[1L], dim[2L]) else array
	x.sym <-
	    (dind[1L] == dind[2L] && all(i1 == i2) &&
	     (lenRepl == 1L || lenV == 1L ||
	      isSymmetric(mkArray(value, dim=dind))))
	## x.sym : result is *still* symmetric
	x <- .M2gen(x) ## but do *not* redefine clx!
    }
    else if(extends(clDx, "triangularMatrix")) {
	xU <- x@uplo == "U"
	r.tri <- ((any(dind == 1) || dind[1L] == dind[2L]) &&
		  if(xU) max(i1) <= min(i2) else max(i2) <= min(i1))
	if(r.tri) { ## result is *still* triangular
	    if(any(i1 == i2)) # diagonal will be changed
		x <- diagU2N(x) # keeps class (!)
	}
	else { # go to "generalMatrix" and (do not redefine clx!) and continue
	    x <- .M2gen(x) # was as(x, paste0(.M.kind(x), "gCMatrix"))
	}
    }
    ## Temporary hack for debugging --- remove eventually -- FIXME :
    ## see also	 MATRIX_SUBASSIGN_VERBOSE in ../src/t_Csparse_subassign.c
    if(!is.null(v <- getOption("Matrix.subassign.verbose")) && v) {
	op <- options(Matrix.verbose = 2); on.exit(options(op))
	## the "hack" to signal "verbose" to the C code:
	if(i1[1L] != 0L)
            i1[1L] <- -i1[1L]
        else warning("i1[1] == 0 ==> C-level verbosity will not happen!")
    }

    if(extends(clDx, "dMatrix")) {
	has.x <- TRUE
	x <- .Call(dCsparse_subassign,
		   if(clx %in% c("dgCMatrix", "dtCMatrix")) x
		   else .M2gen(x), # must get "dgCMatrix"
		   i1, i2,
		   as(value, "sparseVector"))
    }
    else if(extends(clDx, "lMatrix")) {
	has.x <- TRUE
	x <- .Call(lCsparse_subassign,
		   if(clx %in% c("lgCMatrix", "ltCMatrix")) x
		   else .M2gen(x), # must get "lgCMatrix"
		   i1, i2,
		   as(value, "sparseVector"))
    }
    else if(extends(clDx, "nMatrix")) {
	has.x <- FALSE
	x <- .Call(nCsparse_subassign,
		   if(clx %in% c("ngCMatrix", "ntCMatrix"))x
		   else .M2gen(x), # must get "ngCMatrix"
		   i1, i2,
		   as(value, "sparseVector"))
    }
    else if(extends(clDx, "iMatrix")) {
	has.x <- TRUE
	x <- .Call(iCsparse_subassign,
		   if(clx %in% c("igCMatrix", "itCMatrix"))x
		   else .M2gen(x), # must get "igCMatrix"
		   i1, i2,
		   as(value, "sparseVector"))
    }
    else if(extends(clDx, "zMatrix")) {
	has.x <- TRUE
	x <- .Call(zCsparse_subassign,
		   if(clx %in% c("zgCMatrix", "ztCMatrix"))x
		   else .M2gen(x), # must get "zgCMatrix"
		   i1, i2,
		   ## here we only want zsparseVector {to not have to do this in C}:
		   as(value, "zsparseVector"))
    }
    else { ## use "old" code ...
        ## does this happen ? ==>
	if(identical(Sys.getenv("USER"),"maechler"))## does it still happen? __ FIXME __
	    stop("using	 \"old code\" part in  Csparse subassignment")
        ## else
	warning("using\"old code\" part in Csparse subassignment\n >>> please report to Matrix-authors@r-project.org",
		immediate. = TRUE)

	xj <- .Call(Matrix_expand_pointers, x@p)
	sel <- (!is.na(match(x@i, i1)) &
		!is.na(match( xj, i2)))
	has.x <- "x" %in% slotNames(clDx)# === slotNames(x),
	## has.x  <==> *not* nonzero-pattern == "nMatrix"

	if(has.x && sum(sel) == lenRepl) { ## all entries to be replaced are non-zero:
	    ## need indices instead of just 'sel', for, e.g.,  A[2:1, 2:1] <- v
	    non0 <- cbind(match(x@i[sel], i1),
			  match(xj [sel], i2), deparse.level=0L)
	    iN0 <- 1L + .Call(m_encodeInd, non0, di = dind, orig1=TRUE, checkBounds=FALSE)

	    has0 <-
		if(spV) length(value@i) < lenV else any(value[!is.na(value)] == 0)
	    if(lenV < lenRepl)
		value <- rep_len(value, lenRepl)
	    ## Ideally we only replace them where value != 0 and drop the value==0
	    ## ones; FIXME: see Davis(2006) "2.7 Removing entries", p.16, e.g. use cs_dropzeros()
	    ##	     but really could be faster and write something like cs_drop_k(A, k)
	    ## v0 <- 0 == value
	    ## if (lenRepl == 1) and v0 is TRUE, the following is not doing anything
	    ##-	 --> ./Tsparse.R	and its	 replTmat()
	    ## x@x[sel[!v0]] <- value[!v0]
	    x@x[sel] <- as.vector(value[iN0])
	    if(extends(clDx, "compMatrix") && length(x@factors)) # drop cached ones
		x@factors <- list()
	    if(has0) x <- .drop0(x)

	    return(if(x.sym) as_CspClass(x, clx) else x)
	}
	## else go via Tsparse.. {FIXME: a waste! - we already have 'xj' ..}
	## and inside  Tsparse... the above i1, i2,..., sel  are *all* redone!
	## Happens too often {not anymore, I hope!}
	##
	Matrix.message("wasteful C -> T -> C in replCmat(x,i,j,v) for <sparse>[i,j] <- v")
	x <- as(x, "TsparseMatrix")
	if(iMi)
	    x[ ,i2+1L] <- value
	else if(jMi)
	    x[i1+1L, ] <- value
	else
	    x[i1+1L,i2+1L] <- value
	if(extends(clDx, "compMatrix") && length(x@factors)) # drop cached ones
	    x@factors <- list()
    }# else{ not using new memory-sparse  code
    if(has.x && any0(x@x)) ## drop all values that "happen to be 0"
	as_CspClass(drop0(x), clx)
    else as_CspClass(x, clx)
} ## replCmat4

setReplaceMethod("[", signature(x = "CsparseMatrix", i = "index", j = "missing",
                                value = "replValue"),
                 replCmat)

setReplaceMethod("[", signature(x = "CsparseMatrix", i = "missing", j = "index",
                                value = "replValue"),
                 replCmat)

setReplaceMethod("[", signature(x = "CsparseMatrix", i = "index", j = "index",
				value = "replValue"),
                 replCmat)

### When the RHS 'value' is  a sparseVector, now can use  replCmat  as well
setReplaceMethod("[", signature(x = "CsparseMatrix", i = "missing", j = "index",
				value = "sparseVector"),
		 replCmat)

setReplaceMethod("[", signature(x = "CsparseMatrix", i = "index", j = "missing",
				value = "sparseVector"),
		 replCmat)

setReplaceMethod("[", signature(x = "CsparseMatrix", i = "index", j = "index",
				value = "sparseVector"),
		 replCmat)
rm(replCmat)

## A[ ij ] <- value,  where ij is (i,j) 2-column matrix
setReplaceMethod("[", signature(x = "CsparseMatrix", i = "matrix", j = "missing",
				value = "replValue"),
		 function(x, i, j, ..., value)
		 ## goto Tsparse modify and convert back:
		 as(.TM.repl.i.mat(as(x, "TsparseMatrix"), i=i, value=value),
		    "CsparseMatrix"))
## more in ./sparseMatrix.R (and ./Matrix.R )

setReplaceMethod("[", signature(x = "CsparseMatrix", i = "Matrix", j = "missing",
				value = "replValue"),
		 function(x, i, j, ..., value)
		 ## goto Tsparse modify and convert back:
		 as(.TM.repl.i.mat(as(x, "TsparseMatrix"), i=i, value=value),
		    "CsparseMatrix"))


## ==== RsparseMatrix ==================================================

setReplaceMethod("[", signature(x = "RsparseMatrix", i = "index", j = "missing",
				value = "replValue"),
		 function (x, i, j, ..., value)
		 replTmat(.M2T(x), i=i, , value=value))

setReplaceMethod("[", signature(x = "RsparseMatrix", i = "missing", j = "index",
				value = "replValue"),
		 function (x, i, j, ..., value)# extra " , ": want nargs() == 4
		 replTmat(.M2T(x), , j=j, value=value))

setReplaceMethod("[", signature(x = "RsparseMatrix", i = "index", j = "index",
				value = "replValue"),
		 function (x, i, j, ..., value)
		 replTmat(.M2T(x), i=i, j=j, value=value))

setReplaceMethod("[", signature(x = "RsparseMatrix", i = "index", j = "missing",
				value = "sparseVector"),
		 function (x, i, j, ..., value) {
                     if(nargs() == 3L)
                         replTmat(.M2T(x), i=i, value=value) # x[i] <- v
                     else replTmat(.M2T(x), i=i, , value=value) # x[i, ] <- v
                 })

setReplaceMethod("[", signature(x = "RsparseMatrix", i = "missing", j = "index",
				value = "sparseVector"),
		 function (x, i, j, ..., value)# extra " , ": want nargs() == 4
		 replTmat(.M2T(x), , j=j, value=value))

setReplaceMethod("[", signature(x = "RsparseMatrix", i = "index", j = "index",
				value = "sparseVector"),
		 function (x, i, j, ..., value)
		 replTmat(.M2T(x), i=i, j=j, value=value))


setReplaceMethod("[", signature(x = "RsparseMatrix", i = "matrix", j = "missing",
				value = "replValue"),
                 function (x, i, j, ..., value) {
                     if(nargs() == 3L)
                         .TM.repl.i.mat(.M2T(x), i=i, value=value)
                     else replTmat(.M2T(x), i=as.vector(i), , value=value)
                 })


## ==== TsparseMatrix ==================================================

##' a simplified "subset" of  intI() below
int2i <- function(i, n) {
    if(any(i < 0L)) {
	if(any(i > 0L))
	    stop("you cannot mix negative and positive indices")
	seq_len(n)[i]
    } else {
	if(length(i) && max(i, na.rm=TRUE) > n)
	    stop(gettextf("index larger than maximal %d", n), domain=NA)
	if(any(z <- i == 0)) i <- i[!z]
	i
    }
}

intI <- function(i, n, dn, give.dn = TRUE) {
    ## Purpose: translate numeric | logical | character index
    ##		into 0-based integer
    ## ----------------------------------------------------------------------
    ## Arguments: i: index vector (numeric | logical | character)
    ##		  n: array extent		    { ==  dim(.) [margin] }
    ##		 dn: character col/rownames or NULL { == dimnames(.)[[margin]] }
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 23 Apr 2007

    has.dn <- !is.null.DN(dn)
    DN <- has.dn && give.dn
    if(is.numeric(i) || is(i, "numeric")) { # inherits(<integer>, "numeric") is FALSE
	storage.mode(i) <- "integer"
	if(anyNA(i)) stop("'NA' indices are not (yet?) supported for sparse Matrices")
	if(any(i < 0L)) {
	    if(any(i > 0L))
		stop("you cannot mix negative and positive indices")
	    i0 <- (0:(n - 1L))[i]
	} else {
	    if(length(i) && max(i, na.rm=TRUE) > n) # base has "subscript out of bounds":
		stop(gettextf("index larger than maximal %d", n), domain=NA)
	    if(any(z <- i == 0)) i <- i[!z]
	    i0 <- i - 1L		# transform to 0-indexing
	}
	if(DN) dn <- dn[i]
    }
    else if (is.logical(i) || inherits(i, "logical")) {
	if(length(i) > n)
	    stop(gettextf("logical subscript too long (%d, should be %d)",
			  length(i), n), domain=NA)
	if(anyNA(i)) stop("'NA' indices are not (yet?) supported for sparse Matrices")
	i0 <- (0:(n - 1L))[i]
	if(DN) dn <- dn[i]
    } else { ## character
	if(!has.dn)
	    stop("no 'dimnames[[.]]': cannot use character indexing")
	i0 <- match(i, dn)
	if(anyNA(i0)) stop("invalid character indexing")
	if(DN) dn <- dn[i0]
	i0 <- i0 - 1L
    }
    if(!give.dn) i0 else list(i0 = i0, dn = dn)
} ## {intI}

.ind.prep <- function(xi, intIlist, iDup = duplicated(i0), anyDup = any(iDup)) {
    ## Purpose: do the ``common things'' for "*gTMatrix" indexing for 1 dim.
    ##		and return match(.,.) + li = length of corresponding dimension
    ##
    ## xi = "x@i" ; intIlist = intI(i, dim(x)[margin], ....)

    i0 <- intIlist$i0
    stopifnot(is.numeric(i0))# cheap fast check (i0 may have length 0 !)

    m <- match(xi, i0, nomatch=0)
    if(anyDup) { # assuming   anyDup <- any(iDup <- duplicated(i0))
	## i0i: where in (non-duplicated) i0 are the duplicated ones
	i0i <- match(i0[iDup], i0)
	i.x <- which(iDup) - 1L
	jm <- lapply(i0i, function(.) which(. == m))
    }

    c(list(m = m, li = length(i0),
	   i0 = i0, anyDup = anyDup, dn = intIlist$dn),
      ## actually,  iDup  is rarely needed in calling code
      if(anyDup) list(iDup = iDup, i0i = i0i, i.x = i.x,
		      jm = unlist(jm), i.xtra = rep.int(i.x, lengths(jm))))
} ## {.ind.prep}

##' <description>
##' Do the ``common things'' for "*gTMatrix" sub-assignment
##' for 1 dimension, 'margin' ,
##' <details>
##' @title Indexing Preparation
##' @param i "index"
##' @param margin in {1,2};
##' @param di = dim(x)	{ used when i is not character }
##' @param dn = dimnames(x)
##' @return match(.,.) + li = length of corresponding dimension
##' difference to .ind.prep(): use 1-indices; no match(xi,..), no dn at end
##' @author Martin Maechler
.ind.prep2 <- function(i, margin, di, dn)
    intI(i, n = di[margin], dn = dn[[margin]], give.dn = FALSE)

### FIXME: make this `very fast'  for the very very common case of
### -----   M[i,j] <- v  with   i,j = length-1-numeric;  v= length-1 number
###                            *and* M[i,j] == 0 previously
##
## FIXME(2): keep in sync with replCmat() in ./Csparse.R
## FIXME(3): It's terribly slow when used e.g. from diag(M[,-1]) <- value
## -----     which has "workhorse"   M[,-1] <- <dsparseVector>
##
## workhorse for "[<-" :
replTmat <- function (x, i, j, ..., value) {
## NOTE:  need '...', i.e., exact signature such that setMethod()
##	  does not use .local() such that nargs() will work correctly:
    di <- dim(x)
    dn <- dimnames(x)
    iMi <- missing(i)
    jMi <- missing(j)
    ## "FIXME": could pass this (and much ? more) when this function would not *be* a
    ## method but be *called* from methods

    clDv <- getClassDef(class(value))
    spV <- extends(clDv, "sparseVector")
    ## own version of all0() that works both for sparseVector and atomic vectors:
    .all0 <- function(v) if(spV) length(v@i) == 0 else all0(v)
    delayedAssign("value.not.logical",
                  !(if(spV) {
                      extends1of(clDv, "lsparseVector", "nsparseVector")
                  } else {
                      is.logical(value) || is.logical(as.vector(value))
                  }))
    na <- nargs()
    if(na == 3) { ## i = vector indexing  M[i] <- v,  e.g.,  M[TRUE] <- v or M[] <- v !
	Matrix.message("diagnosing replTmat(x,i,j,v): nargs()= 3; ",
		   if(iMi | jMi) sprintf("missing (i,j) = (%d,%d)", iMi,jMi))
	if(iMi) stop("internal bug: missing 'i' in replTmat(): please report")
	if(is.character(i))
	    stop("[ <character> ] indexing not allowed: forgot a \",\" ?")
	if(is.matrix(i))
	    stop("internal bug: matrix 'i' in replTmat(): please report")
	## Now: have  M[i] <- v	 with vector logical or "integer" i :
	## Tmatrix maybe non-unique, have an entry split into a sum of several ones:

	if(!is(x,"generalMatrix")) {
	    cl <- class(x)
	    x <- .M2gen(x)
	    Matrix.message("'sub-optimal sparse 'x[i] <- v' assignment: Coercing class ",
		       cl," to ",class(x))
	}
	nr <- di[1]
    x <- aggregateT(x)
	x.i <- .Call(m_encodeInd2, x@i, x@j, di=di, FALSE, FALSE)

        n <- prod(di)
	i <- if(is.logical(i)) { # full-size logical indexing
	    if(n) {
                if(isTRUE(i)) # shortcut
                    0:(n-1)
                else {
                    if(length(i) < n) i <- rep_len(i, n)
                    (0:(n-1))[i] # -> 0-based index vector as well {maybe LARGE!}
                }
	    } else integer(0)
	} else {
	    ## also works with *negative* indices etc:
	    int2i(as.integer(i), n) - 1L ## 0-based indices [to match m_encodeInd2()]
	}

        clx <- class(x)
        clDx <- getClassDef(clx) # extends(), is() etc all use the class definition
        has.x <- "x" %in% slotNames(clDx) # === slotNames(x)
	if(!has.x && # <==> "n.TMatrix"
	   ((iNA <- any(ina <- is.na(value))) || value.not.logical)) {
            if(value.not.logical) value <- as.logical(value)
	    if(iNA) {
		value[ina] <- TRUE
		warning(
		    gettextf("x[.] <- val: x is %s, val not in {TRUE, FALSE} is coerced; NA |--> TRUE.",
			     dQuote(clx)), domain=NA)
	    }
	    else warning(
		    gettextf("x[.] <- val: x is %s, val not in {TRUE, FALSE} is coerced.",
			     dQuote(clx)), domain=NA)
	}

	## now have 0-based indices   x.i (entries) and	 i (new entries)

	## the simplest case:
	if(.all0(value)) { ## just drop the non-zero entries
	    if(!all(sel <- is.na(match(x.i, i)))) { ## non-zero there
		x@i <- x@i[sel]
		x@j <- x@j[sel]
		if(has.x)
		    x@x <- x@x[sel]
		if(.hasSlot(x, "factors") && length(x@factors)) # drop cached ones
		    x@factors <- list()
	    }
	    return(x)
	}

	m <- length(i)
	if(length(value) != m) { ## use recycling rules
	    if(m %% length(value) != 0)
		warning("number of items to replace is not a multiple of replacement length")
	    value <- rep_len(value, m)
	}

        ## With duplicated entries i, only use the last ones!
        if(id <- anyDuplicated(i, fromLast=TRUE)) {
            i <- i[-id]
            value <- value[-id]
            if(any(id <- duplicated(i, fromLast=TRUE))) {
                nd <- -which(id)
                i <- i[nd]
                value <- value[nd]
            }
        }

	## matching existing non-zeros and new entries; isE := "is Existing"
	##  isE <- i %in% x.i;  mi <- {matching i's}
        isE <- !is.na(mi <- match(i, x.i))
        ## => mi[isE] entries in (i,j,x) to be set to new value[]s

	## 1) Change the matching non-zero entries
	if(has.x)
	    x@x[mi[isE]] <- as(value[isE], class(x@x))
        else if(any0(value[isE])) { ## "n.TMatrix" : remove (i,j) where value is FALSE
            get0 <- !value[isE] ## x[i,j] is TRUE, should become FALSE
            i.rm <- - mi[isE][get0]
            x@i <- x@i[i.rm]
            x@j <- x@j[i.rm]
        }
	## 2) add the new non-zero entries
	i <- i[!isE]
	xv <- value[!isE]
	## --- Be be efficient when  'value' is sparse :
	if(length(notE <- which(isN0(xv)))) { # isN0(): non-0's; NAs counted too
	    xv <- xv[notE]
	    i <- i[notE]
	    if(has.x) {
		x@x <- c(x@x, as(xv, class(x@x)))
	    } else { # n.TMatrix : assign (i,j) only where value is TRUE:
		i <- i[xv]
	    }
	    x@i <- c(x@i, i %%  nr)
	    x@j <- c(x@j, i %/% nr)
	}
	if(.hasSlot(x, "factors") && length(x@factors)) # drop cached ones
	    x@factors <- list()
	return(x)
    } ## {nargs = 3;  x[ii] <- value }

    ## nargs() == 4 :  x[i,j] <- value
    ## --------------------------------------------------------------------------
    lenV <- length(value)
    Matrix.message(".. replTmat(x,i,j,v): nargs()= 4; cl.(x)=",
	       class(x),"; len.(value)=", lenV,"; ",
	       if(iMi | jMi) sprintf("missing (i,j) = (%d,%d)", iMi,jMi),
	       .M.level = 2)# level 1  gives too many messages

    ## FIXME: use  'abIndex' or a better algorithm, e.g.  if(iMi)
    i1 <- if(iMi) 0:(di[1] - 1L) else .ind.prep2(i, 1, di, dn)
    i2 <- if(jMi) 0:(di[2] - 1L) else .ind.prep2(j, 2, di, dn)
    dind <- c(length(i1), length(i2)) # dimension of replacement region
    lenRepl <- prod(dind)
    if(lenV == 0) {
        if(lenRepl != 0)
            stop("nothing to replace with")
        else return(x)
    }
    ## else: lenV := length(value)	 is > 0
    if(lenRepl %% lenV != 0)
	stop("number of items to replace is not a multiple of replacement length")
    if(!spV && lenRepl > 2^16) { # (somewhat arbitrary cutoff)
	value <- as(value, "sparseVector")# so that subsequent rep(.) are fast
        spV <- TRUE
    }
    ## Now deal with duplicated / repeated indices: "last one wins"
    if(!iMi && any(dup <- duplicated(i1, fromLast = TRUE))) { ## duplicated rows
        keep <- !dup
        i1 <- i1[keep]
        ## keep is "internally" recycled below {and that's important: it is dense!}
	lenV <- length(value <- rep_len(value, lenRepl)[keep])
        dind[1] <- length(i1)
        lenRepl <- prod(dind)
    }
    if(!jMi && any(dup <- duplicated(i2, fromLast = TRUE))) { ## duplicated columns
        iDup <- which(dup)
        ## The following is correct, but  rep(keep,..) can be *HUGE*
        ## keep <- !dup
        ## i2 <- i2[keep]
	## lenV <- length(value <- rep_len(value, lenRepl)[rep(keep, each=dind[1])])
        ## solution: sv[-i] is efficient for sparseVector:
        i2 <- i2[- iDup]
        nr <- dind[1]
        iDup <- rep((iDup - 1)*nr, each=nr) + seq_len(nr)
	lenV <- length(value <- rep_len(value, lenRepl)[-iDup])
        dind[2] <- length(i2)
        lenRepl <- prod(dind)
    }
    clx <- class(x)
    clDx <- getClassDef(clx) # extends() , is() etc all use the class definition
    stopifnot(extends(clDx, "TsparseMatrix"))
    ## Tmatrix maybe non-unique, have an entry split into a sum of several ones:
    x <- aggregateT(x)

    toGeneral <- r.sym <- FALSE
    if(extends(clDx, "symmetricMatrix")) {
	## using array() for large dind is a disaster...
	mkArray <- if(spV) # TODO: room for improvement
	    function(v, dim) spV2M(v, dim[1],dim[2]) else array
	r.sym <-
	    (dind[1] == dind[2] && all(i1 == i2) &&
	     (lenRepl == 1 || lenV == 1 ||
	      isSymmetric(mkArray(value, dim=dind))))
	if(r.sym) { ## result is *still* symmetric --> keep symmetry!
	    xU <- x@uplo == "U"
            # later, we will consider only those indices above / below diagonal:
	}
	else toGeneral <- TRUE
    } else if(extends(clDx, "triangularMatrix")) {
        xU <- x@uplo == "U"
	r.tri <- ((any(dind == 1) || dind[1] == dind[2]) &&
		  if(xU) max(i1) <= min(i2) else max(i2) <= min(i1))
	if(r.tri) { ## result is *still* triangular
            if(any(i1 == i2)) # diagonal will be changed
                x <- diagU2N(x) # keeps class (!)
	}
	else toGeneral <- TRUE
    }
    if(toGeneral) { # go to "generalMatrix" and continue
        Matrix.message("M[i,j] <- v :  coercing symmetric M[] into non-symmetric")
        x <- .M2gen(x)
        clDx <- getClassDef(clx <- class(x))
    }

    ## TODO (efficiency): replace  'sel' by 'which(sel)'
    get.ind.sel <- function(ii,ij)
	(match(x@i, ii, nomatch = 0L) & match(x@j, ij, nomatch = 0L))
    ## sel[k] := TRUE iff k-th non-zero entry (typically x@x[k]) is to be replaced
    sel <- get.ind.sel(i1,i2)

    has.x <- "x" %in% slotNames(clDx) # === slotNames(x)

    ## the simplest case: for all Tsparse, even for i or j missing
    if(.all0(value)) { ## just drop the non-zero entries
	if(any(sel)) { ## non-zero there
	    x@i <- x@i[!sel]
	    x@j <- x@j[!sel]
            if(has.x)
		x@x <- x@x[!sel]
	    if(.hasSlot(x, "factors") && length(x@factors)) # drop cached ones
		x@factors <- list()
	}
	return(x)
    }
    ## else --  some( value != 0 ) --
    if(lenV > lenRepl)
        stop("too many replacement values")
    ## now have  lenV <= lenRepl

    if(!has.x && # <==> "n.TMatrix"
       ((iNA <- anyNA(value)) || value.not.logical))
	warning(if(iNA)
		gettextf("x[.,.] <- val: x is %s, val not in {TRUE, FALSE} is coerced NA |--> TRUE.",
			 dQuote(clx))
		else
		gettextf("x[.,.] <- val: x is %s, val not in {TRUE, FALSE} is coerced.",
			 dQuote(clx)), domain=NA)

    ## another simple, typical case:
    if(lenRepl == 1) {
        if(spV && has.x) value <- as(value, "vector")
        if(any(sel)) { ## non-zero there
            if(has.x)
                x@x[sel] <- value
        } else { ## new non-zero
            x@i <- c(x@i, i1)
            x@j <- c(x@j, i2)
            if(has.x)
                x@x <- c(x@x, value)
        }
	if(.hasSlot(x, "factors") && length(x@factors)) # drop cached ones
	    x@factors <- list()
        return(x)
    }

### Otherwise, for large lenRepl, we get into trouble below

    if(lenRepl > 2^20) { # (somewhat arbitrary cutoff)
## FIXME: just for testing !!
## if(identical(Sys.getenv("USER"),"maechler")
##    if(lenRepl > 2) { # __________ ___ JUST for testing! _______________
	if(!isTRUE(getOption("Matrix.quiet")))
	    message(gettextf("x[.,.] <- val : x being coerced from Tsparse* to CsparseMatrix"),
		    domain = NA)
	return(replCmat4(.M2C(x), i1, i2, iMi=iMi, jMi=jMi,
			 value = if(spV) value else as(value, "sparseVector"),
			 spV = TRUE))
    }

    ##     if(r.sym) # value already adjusted, see above
    ##        lenRepl <- length(value) # shorter (since only "triangle")
    if(!r.sym && lenV < lenRepl)
	value <- rep_len(value, lenRepl)

    ## now:  length(value) == lenRepl  {but value is sparseVector if it's "long" !}

    ## value[1:lenRepl]:  which are structural 0 now, which not?
    ## v0 <- is0(value)
    ## - replaced by using isN0(as.vector(.)) on a typical small subset value[.]
    ## --> more efficient for sparse 'value' & large 'lenRepl' :
    ## FIXME [= FIXME(3) above]:
    ## ----- The use of  seq_len(lenRepl) below is *still* inefficient
    ##   (or impossible e.g. when lenRepl == 50000^2)
    ##       and the  vN0 <- isN0(as.vector(value[iI0]))  is even more ...

    ## One idea: use "abIndex", (a very efficient storage of index vectors which are
    ## a concatenation of only a few arithmetic seq()ences
    use.abI <- isTRUE(getOption("Matrix.use.abIndex"))
    ## This 'use.abI' should later depend on the *dimension* of things !
    ##>>> But for that, we need to implement the following abIndex - "methods":
    ##>>>   <abI>[-n],  <value>[ <abIndex> ] , intersect(<abI>, <abI>)
    ## and for intersect(): typically sort(), unique() & similar

    iI0 <- if(use.abI) abIseq1(1L, lenRepl) else seq_len(lenRepl)

    if(any(sel)) {
	## the 0-based indices of non-zero entries -- WRT to submatrix
	iN0 <- 1L + .Call(m_encodeInd2,
			  match(x@i[sel], i1),
			  match(x@j[sel], i2),
			  di = dind, orig1=TRUE, FALSE)

	## 1a) replace those that are already non-zero with non-0 values
	vN0 <- isN0(value[iN0])
	if(any(vN0) && has.x) {
	    vv0 <- which(vN0)
	    x@x[sel][vv0] <- as.vector(value[iN0[vv0]])
	}

	## 1b) replace non-zeros with 0 --> drop entries
	if(!all(vN0)) { ##-> ii will not be empty
	    ii <- which(sel)[which(!vN0)] # <- vN0 may be sparseVector
	    if(has.x)
		x@x <- x@x[-ii]
	    x@i <- x@i[-ii]
	    x@j <- x@j[-ii]
	}
	iI0 <- if(length(iN0) < lenRepl) iI0[-iN0] ## else NULL
                                        # == complementInd(non0, dind)
    }
    if(length(iI0)) {
        if(r.sym) {
	    ## should only set new entries above / below diagonal, i.e.,
            ## subset iI0 such as to contain only  above/below ..
	    iSel <-
		if(use.abI) abIindTri(dind[1], upper=xU, diag=TRUE)
		else	       indTri(dind[1], upper=xU, diag=TRUE)
	    ## select also the corresponding triangle of values
### TODO for "abIndex" -- note we KNOW that both  iI0 and iSel
### are strictly increasing :
	    iI0 <- intersect(iI0, iSel)
        }
        full <- length(iI0) == lenRepl
	vN0 <-
	    if(spV) ## "sparseVector"
		(if(full) value else value[iI0])@i
	    else which(isN0(if(full) value else value[iI0]))
	if(length(vN0)) {
	    ## 2) add those that were structural 0 (where value != 0)
	    iIN0 <- if(full) vN0 else iI0[vN0]
	    ij0 <- decodeInd(iIN0 - 1L, nr = dind[1])
	    x@i <- c(x@i, i1[ij0[,1] + 1L])
	    x@j <- c(x@j, i2[ij0[,2] + 1L])
	    if(has.x)
		x@x <- c(x@x, as.vector(value[iIN0]))
	}
    }
    if(.hasSlot(x, "factors") && length(x@factors)) # drop cached ones
	x@factors <- list()
    x
} ## end{replTmat}

## A[ ij ] <- value,  where ij is a matrix; typically (i,j) 2-column matrix :
## ----------------   ./Matrix.R has a general cheap method
## This one should become as fast as possible -- is also used from Csparse.R --
.TM.repl.i.mat <- function (x, i, j, ..., value) {
    nA <- nargs()
    if(nA != 3)
	stop(gettextf("nargs() = %d should never happen; please report.", nA), domain=NA)

    ## else: nA == 3  i.e.,  M [ cbind(ii,jj) ] <- value or M [ Lmat ] <- value
    if(is.logical(i)) {
	Matrix.message(".TM.repl.i.mat(): drop 'matrix' case ...", .M.level=2)
	## c(i) : drop "matrix" to logical vector
	x[as.vector(i)] <- value
	return(x)
    } else if(extends1of(cli <- getClassDef(class(i)), c("lMatrix", "nMatrix"))) {
	Matrix.message(".TM.repl.i.mat(): \"lMatrix\" case ...", .M.level=2)
	i <- which(as(i, if(extends(cli, "sparseMatrix")) "sparseVector" else "vector"))
	## x[i] <- value ; return(x)
	return(`[<-`(x,i, value=value))
    } else if(extends(cli, "Matrix")) { # "dMatrix" or "iMatrix"
	if(ncol(i) != 2)
	    stop("such indexing must be by logical or 2-column numeric matrix")
	i <- as(i, "matrix")
    } else if(!is.numeric(i) || ncol(i) != 2)
	stop("such indexing must be by logical or 2-column numeric matrix")
    if(!is.integer(i)) storage.mode(i) <- "integer"
    if(any(i < 0))
	stop("negative values are not allowed in a matrix subscript")
    if(anyNA(i))
	stop("NAs are not allowed in subscripted assignments")
    if(any(i0 <- (i == 0))) # remove them
	i <- i[ - which(i0, arr.ind = TRUE)[,"row"], ]
    if(length(attributes(i)) > 1) # more than just 'dim'; simplify: will use identical
	attributes(i) <- list(dim = dim(i))
    ## now have integer i >= 1
    m <- nrow(i)
    if(m == 0)
	return(x)
    if(length(value) == 0)
	stop("nothing to replace with")
    ## mod.x <- .type.kind[.M.kind(x)]
    if(length(value) != m) { ## use recycling rules
	if(m %% length(value) != 0)
	    warning("number of items to replace is not a multiple of replacement length")
	value <- rep_len(value, m)
    }
    clx <- class(x)
    clDx <- getClassDef(clx) # extends() , is() etc all use the class definition
    stopifnot(extends(clDx, "TsparseMatrix"))

    di <- dim(x)
    nr <- di[1]
    nc <- di[2]
    i1 <- i[,1]
    i2 <- i[,2]
    if(any(i1 > nr)) stop(gettextf("row indices must be <= nrow(.) which is %d", nr), domain=NA)
    if(any(i2 > nc)) stop(gettextf("column indices must be <= ncol(.) which is %d", nc), domain=NA)

    ## Tmatrix maybe non-unique, have an entry split into a sum of several ones:
    x <- aggregateT(x)

    toGeneral <- FALSE
    isN <- extends(clDx, "nMatrix")
    if(r.sym <- extends(clDx, "symmetricMatrix")) {
	## Tests to see if the assignments are symmetric as well
	r.sym <- all(i1 == i2)
	if(!r.sym) { # do have *some* Lower or Upper entries
	    iL <- i1 > i2
	    iU <- i1 < i2
	    r.sym <- sum(iL) == sum(iU) # same number
	    if(r.sym) {
		iLord <- order(i1[iL], i2[iL])
		iUord <- order(i2[iU], i1[iU]) # row <-> col. !
		r.sym <- {
		    identical(i[iL,    , drop=FALSE][iLord,],
			      i[iU, 2:1, drop=FALSE][iUord,]) &&
		    all(value[iL][iLord] ==
			value[iU][iUord])
		}
	    }
	}
	if(r.sym) { ## result is *still* symmetric --> keep symmetry!
	    ## now consider only those indices above / below diagonal:
	    useI <- if(x@uplo == "U") i1 <= i2 else i2 <= i1
	    i <- i[useI, , drop=FALSE]
	    value <- value[useI]
	}
	else toGeneral <- TRUE
    }
    else if(extends(clDx, "triangularMatrix")) {
	r.tri <- all(if(x@uplo == "U") i1 <= i2 else i2 <= i1)
	if(r.tri) { ## result is *still* triangular
	    if(any(ieq <- i1 == i2)) { # diagonal will be changed
		if(x@diag == "U" && all(ieq) &&
		   all(value == if(isN) TRUE else as1(x@x)))
		    ## only diagonal values are set to 1 -- i.e. unchanged
		    return(x)
		x <- diagU2N(x) # keeps class (!)
	    }
	}
	else toGeneral <- TRUE
    }
    if(toGeneral) { # go to "generalMatrix" and continue
	Matrix.message("M[ij] <- v :  coercing symmetric M[] into non-symmetric")
	x <- .M2gen(x)
	clDx <- getClassDef(clx <- class(x))
    }

    ii.v <- .Call(m_encodeInd, i, di, orig1=TRUE, checkBounds = TRUE)
    if(id <- anyDuplicated(ii.v, fromLast=TRUE)) {
        Matrix.message("M[ij] <- v :  duplicate ij-entries; using last")
        ii.v  <- ii.v [-id]
	value <- value[-id]
        if(any(id <- duplicated(ii.v, fromLast=TRUE))) {
            nd <- -which(id)
            ii.v  <- ii.v [nd]
            value <- value[nd]
        }
    }
    ii.x <- .Call(m_encodeInd2, x@i, x@j, di, FALSE, FALSE)
    m1 <- match(ii.v, ii.x)
    i.repl <- !is.na(m1) # those that need to be *replaced*

    if(isN) { ## no 'x' slot
	isN <- is.logical(value) # will result remain  "nMatrix" ?
	if(!isN)
            x <- .M2kind(x, "d")
    }
    has.x <- !isN ## isN  <===> "remains pattern matrix" <===> has no 'x' slot

    if(any(i.repl)) { ## some to replace at matching (@i, @j)
	if(has.x)
	    x@x[m1[i.repl]] <- value[i.repl]
	else { # nMatrix ; eliminate entries that are set to FALSE; keep others
	    if(any(isF <- is0(value[i.repl])))  {
		ii <- m1[i.repl][isF]
		x@i <- x@i[ -ii]
		x@j <- x@j[ -ii]
	    }
	}
    }
    if(any(i.new <- !i.repl & isN0(value))) { ## some new entries
	i.j <- decodeInd(ii.v[i.new], nr)
	x@i <- c(x@i, i.j[,1])
	x@j <- c(x@j, i.j[,2])
	if(has.x)
	    x@x <- c(x@x, value[i.new])
    }

    if(.hasSlot(x, "factors") && length(x@factors)) # drop cached ones
	x@factors <- list()
    x
} ## end{.TM.repl.i.mat}

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "index", j = "missing",
				value = "replValue"),
		 replTmat)

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "missing", j = "index",
				value = "replValue"),
		 replTmat)

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "index", j = "index",
				value = "replValue"),
		 replTmat)

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "matrix", j = "missing",
				value = "replValue"),
		 .TM.repl.i.mat)
setReplaceMethod("[", signature(x = "TsparseMatrix", i = "Matrix", j = "missing",
				value = "replValue"),
		 .TM.repl.i.mat)


### When the RHS 'value' is  a sparseVector, now can use  replTmat  as well
setReplaceMethod("[", signature(x = "TsparseMatrix", i = "missing", j = "index",
				value = "sparseVector"),
		 replTmat)

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "index", j = "missing",
				value = "sparseVector"),
		 replTmat)

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "index", j = "index",
				value = "sparseVector"),
		 replTmat)


## ==== diagonalMatrix =================================================

## When you assign to a diagonalMatrix, the result should be
## diagonal or sparse ---
replDiag <- function(x, i, j, ..., value) {
## FIXME: if   (i == j)  &&  isSymmetric(value) then -- want symmetricMatrix result! -- or diagMatrix
    x <- .diag2sparse(x, ".", "g", "C") # was ->TsparseMatrix till 2012-07
    if(missing(i))
        x[, j] <- value
    else if(missing(j)) { ##  x[i , ] <- v  *OR*   x[i] <- v
        na <- nargs()
        ##         message("diagnosing replDiag() -- nargs()= ", na)
        if(na == 4L)
            x[i, ] <- value
        else if(na == 3L)
            x[i] <- value
        else stop(gettextf("Internal bug: nargs()=%d; please report",
                           na), domain=NA)
    } else
        x[i,j] <- value
    ## TODO: the following is a bit expensive; have cases above e.g. [i,] where
    ## ----- we could check *much* faster :
    if(isDiagonal(x))
        forceDiagonal(x)
    else if(isSymmetric(x))
        forceSymmetric(x)
    else if(!(it <- isTriangular(x)))
        x
    else if(attr(it, "kind") == "U")
        triu(x)
    else tril(x)
}

setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index",
                                j = "index", value = "replValue"), replDiag)

setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index",
                                j = "missing", value = "replValue"),
                 function(x,i,j, ..., value) {
                     ## message("before replDiag() -- nargs()= ", nargs())
                     if(nargs() == 3L)
                         replDiag(x, i=i, value=value)
                     else ## nargs() == 4 :
                         replDiag(x, i=i, , value=value)
                 })

setReplaceMethod("[", signature(x = "diagonalMatrix", i = "missing",
                                j = "index", value = "replValue"),
                 function(x,i,j, ..., value) replDiag(x, j=j, value=value))

## x[] <- value :
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "missing",
                                j = "missing", value = "ANY"),
                 function(x,i,j, ..., value) {
                     if(all0(value)) { # be faster
                         r <- new(paste0(.M.kind(x), "tTMatrix")) # of all "0"
                         r@Dim <- x@Dim
                         r@Dimnames <- x@Dimnames
                         r
                     } else {
                         ## typically non-sense: assigning to full sparseMatrix
                         x[TRUE] <- value
                         x
                     }
                 })


setReplaceMethod("[", signature(x = "diagonalMatrix",
                                i = "matrix", # 2-col.matrix
                                j = "missing", value = "replValue"),
                 function(x,i,j, ..., value) {
                     if(ncol(i) == 2L) {
                         if(all((ii <- i[,1L]) == i[,2L])) {
                             ## replace in diagonal only
                             if(x@diag == "U") {
                                 one <- as1(x@x)
                                 if(any(value != one | is.na(value))) {
                                     x@diag <- "N"
                                     x@x <- rep.int(one, x@Dim[1L])
                                 } else return(x)
                             }
                             x@x[ii] <- value
                             x
                         } else { ## no longer diagonal, but remain sparse:
### FIXME:  use  uplo="U" or uplo="L"  (or *not* "triangularMatrix")
### depending on LE <- i <= j
### all(LE) //  all(!LE) // remaining cases
                             x <- .diag2sparse(x, ".", "t", "C") # was ->TsparseMatrix
                             x[i] <- value
                             x
                         }
                     }
                     else { # behave as "base R": use as if vector
                         x <- as(x, "matrix")
                         x[i] <- value
                         Matrix(x)
                     }
                 })


## value = "sparseMatrix":
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "missing", j = "index",
                                value = "sparseMatrix"),
                 function (x, i, j, ..., value)
                     callGeneric(x=x, , j=j, value=as(value, "sparseVector")))

setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index", j = "missing",
                                value = "sparseMatrix"),
                 function (x, i, j, ..., value)
                     callGeneric(x=x, i=i, , value=as(value, "sparseVector")))
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index", j = "index",
                                value = "sparseMatrix"),
                 function (x, i, j, ..., value)
                     callGeneric(x=x, i=i, j=j, value=as(value, "sparseVector")))

## value = "sparseVector":
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "missing", j = "index",
                                value = "sparseVector"),
                 replDiag)
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index", j = "missing",
                                value = "sparseVector"),
                 replDiag)
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index", j = "index",
                                value = "sparseVector"),
                 replDiag)


## ==== indMatrix ======================================================

.indMatrix.sub <- function(x, i, j, ..., value) {
    x <- as(x, "TsparseMatrix")
    callGeneric()
}
for (.i in c("missing", "index"))
for (.j in c("missing", "index"))
setReplaceMethod("[", signature(x = "indMatrix", i = .i, j = .j, value = "ANY"),
                 .indMatrix.sub)
rm(.indMatrix.sub, .i, .j)


## ==== sparseVector ===================================================

## This is a simplified intI() -- for sparseVector indexing:
intIv <- function(i, n, cl.i = getClassDef(class(i))) {
### Note: undesirable to use this for negative indices;
### ----  using seq_len(n) below means we are  NON-sparse ...
### Fixed, for "x[i] with negative i" at least.

    ## Purpose: translate numeric | logical index     into  1-based integer
    ## --------------------------------------------------------------------
    ## Arguments: i: index vector (numeric | logical) *OR* sparseVector
    ##		  n: array extent { ==	length(.) }
    if(missing(i))
	seq_len(n)
    else if(extends(cl.i, "numeric")) {
        ## not ok, when max(i) > .Machine$integer.max !  storage.mode(i) <- "integer"
        int2i(i,n) ##-> ./Tsparse.R
    }
    else if (extends(cl.i, "logical")) {
	seq_len(n)[i]
    } else if(extends(cl.i, "nsparseVector")) {
	i@i # the indices are already there !
    } else if(extends(cl.i, "lsparseVector")) {
	i@i[i@x] # "drop0", i.e. FALSE; NAs ok
    } else if (extends(cl.i, "sparseVector")) { ## 'i'sparse, 'd'sparse	 (etc)
	as.integer(i@x[i@i])
    }
    else
        stop("index must be numeric, logical or sparseVector for indexing sparseVectors")
} ## intIv()

replSPvec <- function (x, i, value) {
    n <- x@length
    ii <- intIv(i, n)
    lenRepl <- length(ii)
    if(!lenRepl) return(x)
    ## else:  lenRepl = length(ii) > 0
    lenV <- length(value)
    if(lenV == 0)
        stop("nothing to replace with")
    ## else: lenV := length(value) > 0
    if(lenRepl %% lenV != 0)
	stop("number of items to replace is not a multiple of replacement length")
    if(anyDuplicated(ii)) { ## multiple *replacement* indices: last one wins
	## TODO: in R 2.6.0 use	 duplicate(*, fromLast=TRUE)
	ir <- lenRepl:1
	keep <- match(ii, ii[ir]) == ir
	ii <- ii[keep]
	lenV <- length(value <- rep(value, length.out = lenRepl)[keep])
	lenRepl <- length(ii)
    }

    has.x <- .hasSlot(x, "x")## has "x" slot
    m <- match(x@i, ii, nomatch = 0)
    sel <- m > 0L

    ## the simplest case
    if(all0(value)) { ## just drop the non-zero entries
	if(any(sel)) { ## non-zero there
	    x@i <- x@i[!sel]
	    if(has.x)
		x@x <- x@x[!sel]
	}
	return(x)
    }
    ## else --	some( value != 0 ) --
    if(lenV > lenRepl)
	stop("too many replacement values")
    else if(lenV < lenRepl)
	value <- rep(value, length.out = lenRepl)
    ## now:  length(value) == lenRepl > 0

    v0 <- is0(value)
    ## value[1:lenRepl]:  which are structural 0 now, which not?
    v.sp <- inherits(value, "sparseVector")

    if(any(sel)) {
	## indices of non-zero entries -- WRT to subvector
	iN0 <- m[sel] ## == match(x@i[sel], ii)

	## 1a) replace those that are already non-zero with new val.
	vN0 <- !v0[iN0]
	if(any(vN0) && has.x) {
	    vs <- value[iN0[vN0]]
	    x@x[sel][vN0] <- if(v.sp) sp2vec(vs, mode=typeof(x@x)) else vs
	}
	## 1b) replace non-zeros with 0 --> drop entries
	if(any(!vN0)) {
	    i <- which(sel)[!vN0]
	    if(has.x)
		x@x <- x@x[-i]
	    x@i <- x@i[-i]
	}
	iI0 <- if(length(iN0) < lenRepl) seq_len(lenRepl)[-iN0] # else NULL
    } else iI0 <- seq_len(lenRepl)

    if(length(iI0) && any(vN0 <- !v0[iI0])) {
	## 2) add those that were structural 0 (where value != 0)
	ij0 <- iI0[vN0]
	ii <- c(x@i, ii[ij0]) # new x@i, must be sorted:
	iInc <- sort.list(ii)
	x@i <- ii[iInc]
	if(has.x) # new @x, sorted along '@i':
	    x@x <- c(x@x, if(v.sp)
			      sp2vec(value[ij0], mode=typeof(x@x))
			  else value[ij0]
		     )[iInc]
    }
    x
}

setReplaceMethod("[", signature(x = "sparseVector", i = "index", j = "missing",
				value = "replValueSp"),
		 replSPvec)

setReplaceMethod("[", signature(x = "sparseVector",
                                i = "sparseVector", j = "missing",
				value = "replValueSp"),
                 ## BTW, the important case: 'i' a *logical* sparseVector
		 replSPvec)
