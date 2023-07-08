## METHODS FOR CLASS: CsparseMatrix (virtual)
## sparse matrices in compressed sparse column (CSC) format
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.validateCsparse <- function(x, sort.if.needed = FALSE)
    .Call(Csparse_validate2, x, sort.if.needed)
##-> to be used in sparseMatrix(.), e.g. --- but is unused currently
## NB: 'sort.if.needed' is called 'maybe_modify' in C -- so be careful
## more useful:
.sortCsparse <- function(x) .Call(Csparse_sort, x) ## modifies 'x' !!

## workhorse for "[<-" -- for d*, l*, and n..C-sparse matrices :
## ---------     -----
replCmat <- function (x, i, j, ..., value)
{
    di <- dim(x)
    dn <- dimnames(x)
    iMi <- missing(i)
    jMi <- missing(j)
    na <- nargs()
    Matrix.msg("replCmat[x,i,j,..,val] : nargs()=", na, "; ",
	       if(iMi || jMi) sprintf("missing (i,j) = (%d,%d)", iMi, jMi),
	       .M.level = 2)
    if(na == 3L) { ## vector (or 2-col) indexing M[i] <- v : includes M[TRUE] <- v or M[] <- v !
	x <- .CR2T(x)
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
                      spV = is(value, "sparseVector"))
{
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
	x <- .sparse2g(x) ## but do *not* redefine clx!
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
	    x <- .sparse2g(x) # was as(x, paste0(.M.kind(x), "gCMatrix"))
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
		   else .sparse2g(x), # must get "dgCMatrix"
		   i1, i2,
		   as(value, "sparseVector"))
    }
    else if(extends(clDx, "lMatrix")) {
	has.x <- TRUE
	x <- .Call(lCsparse_subassign,
		   if(clx %in% c("lgCMatrix", "ltCMatrix")) x
		   else .sparse2g(x), # must get "lgCMatrix"
		   i1, i2,
		   as(value, "sparseVector"))
    }
    else if(extends(clDx, "nMatrix")) {
	has.x <- FALSE
	x <- .Call(nCsparse_subassign,
		   if(clx %in% c("ngCMatrix", "ntCMatrix"))x
		   else .sparse2g(x), # must get "ngCMatrix"
		   i1, i2,
		   as(value, "sparseVector"))
    }
    else if(extends(clDx, "iMatrix")) {
	has.x <- TRUE
	x <- .Call(iCsparse_subassign,
		   if(clx %in% c("igCMatrix", "itCMatrix"))x
		   else .sparse2g(x), # must get "igCMatrix"
		   i1, i2,
		   as(value, "sparseVector"))
    }
    else if(extends(clDx, "zMatrix")) {
	has.x <- TRUE
	x <- .Call(zCsparse_subassign,
		   if(clx %in% c("zgCMatrix", "ztCMatrix"))x
		   else .sparse2g(x), # must get "zgCMatrix"
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
	    if(has0) x <- .Call(Csparse_drop, x, 0)

	    return(if(x.sym) as_CspClass(x, clx) else x)
	}
	## else go via Tsparse.. {FIXME: a waste! - we already have 'xj' ..}
	## and inside  Tsparse... the above i1, i2,..., sel  are *all* redone!
	## Happens too often {not anymore, I hope!}
	##
	Matrix.msg("wasteful C -> T -> C in replCmat(x,i,j,v) for <sparse>[i,j] <- v")
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

setMethod("writeMM", "CsparseMatrix",
	  function(obj, file, ...)
	  .Call(Csparse_MatrixMarket, obj, path.expand(as.character(file))))

dmperm <- function(x, nAns = 6L, seed = 0L) {
    stopifnot(length(nAns <- as.integer(nAns)) == 1L, nAns %in% c(2L, 4L, 6L),
              length(seed <- as.integer(seed)) == 1L, seed %in% -1:1)
    if(isS4(x)) {
        cld <- getClassDef(class(x))
        if(!extends(cld, "CsparseMatrix"))
            cld <- getClassDef(class(x <- as(x, "CsparseMatrix")))
        if(extends(cld, "symmetricMatrix"))
            cld <- getClassDef(class(x <- .sparse2g(x)))
        if(!(extends(cld, "dMatrix") || extends(cld, "nMatrix")))
            x <- ..sparse2d(x)
    } else { # typically a traditional matrix
        x <- .m2sparse(x, "dgC", NULL, NULL)
    }
    .Call(Csparse_dmperm, x, seed, nAns) # tolerating only [dn][gt]CMatrix 'x'
}
