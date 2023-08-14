## METHODS FOR CLASS: sparseMatrix (virtual)
## sparse matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~ CONSTRUCTORS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

spMatrix <- function(nrow, ncol, i = integer(), j = integer(), x = double())
    new(paste0(.M.kind(x), "gTMatrix"), # rely on new() to check validity
        Dim = c(as.integer(nrow), as.integer(ncol)),
        i = as.integer(i) - 1L,
        j = as.integer(j) - 1L,
        x = if(is.integer(x)) as.double(x) else x)

sparseMatrix <- function(i, j, p, x, dims, dimnames,
                         symmetric = FALSE, triangular = FALSE, index1 = TRUE,
                         repr = c("C", "R", "T"), giveCsparse,
                         check = TRUE, use.last.ij = FALSE)
{
    if((m.i <- missing(i)) + (m.j <- missing(j)) + (m.p <- missing(p)) != 1L)
        stop("exactly one of 'i', 'j', and 'p' must be missing from call")
    if(symmetric && triangular)
        stop("use Diagonal() to construct diagonal (symmetric && triangular) sparse matrices")
    index1 <- as.logical(index1) # allowing {0,1}

    repr <-
        ## NB: prior to 2020-05, we had 'giveCsparse' {T->"C" [default], F->"T"}
        ##     but no 'repr' ... the following is to remain backwards compatible
        if(missing(giveCsparse))
            match.arg(repr)
        else if(!missing(repr)) {
            warning("'giveCsparse' is deprecated; using 'repr' instead")
            match.arg(repr)
        ## } else {
        ##     repr <- if(giveCsparse) "C" else "T"
        ##     warning(gettextf("'giveCsparse' is deprecated; setting repr=\"%s\" for you", repr),
        ##             domain = NA)
        ## }
        } else if(giveCsparse) {
            ## NOT YET:
            ## warning("'giveCsparse' is deprecated; setting repr=\"C\" for you")
            "C"
        } else {
            warning("'giveCsparse' is deprecated; setting repr=\"T\" for you")
            "T"
        }

    if(!m.p) {
        p <- as.integer(p)
        if((n.p <- length(p)) == 0L || anyNA(p) || p[1L] != 0L ||
           any((dp <- p[-1L] - p[-n.p]) < 0L))
            stop("'p' must be a nondecreasing vector c(0, ...)")
        if((n.dp <- length(dp)) > .Machine$integer.max)
            stop("dimensions cannot exceed 2^31-1")
        i. <- rep.int(seq.int(from = 0L, length.out = n.dp), dp)
        if(m.i) i <- i. else j <- i.
    }

    if(!m.i)
        i <- if(index1) as.integer(i) - 1L else as.integer(i) # need 0-index
    if(!m.j)
        j <- if(index1) as.integer(j) - 1L else as.integer(j) # need 0-index

    rij <- cbind(if(n.i <- length(i)) range(i) else 0:-1,
                 if(n.j <- length(j)) range(j) else 0:-1,
                 deparse.level = 0L)
    if(anyNA(rij))
        stop("'i' and 'j' must not contain NA") # and not overflow
    if(any(rij[1L, ] < 0L))
        stop("'i' and 'j' must be ", if(index1) "positive" else "non-negative")
    dims <-
        if(!missing(dims)) {
            if(length(dims) != 2L ||
               any(is.na(dims) | dims < 0L | dims >= .Machine$integer.max + 1))
                stop("invalid 'dims'")
            if(any(dims - 1L < rij[2L, ]))
                stop("'dims' must contain all (i,j) pairs")
            as.integer(dims)
        } else if(symmetric || triangular)
            rep.int(max(rij), 2L) + 1L
        else rij[2L, ] + 1L

    kind <- if(m.x <- missing(x)) "n" else .M.kind(x)
    shape <-
        if(symmetric) {
            if(dims[1L] != dims[2L])
                stop("symmetric matrix must be square")
            "s"
        } else if(triangular) {
            if(dims[1L] != dims[2L])
                stop("triangular matrix must be square")
            "t"
        } else "g"

    r <- new(paste0(kind, shape, "TMatrix"))
    r@Dim <- dims
    if(!missing(dimnames) && !is.null(dimnames))
        r@Dimnames <-
            if(is.character(validDN(dimnames, dims)))
                dimnames
            else fixupDN(dimnames) # needs a valid argument
    if((symmetric || triangular) && all(i >= j))
        r@uplo <- "L" # else "U", the prototype
    if(!m.x) {
        if(is.integer(x))
            x <- as.double(x)
        if((n.x <- length(x)) > 0L && n.x != n.i) {
            if(n.x < n.i) {
                if(n.i %% n.x != 0L)
                    warning(if(m.i) "p[length(p)] " else "length(i) ",
                            "is not an integer multiple of length(x)")
                x <- rep_len(x, n.i) # recycle
            } else if(n.x == 1L)
                x <- x[0L] # tolerate length(i) = 0, length(x) = 1
            else stop("length(x) must not exceed ",
                      if(m.i) "p[length(p)]" else "length(i)")
        }
        if(use.last.ij && n.i == n.j &&
           anyDuplicated.matrix(ij <- cbind(i, j, deparse.level = 0L),
                                fromLast = TRUE)) {
            which.not.dup <- which(!duplicated(ij, fromLast = TRUE))
            i <- i[which.not.dup]
            j <- j[which.not.dup]
            x <- x[which.not.dup]
        }
        r@x <- x
    }
    r@i <- i
    r@j <- j

    if(check)
        validObject(r)
    switch(repr, "C" = .M2C(r), "T" = r, "R" = .M2R(r),
           ## should never happen:
           stop("invalid 'repr'; must be \"C\", \"R\", or \"T\""))
}


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("mean", signature(x = "sparseMatrix"),
          function(x, ...) mean(as(x, "sparseVector"), ...))

### "[<-" : -----------------

## setReplaceMethod("[", .........)
## -> ./Tsparse.R
## &  ./Csparse.R  & ./Rsparse.R {those go via Tsparse}

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

### --- print() and show() methods ---

.formatSparseSimple <- function(m, asLogical=FALSE, digits=NULL,
                                col.names, note.dropping.colnames = TRUE,
                                dn=dimnames(m))
{
    stopifnot(is.logical(asLogical))
    if(asLogical)
        cx <- array("N", dim(m), dimnames=dn)
    else { ## numeric (or --not yet implemented-- complex):
        cx <- apply(m, 2, format, digits=digits)
        if(is.null(dim(cx))) {# e.g. in	1 x 1 case
            dim(cx) <- dim(m)
            dimnames(cx) <- dn
        }
    }
    if (missing(col.names))
        col.names <- {
            if(!is.null(cc <- getOption("sparse.colnames")))
                cc
            else if(is.null(dn[[2]]))
                FALSE
            else { # has column names == dn[[2]]
                ncol(m) < 10
            }
        }
    if(identical(col.names, FALSE))
        cx <- emptyColnames(cx, msg.if.not.empty = note.dropping.colnames)
    else if(is.character(col.names)) {
        stopifnot(length(col.names) == 1)
        cn <- col.names
        switch(substr(cn, 1,3),
               "abb" = {
            iarg <- as.integer(sub("^[^0-9]*", '', cn))
            colnames(cx) <- abbreviate(colnames(cx), minlength = iarg)
        },
        "sub" = {
            iarg <- as.integer(sub("^[^0-9]*", '', cn))
            colnames(cx) <- substr(colnames(cx), 1, iarg)
        },
        stop(gettextf("invalid 'col.names' string: %s", cn), domain=NA))
    }
    ## else: nothing to do for col.names == TRUE
    cx
}## .formatSparseSimple


### NB: Want this to work also for logical or numeric traditional matrix 'x':
formatSparseM <- function(x, zero.print = ".", align = c("fancy", "right"),
                          m = as(x,"matrix"), asLogical=NULL, uniDiag=NULL,
                          digits=NULL, cx, iN0, dn = dimnames(m))
{
    cld <- getClassDef(class(x))
    if(is.null(asLogical)) {
        asLogical <-
            extends1of(cld, c("nsparseMatrix", "indMatrix", "lsparseMatrix")) ||
                                        # simple TRUE/FALSE
            (extends(cld, "matrix") && is.logical(x))
                                        # has NA and (non-)structural FALSE
    }
    if(missing(cx))
        cx <- .formatSparseSimple(m, asLogical=asLogical, digits=digits, dn=dn)
    if(is.null(d <- dim(cx))) {# e.g. in 1 x 1 case
        d <- dim(cx) <- dim(m)
        dimnames(cx) <- dn
    }
    if(missing(iN0))
        iN0 <- 1L + .Call(m_encodeInd, non0ind(x, cld), di = d, FALSE, FALSE)
    ## ne <- length(iN0)
    if(asLogical) {
        cx[m] <- "|"
        if(!extends(cld, "sparseMatrix"))
            x <- as(x,"sparseMatrix")
        if(anyFalse(x@x)) { ## any (x@x == FALSE)
            ## Careful for *non-sorted* Tsparse, e.g. from U-diag
            if(extends(cld, "TsparseMatrix")) {
                ## have no "fast  uniqTsparse():
                x <- as(x, "CsparseMatrix")
                cld <- getClassDef(class(x))
            }
            F. <- is0(x@x)              # the 'FALSE' ones
### FIXME: have  iN0 already above -- *really* need the following ??? --FIXME--
            ij <- non0.i(x, cld, uniqT=FALSE)
            if(extends(cld, "symmetricMatrix")) {
                ## also get "other" triangle
                notdiag <- ij[,1] != ij[,2] # but not the diagonals again
                ij <- rbind(ij, ij[notdiag, 2:1], deparse.level=0)
                F. <-	  c(F., F.[notdiag])
            }
            iN0 <- 1L + .Call(m_encodeInd, ij, di = d, FALSE, FALSE)
            cx[iN0[F.]] <- ":" # non-structural FALSE (or "o", "," , "-" or "f")?
        }
    }
    else if(match.arg(align) == "fancy" && !is.integer(m)) {
        fi <- apply(m, 2, format.info) ## fi[3,] == 0  <==> not expo.

        ## now 'format' the zero.print by padding it with ' ' on the right:
        ## case 1: non-exponent:  fi[2,] + as.logical(fi[2,] > 0)
        ## the column numbers of all 'zero' entries -- (*large*)
        cols <- 1L + (0:(prod(d)-1L))[-iN0] %/% d[1]
        pad <-
            ifelse(fi[3,] == 0,
                   fi[2,] + as.logical(fi[2,] > 0),
                   ## exponential:
                   fi[2,] + fi[3,] + 4)
        ## now be efficient ; sprintf() is relatively slow
        ## and pad is much smaller than 'cols'; instead of "simply"
        ## zero.print <- sprintf("%-*s", pad[cols] + 1, zero.print)
        if(any(doP <- pad > 0)) {       #
            ## only pad those that need padding - *before* expanding
            z.p.pad <- rep.int(zero.print, length(pad))
            z.p.pad[doP] <- sprintf("%-*s", pad[doP] + 1, zero.print)
            zero.print <- z.p.pad[cols]
        }
        else
            zero.print <- rep.int(zero.print, length(cols))
    } ## else "right" : nothing to do
    if(!asLogical && isTRUE(uniDiag)) { ## use "I" in diagonal -- pad correctly
        if(any(diag(x) != 1))
            stop("uniDiag=TRUE, but not all diagonal entries are 1")
        D <- diag(cx) # use
        if(any((ir <- regexpr("1", D)) < 0)) {
            warning("uniDiag=TRUE, not all entries in diagonal coded as 1")
        } else {
            ir <- as.vector(ir)
            nD <- nchar(D, "bytes")
            ## replace "1..." by "I  " (I plus blanks)
            substr(D, ir, nD) <- sprintf("I%*s", nD - ir, "")
            diag(cx) <- D
        }
    }
    cx[-iN0] <- zero.print
    cx
}## formatSparseM()

##' The `format()` method for sparse Matrices; also used inside sparseMatrix print()ing,
##' exported as it might be useful directly.
formatSpMatrix <- function(x, digits = NULL, # getOption("digits"),
                           maxp = 1e9, # ~ 1/2 * .Machine$integer.max, ## getOption("max.print"),
                           cld = getClassDef(class(x)), zero.print = ".",
                           col.names, note.dropping.colnames = TRUE, uniDiag = TRUE,
                           align = c("fancy", "right"))
{
    stopifnot(extends(cld, "sparseMatrix"))
    validObject(x) # have seen seg.faults for invalid objects
    d <- dim(x)
    unitD <- extends(cld, "triangularMatrix") && x@diag == "U"
    ## Will note it is *unit*-diagonal by using "I" instead of "1"
    if(unitD)
        x <- .Call(R_sparse_diag_U2N, x)

    if(maxp < 100) maxp <- 100L # "stop gap"
    if(prod(d) > maxp) { # "Large" => will be "cut"
        ## only coerce to dense that part which won't be cut :
        nr <- maxp %/% d[2]
        m <- as(x[1:max(1, nr), ,drop=FALSE], "matrix")
    } else {
        m <- as(x, "matrix")
    }
    dn <- dimnames(m) ## will be === dimnames(cx)
    binary <- extends(cld,"nsparseMatrix") || extends(cld, "indMatrix") # -> simple T / F
    logi <- binary || extends(cld,"lsparseMatrix") # has NA and (non-)structural FALSE
    cx <- .formatSparseSimple(m, asLogical = logi, digits=digits,
                              col.names=col.names,
                              note.dropping.colnames=note.dropping.colnames, dn=dn)
    if(is.logical(zero.print))
        zero.print <- if(zero.print) "0" else " "
    if(binary) {
        cx[!m] <- zero.print
        cx[m] <- "|"
    } else { # non-binary ==> has 'x' slot
        ## show only "structural" zeros as 'zero.print', not all of them..
        ## -> cannot use 'm' alone
        d <- dim(cx)
        ne <- length(iN0 <- 1L + .Call(m_encodeInd, non0ind(x, cld),
                                       di = d, FALSE, FALSE))
        if(0 < ne && (logi || ne < prod(d))) {
            cx <- formatSparseM(x, zero.print, align, m=m,
                                asLogical = logi, uniDiag = unitD & uniDiag,
                                digits=digits, cx=cx, iN0=iN0, dn=dn)
        } else if (ne == 0)# all zeroes
            cx[] <- zero.print
    }
    cx
}## formatSpMatrix()


## FIXME(?) -- ``merge this'' (at least ``synchronize'') with
## - - -   prMatrix() from ./Auxiliaries.R
## FIXME: prTriang() in ./Auxiliaries.R  should also get  align = "fancy"
##
printSpMatrix <- function(x,
                          digits = NULL, # getOption("digits"),
                          maxp = max(100L, getOption("max.print")),
                          cld = getClassDef(class(x)),
                          zero.print = ".",
                          col.names,
                          note.dropping.colnames = TRUE,
                          uniDiag = TRUE,
                          col.trailer = "",
                          align = c("fancy", "right"))
{
    stopifnot(extends(cld, "sparseMatrix"))
    cx <- formatSpMatrix(x, digits=digits, maxp=maxp, cld=cld,
                         zero.print=zero.print, col.names=col.names,
                         note.dropping.colnames=note.dropping.colnames,
                         uniDiag=uniDiag, align=align)
    if(col.trailer != '')
        cx <- cbind(cx, col.trailer, deparse.level = 0)
    ## right = TRUE : cheap attempt to get better "." alignment
    print(cx, quote = FALSE, right = TRUE, max = maxp)
    invisible(x)
} ## printSpMatrix()

##' The "real" show() / print() method, calling the above printSpMatrix():
printSpMatrix2 <- function(x, digits = NULL, # getOption("digits"),
                           maxp = max(100L, getOption("max.print")), zero.print = ".",
                           col.names, note.dropping.colnames = TRUE, uniDiag = TRUE,
                           suppRows = NULL, suppCols = NULL,
                           col.trailer = if(suppCols) "......" else "",
                           align = c("fancy", "right"),
                           width = getOption("width"), fitWidth = TRUE)
{
    d <- dim(x)
    cl <- class(x)
    cld <- getClassDef(cl)
    xtra <- if(extends(cld, "triangularMatrix") && x@diag == "U")
                " (unitriangular)" else ""
    cat(sprintf('%d x %d sparse Matrix of class "%s"%s\n',
                d[1], d[2], cl, xtra))
    setW <-  !missing(width) && width > getOption("width")
    if(setW) {
        op <- options(width = width) ; on.exit( options(op) ) }
    if((isFALSE(suppRows) && isFALSE(suppCols)) ||
       (!isTRUE(suppRows) && !isTRUE(suppCols) && prod(d) <= maxp))
    { ## "small matrix" and supp* not TRUE : no rows or columns are suppressed
        if(missing(col.trailer) && is.null(suppCols))
            suppCols <- FALSE # for default 'col.trailer'
        printSpMatrix(x, cld=cld,
                      digits=digits,
                      maxp=maxp,
                      zero.print=zero.print,
                      col.names=col.names,
                      note.dropping.colnames=note.dropping.colnames,
                      uniDiag=uniDiag,
                      col.trailer=col.trailer,
                      align=align)
    }
    else { ## d[1] > maxp / d[2] >= nr : -- this needs [,] working:
        validObject(x)
        sTxt <- c(" ", gettext("in show(); maybe adjust 'options(max.print= *, width = *)'"),
                  "\n ..............................\n")
        useW <- width - (format.info(d[1], digits=digits)[1] + 3+1)
        ##  ==  width - space for the largest row label : "[<last>,] "
        ## Suppress rows and/or columns in printing ...
        ## ---------------------------------------- but which exactly depends on format
        ## Determining number of columns - first assuming all zeros : ". . "..: 2 chars/column
        ## i.e., we get the *maximal* numbers of columns to keep, nc :
        if(is.null(suppCols)) # i.e., "it depends" ..
            suppCols <- (d[2] * 2 > useW) # used in 'col.trailer' default
        nCc <- 1 + nchar(col.trailer, "width")
        if(suppCols) {
            nc <- (useW - nCc) %/% 2
            x <- x[ , 1:nc, drop = FALSE]
        } else
            nc <- d[2]
        nr <- maxp %/% nc # if nc becomes smaller,  nr will become larger (!)
        if(is.null(suppRows)) suppRows <- (nr < d[1])
        if(suppRows) {
            n2 <- ceiling(nr / 2)
            nr1 <- min(d[1], max(1L, n2)) #{rows} in 1st part
            nr2 <- max(1L, nr-n2)         #{rows} in 2nd part
            nr <- nr1+nr2 # total #{rows} to be printed
            if(fitWidth) {
                ## one iteration of improving the width, by "fake printing" :
                cM <- formatSpMatrix(x[seq_len(nr1), , drop = FALSE],
                                     digits=digits,
                                     maxp=maxp,
                                     zero.print=zero.print,
                                     col.names=col.names,
                                     align=align,
                                     note.dropping.colnames=note.dropping.colnames,
                                     uniDiag=FALSE)
                ## width needed (without the 'col.trailer's  'nCc'):
                matW <- nchar(capture.output(print(cM, quote=FALSE, right=FALSE))[[1]])
                needW <- matW + (if(suppCols) nCc else 0)
                if(needW > useW) { ## need more width
                    op <- options(width = width+(needW-useW))
                    if(!setW) on.exit( options(op) )
                }
            }
            printSpMatrix(x[seq_len(nr1), , drop=FALSE],
                          digits=digits,
                          maxp=maxp,
                          zero.print=zero.print,
                          col.names=col.names,
                          note.dropping.colnames=note.dropping.colnames,
                          uniDiag=uniDiag,
                          col.trailer = col.trailer, align=align)
            suppTxt <- if(suppCols)
                           gettextf("suppressing %d columns and %d rows", d[2]-nc , d[1]-nr)
                       else gettextf("suppressing %d rows", d[1]-nr)
            cat("\n ..............................",
                "\n ........", suppTxt, sTxt, sep='')
            ## tail() automagically uses "[..,]" rownames:
            printSpMatrix(tail(x, nr2),
                          digits=digits,
                          maxp=maxp,
                          zero.print=zero.print,
                          col.names=col.names,
                          note.dropping.colnames=note.dropping.colnames,
                          uniDiag=FALSE,
                          col.trailer = col.trailer,
                          align=align)
        }
        else if(suppCols) {
            printSpMatrix(x[ , 1:nc , drop = FALSE],
                          digits=digits,
                          maxp=maxp,
                          zero.print=zero.print,
                          col.names=col.names,
                          note.dropping.colnames=note.dropping.colnames,
                          uniDiag=uniDiag,
                          col.trailer = col.trailer,
                          align=align)
            cat("\n .....", gettextf("suppressing %d columns", d[2]-nc), sTxt, sep='')
        }
        else stop("logic programming error in printSpMatrix2(), please report")
        invisible(x)
    }
} ## printSpMatrix2 ()

setMethod("format", signature(x = "sparseMatrix"), formatSpMatrix)

setMethod("print", signature(x = "sparseMatrix"), printSpMatrix2)

setMethod("show", signature(object = "sparseMatrix"),
          function(object) printSpMatrix2(object))



## For very large and very sparse matrices,  the above show()
## is not really helpful;  Use  summary() showing "triplet" as an alternative:

mat2triplet <- function(x, uniqT = FALSE) {
    T <- as(x, "TsparseMatrix")
    if(uniqT && anyDuplicatedT(T)) T <- .uniqTsparse(T)
    if(is(T, "nsparseMatrix"))
         list(i = T@i + 1L, j = T@j + 1L)
    else list(i = T@i + 1L, j = T@j + 1L, x = T@x)
}

setMethod("summary", signature(object = "sparseMatrix"),
          function(object, uniqT = FALSE, ...) {
              d <- dim(object)
              ## return a data frame (int, int,	 {double|logical|...})	:
              r <- as.data.frame(mat2triplet(object, uniqT=uniqT))
              attr(r, "header") <-
                  sprintf('%d x %d sparse Matrix of class "%s", with %d entries',
                          d[1], d[2], class(object), nrow(r))
              ## use ole' S3 technology for such a simple case
              class(r) <- c("sparseSummary", class(r))
              r
          })

print.sparseSummary <- function (x, ...) {
    cat(attr(x, "header"),"\n")
    print.data.frame(x, ...)
    invisible(x)
}

setMethod("dim<-", signature(x = "sparseMatrix"),
          function(x, value) {
              if(!is.numeric(value) || length(value) != 2L)
                  stop("dimensions must be numeric of length 2")
              if(anyNA(value))
                  stop("dimensions cannot contain NA")
              if(any(value < 0))
                  stop("dimensions cannot contain negative values")
              if(!is.integer(value)) {
                  if(any(value > .Machine$integer.max))
                      stop("dimensions cannot exceed 2^31-1")
                  value <- as.integer(value)
              }
              if(all(value == (d <- x@Dim)))
                  return(x)
              if((pv <- prod(value)) != (pd <- prod(d)))
                  stop(gettextf("assigned dimensions [product %.0f] do not match object length [%.0f]",
                                pv, pd, domain = NA))
              r <- spV2M(as(x, "sparseVector"),
                         nrow = value[1L], ncol = value[2L])
              ## 'r' is a TsparseMatrix
              if(extends(cd <- getClassDef(class(x)) , "CsparseMatrix"))
                  as(r, "CsparseMatrix")
              else if(extends(cd, "RsparseMatrix"))
                  as(r, "RsparseMatrix")
              else r
          })

setMethod("rep", "sparseMatrix",
          function(x, ...) rep(as(x, "sparseVector"), ...))

setMethod("cov2cor", signature(V = "sparseMatrix"),
          function(V) {
              ## like stats::cov2cor() but making sure all matrices stay sparse
              p <- (d <- dim(V))[1]
              if (p != d[2])
                  stop("'V' is not a *square* matrix")
              if(!is(V, "dMatrix"))
                  V <- as(V, "dMatrix")# actually "dsparseMatrix"
              Is <- sqrt(1/diag(V))
              if (any(!is.finite(Is))) ## original had 0 or NA
                  warning("diag(.) had 0 or NA entries; non-finite result is doubtful")
              Is <- Diagonal(x = Is)# , names = TRUE
              r <- Is %*% V %*% Is
              r[cbind(1:p,1:p)] <- 1 # exact in diagonal
              as(`dimnames<-`(r, symmDN(dimnames(V))), "symmetricMatrix")
              ## as(r, "symmetricMatrix")
          })

## all.equal(): similar to all.equal_Mat() in ./Matrix.R ;
## -----------	eventually defer to  "sparseVector" methods:
setMethod("all.equal", c(target = "sparseMatrix", current = "sparseMatrix"),
          function(target, current, check.attributes = TRUE, ...) {
              msg <- attr.all_Mat(target, current, check.attributes=check.attributes, ...)
              if(is.list(msg)) msg[[1]]
              else .a.e.comb(msg,
                             all.equal(as(target, "sparseVector"), as(current, "sparseVector"),
                                       check.attributes=check.attributes, ...))
          })
setMethod("all.equal", c(target = "sparseMatrix", current = "ANY"),
          function(target, current, check.attributes = TRUE, ...) {
              msg <- attr.all_Mat(target, current, check.attributes=check.attributes, ...)
              if(is.list(msg)) msg[[1]]
              else .a.e.comb(msg,
                             all.equal(as(target, "sparseVector"), current,
                                       check.attributes=check.attributes, ...))
          })
setMethod("all.equal", c(target = "ANY", current = "sparseMatrix"),
          function(target, current, check.attributes = TRUE, ...) {
              msg <- attr.all_Mat(target, current, check.attributes=check.attributes, ...)
              if(is.list(msg)) msg[[1]]
              else .a.e.comb(msg,
                             all.equal(target, as(current, "sparseVector"),
                                       check.attributes=check.attributes, ...))
          })


setMethod("writeMM", "sparseMatrix",
          function(obj, file, ...)
              writeMM(as(obj, "CsparseMatrix"), as.character(file), ...))

### --- sparse model matrix,  fac2sparse, etc ----> ./spModels.R

###  xtabs(*, sparse = TRUE) ---> part of standard package 'stats' since R 2.10.0

##' @title Random Sparse Matrix
##' @param nrow,
##' @param ncol number of rows and columns, i.e., the matrix dimension
##' @param nnz number of non-zero entries
##' @param rand.x random number generator for 'x' slot
##' @param ... optionally further arguments passed to sparseMatrix()
##' @return a sparseMatrix of dimension (nrow, ncol)
##' @author Martin Maechler
##' @examples M1 <- rsparsematrix(1000, 20, nnz = 200)
##'           summary(M1)
if(FALSE) ## better version below
rsparsematrix <- function(nrow, ncol, nnz,
                          rand.x = function(n) signif(rnorm(nnz), 2),
                          warn.nnz = TRUE, ...) {
    maxi.sample <- 2^31 # maximum n+1 for which sample(n) returns integer
    stopifnot((nnz <- as.integer(nnz)) >= 0,
              nrow >= 0, ncol >= 0, nnz <= nrow * ncol,
              nrow < maxi.sample, ncol < maxi.sample)
    ## to ensure that nnz is strictly followed, must act on duplicated (i,j):
    i <- sample.int(nrow, nnz, replace = TRUE)
    j <- sample.int(ncol, nnz, replace = TRUE)
    dim <- c(nrow, ncol)
    it <- 0
    while((it <- it+1) < 100 &&
          anyDuplicated(n.ij <- encodeInd2(i, j, dim, checkBnds = FALSE))) {
        m <- length(k.dup <- which(duplicated(n.ij)))
        Matrix.msg(sprintf("%3g duplicated (i,j) pairs", m), .M.level = 2)
        if(runif(1) <= 1/2)
            i[k.dup] <- sample.int(nrow, m, replace = TRUE)
        else
            j[k.dup] <- sample.int(ncol, m, replace = TRUE)
    }
    if(warn.nnz && it == 100 &&
       anyDuplicated(encodeInd2(i, j, dim, checkBnds = FALSE)))
        warning("number of non zeros is smaller than 'nnz' because of duplicated (i,j)s")
    sparseMatrix(i = i, j = j, x = rand.x(nnz), dims = dim, ...)
}

## No warn.nnz needed, as we sample the encoded (i,j) with*out* replacement:
rsparsematrix <- function(nrow, ncol, density,
                          nnz = round(density * maxE), symmetric = FALSE,
                          rand.x = function(n) signif(rnorm(n), 2), ...) {
    maxE <- if(symmetric) nrow*(nrow+1)/2 else nrow*ncol
    stopifnot((nnz <- as.integer(nnz)) >= 0,
              nrow >= 0, ncol >= 0, nnz <= maxE)
    ## sampling with*out* replacement (replace=FALSE !):
    ijI <- -1L +
        if(symmetric) sample(indTri(nrow, diag=TRUE), nnz)
        else sample.int(maxE, nnz)
    ## i,j below correspond to  ij <- decodeInd(code, nr) :
    if(is.null(rand.x))
        sparseMatrix(i = ijI  %% nrow,
                     j = ijI %/% nrow,
                     index1 = FALSE,
                     symmetric = symmetric,
                     dims = c(nrow, ncol), ...)
    else
        sparseMatrix(i = ijI  %% nrow,
                     j = ijI %/% nrow,
                     x = rand.x(nnz),
                     index1 = FALSE,
                     symmetric = symmetric,
                     dims = c(nrow, ncol), ...)
}

if(FALSE) ### FIXME: This would *NOT* be needed, if    as.matrix(<sparseMatrix>) was a no-op ;
          ### -----  and then,  base::scale() -> base::scale.default() would work "magically" already..
## scale() is S3 generic in base
scale.sparseMatrix <- function(x, center = FALSE, scale = TRUE) {
    if(center) warning("a sparseMatrix should rarely be centered: will not be sparse anymore")
    ## x <- as.matrix(x)

    ## This rest is *identically*  == base :: scale.default :
    nc <- ncol(x)
    if (is.logical(center)) {
        if (center) {
            center <- colMeans(x, na.rm=TRUE)
            x <- sweep(x, 2L, center, check.margin=FALSE)
        }
    }
    else if (is.numeric(center) && (length(center) == nc))
        x <- sweep(x, 2L, center, check.margin=FALSE)
    else
        stop("length of 'center' must equal the number of columns of 'x'")
    if (is.logical(scale)) {
        if (scale) {
            f <- function(v) {
                v <- v[!is.na(v)]
                sqrt(sum(v^2) / max(1, length(v) - 1L))
            }
            scale <- apply(x, 2L, f)
            x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
        }
    }
    else if (is.numeric(scale) && length(scale) == nc)
        x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
    else
        stop("length of 'scale' must equal the number of columns of 'x'")
    if(is.numeric(center)) attr(x, "scaled:center") <- center
    if(is.numeric(scale)) attr(x, "scaled:scale") <- scale
    x
}

.sparse.diag.get <- function(x, nrow, ncol, names = TRUE)
    .Call(R_sparse_diag_get, x, names)
.sparse.diag.set <- function(x, value)
    .Call(R_sparse_diag_set, x, value)
.sparse.band <- function(x, k1, k2, ...) .Call(R_sparse_band, x, k1, k2)
.sparse.triu <- function(x, k = 0L, ...) .Call(R_sparse_band, x, k, NULL)
.sparse.tril <- function(x, k = 0L, ...) .Call(R_sparse_band, x, NULL, k)
.sparse.t    <- function(x)              .Call(R_sparse_transpose, x, FALSE)
.sparse.fS1  <- function(x, uplo) .Call(R_sparse_force_symmetric, x, NULL)
.sparse.fS2  <- function(x, uplo) .Call(R_sparse_force_symmetric, x, uplo)
.sparse.symmpart <- function(x) .Call(R_sparse_symmpart, x)
.sparse.skewpart <- function(x) .Call(R_sparse_skewpart, x)

.C.is.di <- function(object)
    .Call(Csparse_is_diagonal, object)
.R.is.di <- function(object)
    .Call(Rsparse_is_diagonal, object)
.T.is.di <- function(object)
    .Call(Tsparse_is_diagonal, object)

.C.is.tr <- function(object, upper = NA, ...)
    .Call(Csparse_is_triangular, object, upper)
.R.is.tr <- function(object, upper = NA, ...)
    .Call(Rsparse_is_triangular, object, upper)
.T.is.tr <- function(object, upper = NA, ...)
    .Call(Tsparse_is_triangular, object, upper)

.C.is.sy <- function(object, checkDN = TRUE, ...) {
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    .Call(Csparse_is_symmetric, object, checkDN)
}
.R.is.sy <- function(object, checkDN = TRUE, ...) {
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    .Call(Rsparse_is_symmetric, object, checkDN)
}
.T.is.sy <- function(object, checkDN = TRUE, ...) {
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    .Call(Csparse_is_symmetric, as(object, "CsparseMatrix"), checkDN)
}
.sparse.is.sy.dz <- function(object, tol = 100 * .Machine$double.eps,
                             checkDN = TRUE, ...) {
    ## backwards compatibility: don't check DN if check.attributes=FALSE
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    ## be very fast when requiring exact symmetry
    if(tol <= 0) {
        if(!.hasSlot(object, "p"))
            return(.Call(Csparse_is_symmetric, as(object, "CsparseMatrix"), checkDN))
        else if(.hasSlot(object, "i"))
            return(.Call(Csparse_is_symmetric, object, checkDN))
        else
            return(.Call(Rsparse_is_symmetric, object, checkDN))
    }
    ## pretest: is it square?
    d <- object@Dim
    if((n <- d[1L]) != d[2L])
        return(FALSE)
    ## pretest: are DN symmetric in the sense of validObject(<symmetricMatrix>)?
    if(checkDN && !isSymmetricDN(object@Dimnames))
        return(FALSE)
    if(n <= 1L)
        return(TRUE)
    ## now handling an n-by-n [CRT]sparseMatrix, n >= 2:
    x  <- as(  object,  "sparseVector")
    tx <- as(t(object), "sparseVector")
    if(is(tx, "zsparseVector"))
        tx@x <- Conj(tx@x)
    ae <- function(check.attributes, ...) {
        ## discarding possible user-supplied check.attributes:
        all.equal(..., check.attributes = FALSE)
    }
    isTRUE(ae(target = x, current = tx, tolerance = tol, ...))
}

.sparse.subclasses <- names(getClassDef("sparseMatrix")@subclasses)

for (.cl in grep("^[CRT]sparseMatrix$", .sparse.subclasses, value = TRUE)) {
    setMethod("diag",   signature(x = .cl), .sparse.diag.get)
    setMethod("diag<-", signature(x = .cl), .sparse.diag.set)
    setMethod("band", signature(x = .cl), .sparse.band)
    setMethod("triu", signature(x = .cl), .sparse.triu)
    setMethod("tril", signature(x = .cl), .sparse.tril)
    setMethod("t",    signature(x = .cl), .sparse.t)
    setMethod("forceSymmetric", signature(x = .cl, uplo = "missing"),
              .sparse.fS1)
    setMethod("forceSymmetric", signature(x = .cl, uplo = "character"),
              .sparse.fS2)
    setMethod("symmpart", signature(x = .cl), .sparse.symmpart)
    setMethod("skewpart", signature(x = .cl), .sparse.skewpart)
    setMethod("isDiagonal", signature(object = .cl),
              get(paste0(".", substr(.cl, 1L, 1L), ".is.di"),
                  mode = "function", inherits = FALSE))
}

for (.cl in grep("^.g[CRT]Matrix$", .sparse.subclasses, value = TRUE))
    setMethod("isTriangular", signature(object = .cl),
              get(paste0(".", substr(.cl, 3L, 3L), ".is.tr"),
                  mode = "function", inherits = FALSE))
for (.cl in grep("^[lni]g[CRT]Matrix$", .sparse.subclasses, value = TRUE))
    setMethod("isSymmetric", signature(object = .cl),
              get(paste0(".", substr(.cl, 3L, 3L), ".is.sy"),
                  mode = "function", inherits = FALSE))
for (.cl in grep("^[dz][gt][CRT]Matrix$", .sparse.subclasses, value = TRUE))
    setMethod("isSymmetric", signature(object = .cl), .sparse.is.sy.dz)

rm(.cl, .sparse.subclasses, .sparse.is.sy.dz,
   list = c(grep("^[.]sparse[.](band|tri[ul]|t|fS[21]|symmpart|skewpart)$",
                 ls(all.names = TRUE), value = TRUE),
            grep("^[.][CRT][.]is[.](di|tr|sy)$",
                 ls(all.names = TRUE), value = TRUE)))
