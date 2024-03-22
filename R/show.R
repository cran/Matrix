prMatrix <- function(x,
                     digits = getOption("digits"),
                     maxp = getOption("max.print"),
                     ...) {
    d <- dim(x)
    cl <- class(x) ## cld <- getClassDef(cl)
    tri <- extends(cl, "triangularMatrix")
    xtra <- if(tri && x@diag == "U") " (unitriangular)" else ""
    cat(sprintf('%d x %d Matrix of class "%s"%s\n',
                d[1], d[2], cl, xtra))
    if(prod(d) <= maxp) {
        if(tri)
            prTriang(x, digits = digits, maxp = maxp)
        else
            print(as(x, "matrix"), digits = digits, max = maxp)
    }
    else { ## d[1] > maxp / d[2] >= nr :
        m <- as(x, "matrix")
        nr <- maxp %/% d[2]
        n2 <- ceiling(nr / 2)
        print(head(m, max(1, n2)))
        cat("\n ..........\n\n")
        print(tail(m, max(1, nr - n2)))
        cat("\n ..........\n\n")

    }
    ## DEBUG: cat("str(.):\n") ; str(x)
    invisible(x)# as print() S3 methods do
}

prTriang <- function(x,
                     digits = getOption("digits"),
                     maxp = getOption("max.print"),
                     justify = "none",
                     right = TRUE,
                     ...) {
    ## modeled along stats:::print.dist
    upper <- x@uplo == "U"
    m <- as(x, "matrix")
    cf <- format(m, digits = digits, justify = justify)
    cf[if(upper) row(cf) > col(cf)
       else row(cf) < col(cf)] <- "."
    print(cf, quote = FALSE, right = right, max = maxp)
    invisible(x)
}

prDiag <- function(x,
                   digits = getOption("digits"),
                   justify = "none",
                   right = TRUE,
                   ...) {
    cf <- array(".", dim = x@Dim, dimnames = x@Dimnames)
    cf[row(cf) == col(cf)] <-
        vapply(diag(x), format, "", digits = digits, justify = justify)
    print(cf, quote = FALSE, right = right)
    invisible(x)
}

emptyColnames <- function(x, msg.if.not.empty = FALSE) {
    ## Useful for compact printing of (parts) of sparse matrices
    ## possibly	 dimnames(x) "==" NULL :
    if((nd <- length(d <- dim(x))) < 2L)
        return(x)
    nc <- d[2L]
    if(is.null(dn <- dimnames(x)))
        dn <- vector("list", nd)
    else if(msg.if.not.empty &&
            is.character(cn <- dn[[2L]]) &&
            any(nzchar(cn)))
        message(gettextf("  [[ suppressing %d column name%s %s ... ]]",
                         nc,
                         if(nc == 1L) "" else "s",
                         paste0(sQuote(if(nc <= 3L) cn else cn[1:3]),
                                collapse = ", ")),
                domain = NA)
    dn[[2L]] <- character(nc)
    dimnames(x) <- dn
    x
}

.formatSparseSimple <- function(m,
                                asLogical = FALSE,
                                digits = NULL,
                                col.names,
                                note.dropping.colnames = TRUE,
                                dn = dimnames(m)) {
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
}

## NB: Want this to work also for logical or numeric traditional matrix 'x':
formatSparseM <- function(x,
                          zero.print = ".",
                          align = c("fancy", "right"),
                          m = as(x, "matrix"),
                          asLogical = NULL,
                          uniDiag = NULL,
                          digits = NULL,
                          cx,
                          iN0,
                          dn = dimnames(m)) {
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
}

## The `format()` method for sparse Matrices;
## also used inside sparseMatrix print()ing,
## exported as it might be useful directly.
formatSpMatrix <- function(x,
                           digits = NULL,
                           maxp = 1e+09, # ~ 0.5 * .Machine$integer.max
                           cld = getClassDef(class(x)), zero.print = ".",
                           col.names,
                           note.dropping.colnames = TRUE,
                           uniDiag = TRUE,
                           align = c("fancy", "right"),
                           ...) {
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
}


## FIXME(?) -- ``merge this'' (at least ``synchronize'') with
## - - -   prMatrix() from ./Auxiliaries.R
## FIXME: prTriang() in ./Auxiliaries.R  should also get  align = "fancy"
printSpMatrix <- function(x,
                          digits = NULL,
                          maxp = max(100L, getOption("max.print")),
                          cld = getClassDef(class(x)),
                          zero.print = ".",
                          col.names,
                          note.dropping.colnames = TRUE,
                          uniDiag = TRUE,
                          col.trailer = "",
                          align = c("fancy", "right"),
                          ...) {
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
}

##' The "real" show() / print() method, calling the above printSpMatrix():
printSpMatrix2 <- function(x,
                           digits = NULL,
                           maxp = max(100L, getOption("max.print")),
                           zero.print = ".",
                           col.names,
                           note.dropping.colnames = TRUE,
                           uniDiag = TRUE,
                           suppRows = NULL,
                           suppCols = NULL,
                           col.trailer = if(suppCols) "......" else "",
                           align = c("fancy", "right"),
                           width = getOption("width"),
                           fitWidth = TRUE,
                           ...) {
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
        sTxt <- c(" ", gettext("in show(); maybe adjust options(max.print=, width=)"),
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
}

prSpVector <- function(x,
                       digits = getOption("digits"),
                       maxp = getOption("max.print"),
                       zero.print = ".",
                       ...)
{
    cld <- getClassDef(class(x))
    stopifnot(extends(cld, "sparseVector"), maxp >= 1)
    if(is.logical(zero.print))
        zero.print <- if(zero.print) "0" else " "
    ## kind <- .M.kindC(cld)
    ## has.x <- kind != "n"
    n <- x@length
    if(n > 0) {
        if(n > maxp) {
            ## n > maxp =: nn : will cut length of what we'll display :
            x <- head(x, maxp)
            n <- maxp
        }
        xi <- x@i
        is.n <- extends(cld, "nsparseVector")
        logi <- is.n || extends(cld, "lsparseVector")
        cx <- if(logi) rep.int("N", n) else character(n)
        cx[if(length(xi)) -xi else TRUE] <- zero.print
        cx[xi] <-
            if(is.n)
                "|"
            else if(logi)
                c(":", "|")[x@x + 1L]
            else
                ## numeric (or --not yet-- complex): 'has.x' in any cases
                format(x@x, digits = digits)
        ## right = TRUE : cheap attempt to get better "." alignment
        print(cx, quote = FALSE, right = TRUE, max = maxp)
    }
    invisible(x) # TODO? in case of n > maxp, "should" return original x
}


## METHODS FOR GENERIC: show
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("show", c(object = "denseMatrix"),
          function(object) prMatrix(object))

setMethod("show", c(object = "sparseMatrix"),
          function(object) printSpMatrix2(object))

setMethod("show", c(object = "diagonalMatrix"),
          function(object) {
              d <- dim(object)
              cl <- class(object)
              cat(sprintf('%d x %d diagonal matrix of class "%s"',
                          d[1L], d[2L], cl))
              if(d[1L] < 50) {
                  cat("\n")
                  prDiag(object)
              } else {
                  cat(", with diagonal entries\n")
                  show(diag(object))
                  invisible(object)
              }
          })

setMethod("show", c(object = "MatrixFactorization"),
          function(object) {
              cat("matrix factorization of ")
              str(object)
          })

setMethod("show", c(object = "CholeskyFactorization"),
          function(object) {
              cat("Cholesky factorization of ")
              str(object)
          })

setMethod("show", c(object = "BunchKaufmanFactorization"),
          function(object) {
              cat("Bunch-Kaufman factorization of ")
              str(object)
          })

setMethod("show", c(object = "SchurFactorization"),
          function(object) {
              cat("Schur factorization of ")
              str(object)
          })

setMethod("show", c(object = "LU"),
          function(object) {
              cat("LU factorization of ")
              str(object)
          })

setMethod("show", c(object = "QR"),
          function(object) {
              cat("QR factorization of ")
              str(object)
          })

setMethod("show", c(object = "sparseVector"),
          function(object) {
              n <- object@length
              cl <- class(object)
              cat(sprintf("sparse vector (nnz/length = %d/%.0f) of class \"%s\"\n",
                          length(object@i), as.double(n), cl))
              maxp <- max(1, getOption("max.print"))
              if(n <= maxp)
                  prSpVector(object, maxp = maxp)
              else {
                  ## n > maxp : will cut length of what we'll display :
                  ## cannot easily show head(.) & tail(.) because of
                  ## "[1] .." printing of tail
                  prSpVector(head(object, maxp), maxp = maxp)
                  cat(" ............................\n",
                      " ........suppressing ", n - maxp,
                      " entries in show(); maybe adjust options(max.print=)\n",
                      " ............................\n\n",
                      sep = "")
              }
              invisible(object)
          })


## METHODS FOR GENERIC: print
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("print", c(x = "sparseMatrix"),
          printSpMatrix2)

setMethod("print", c(x = "diagonalMatrix"),
          prDiag)


## METHODS FOR GENERIC: format
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("format", c(x = "sparseMatrix"),
          formatSpMatrix)


## METHODS FOR GENERIC: summary
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("summary", c(object = "sparseMatrix"),
          function(object, uniqT = FALSE, ...) {
              d <- object@Dim
              ## return a data frame (int, int, {double|logical}) :
              r <- as.data.frame(mat2triplet(object, uniqT = uniqT))
              attr(r, "header") <-
                  sprintf("%d x %d sparse Matrix of class \"%s\", with %d entries",
                          d[1L], d[2L], class(object), nrow(r))
              class(r) <- c("sparseSummary", oldClass(r))
              r
          })

setMethod("summary", c(object = "diagonalMatrix"),
          function(object, ...) {
              d <- object@Dim
              r <- summary(object@x, ...)
              attr(r, "header") <-
                  sprintf("%d x %d diagonal Matrix of class \"%s\"",
                          d[1L], d[2L], class(object))
              class(r) <- c("diagSummary", class(r))
              r
          })

print.sparseSummary <- print.diagSummary <-
function (x, ...) {
    cat(attr(x, "header"), "\n", sep = "")
    NextMethod()
    invisible(x)
}
