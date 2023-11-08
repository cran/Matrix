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
                           .diag2sparse(data, ".", "s", "C", "U")
                       else .diag2dense(data, ".", "s", FALSE, "U"))
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

sparseMatrix <- function(i, j, p, x, dims, dimnames,
                         symmetric = FALSE,
                         triangular = FALSE,
                         index1 = TRUE,
                         repr = c("C", "R", "T"),
                         giveCsparse,
                         check = TRUE,
                         use.last.ij = FALSE)
{
    if((m.i <- missing(i)) + (m.j <- missing(j)) + (m.p <- missing(p)) != 1L)
        stop("exactly one of 'i', 'j', and 'p' must be missing from call")
    if(symmetric && triangular)
        stop("use Diagonal() to construct diagonal (symmetric && triangular) sparse matrices")
    index1 <- as.logical(index1) # allowing {0,1}

    repr <- # keep in sync with toeplitz(<sparseVector>)
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

    kind <- if(m.x <- missing(x)) "n" else if(is.integer(x)) "d" else .M.kind(x)
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

spMatrix <- function(nrow, ncol,
                     i = integer(0L), j = integer(0L), x = double(0L))
    new(paste0(if(is.integer(x)) "d" else .M.kind(x), "gTMatrix"),
        Dim = c(as.integer(nrow), as.integer(ncol)),
        i = as.integer(i) - 1L,
        j = as.integer(j) - 1L,
        x = if(is.integer(x)) as.double(x) else x)

Diagonal <- function(n, x = NULL, names = FALSE)
{
    nx <- length(x)
    if(missing(n))
        n <- nx
    else if(!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0L)
        stop("'n' must be a non-negative integer")
    if(is.double(n) && n >= .Machine$integer.max + 1)
        stop("dimensions cannot exceed 2^31-1")
    n <- as.integer(n) # discarding attributes
    if(is.null(x)) {
        r <- new("ddiMatrix")
        r@diag <- "U"
        if(n > 0L) {
            r@Dim <- c(n, n)
            if(is.character(names) && length(names) == n)
                r@Dimnames <- list(names, names)
        }
        return(r)
    }
    if(is.object(x))
        stop(gettextf("'x' has unsupported class \"%s\"", class(x)[1L]),
             domain = NA)
    names.x <- names(x) # keeping for later
    r <- new(switch(typeof(x),
                    ## discarding attributes, incl. 'dim' and 'names'
                    logical = { x <- as.logical(x); "ldiMatrix" },
                    integer =,
                    double = { x <- as.double(x); "ddiMatrix" },
                    stop(gettextf("'x' has unsupported type \"%s\"", typeof(x)),
                         domain = NA)))
    if(n == 0L)
        return(r)
    if(nx != 1L)
        r@x <-
            if(nx == n)
                x
            else if(nx > 0L)
                rep_len(x, n)
            else stop("attempt to recycle 'x' of length 0 to length 'n' (n > 0)")
    else if(is.na(x) || x != 1)
        r@x <- rep.int(x, n)
    else r@diag <- "U"
    r@Dim <- c(n, n)
    if(is.character(names)) {
        if(length(names) == n)
            r@Dimnames <- list(names, names)
    } else if(isTRUE(names) && !is.null(names.x)) {
        names.x <- rep_len(names.x, n) # we know length(names.x) > 0L
        r@Dimnames <- list(names.x, names.x)
    }
    r
}

.sparseDiagonal <- function(n, x = NULL, uplo = "U", shape = "t",
                            unitri = TRUE, kind, cols)
{
    if(missing(n))
        n <- length(x)
    else if(!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0L)
        stop("'n' must be a non-negative integer")
    if(is.double(n) && n >= .Machine$integer.max + 1)
        stop("dimensions cannot exceed 2^31-1")
    n <- nj <- as.integer(n) # stripping attributes

    if(!(missing(shape) ||
         (is.character(shape) && length(shape) == 1L && !is.na(shape) &&
          any(shape == c("g", "t", "s")))))
        stop("'shape' must be one of \"g\", \"t\", \"s\"")

    if(!((m.kind <- missing(kind)) ||
         (is.character(kind) && length(kind) == 1L && !is.na(kind) &&
          any(kind == c("d", "l", "n")))))
        stop("'kind' must be one of \"d\", \"l\", \"n\"")

    if(m.kind || kind != "n") {
        if(is.null(x))
           x <- if(m.kind) { kind <- "d"; 1 } else switch(kind, d = 1, l = TRUE)
        else if(is.object(x))
            stop(gettextf("'x' has unsupported class \"%s\"",
                          class(x)[1L]),
                 domain = NA)
        else {
            kind. <- switch(typeof(x),
                            ## discarding attributes, incl. 'dim' in array case
                            logical = { x <- as.logical(x); "l" },
                            integer =,
                            double = { x <- as.double(x); "d" },
                            stop(gettextf("'x' has unsupported type \"%s\"",
                                          typeof(x)),
                                 domain = NA))
            if(m.kind)
                kind <- kind.
            else if(kind != kind.) {
                warning(gettextf("mismatch between typeof(x)=\"%s\" and kind=\"%s\"; using kind=\"%s\"",
                                 typeof(x), kind, kind.),
                        domain = NA)
                kind <- kind.
            }
        }
    }

    if(!(m.cols <- missing(cols))) {
        if(!is.numeric(cols))
            stop("'cols' must be numeric")
        else if((nj <- length(cols)) > 0L &&
                (n == 0L || anyNA(rj <- range(cols)) ||
                 rj[1L] < 0L || rj[2L] >= n))
            stop("'cols' has elements not in seq(0, length.out = n)")
        else {
            cols <- as.integer(cols)
            shape <- "g"
        }
    }

    r <- new(paste0(kind, shape, "CMatrix"))
    r@Dim <- c(n, nj)
    if(shape != "g") {
        if(!missing(uplo)) {
            if(is.character(uplo) && length(uplo) == 1L && !is.na(uplo) &&
               any(uplo == c("U", "L")))
                r@uplo <- uplo
            else stop("'uplo' must be \"U\" or \"L\"")
        }
        if(shape == "t" && unitri &&
           (kind == "n" || (!anyNA(x) && all(if(kind == "l") x else x == 1)))) {
            r@diag <- "U"
            r@p <- integer(nj + 1)
            return(r)
        }
    }
    if(nj > 0L) {
        r@p <- 0:nj
        r@i <- if(m.cols) 0:(nj - 1L) else cols
        if(kind != "n") {
            x <-
                if((nx <- length(x)) == n)
                    x
                else if(nx > 0L)
                    rep_len(x, n)
                else stop("attempt to recycle 'x' of length 0 to length 'n' (n > 0)")
            r@x <- if(m.cols) x else x[1L + cols]
        }
    }
    r
}

.trDiagonal <- function(n, x = NULL, uplo = "U", unitri = TRUE, kind)
    .sparseDiagonal(n, x, uplo, shape = "t", unitri = unitri, kind = kind)

.symDiagonal <- function(n, x = NULL, uplo = "U", kind)
    .sparseDiagonal(n, x, uplo, shape = "s", kind = kind)

.bdiag <- function(lst)
{
    if(!is.list(lst))
        stop("'lst' must be a list")
    if((n <- length(lst)) == 0L)
        return(new("dgTMatrix"))
    if(n == 1L)
        return(.M2T(asCspN(lst[[1L]])))

### FIXME? this is _slow_ when 'lst' is list of 75000 3-by-3 dense matrices
    lst <- unname(lapply(lst, function(x) .M2T(asCspN(x))))

    cl <- vapply(lst, class, "")
    kind  <- substr(cl, 1L, 1L) # "n", "l", or "d"
    shape <- substr(cl, 2L, 2L) # "g", "s", or "t"

    if(!(any(kind == (kind. <- "d")) || any(kind == (kind. <- "l"))))
        kind. <- "n"
    else if(any(z <- kind == "n"))
        lst[z] <- lapply(lst[z], .sparse2kind, kind.)

    shape. <-
        if(all(symmetric <- shape == "s"))
            "s"
        else if(all(shape == "t"))
            "t"
        else "g"

    if(shape. != "g") {
        uplo <- vapply(lst, slot, "", "uplo") # "U" or "L"
        if(shape. == "s")
            uplo. <-
                if(all(z <- uplo == "U"))
                    "U"
                else if(!any(z))
                    "L"
                else {
                    uplo.. <- if(2 * sum(z) >= n) { z <- !z; "U" } else "L"
                    lst[z] <- lapply(lst[z], .tCRT)
                    uplo..
                }
        else if(any(uplo != (uplo. <- uplo[1L])))
            shape. <- "g"
    }

    i_off <- c(0L, cumsum(vapply(lst, function(x) x@Dim[1L], 0L)))
    j_off <- c(0L, cumsum(vapply(lst, function(x) x@Dim[2L], 0L)))

    r <- new(paste0(kind., shape., "TMatrix"))
    r@Dim <- r@Dim <- c(i_off[n + 1L], j_off[n + 1L])
    if(shape. == "g")
        lst[symmetric] <- lapply(lst[symmetric], .sparse2g)
    else r@uplo <- uplo.
    r@i <- unlist(lapply(seq_len(n), function(k) i_off[k] + lst[[k]]@i),
                  FALSE, FALSE)
    r@j <- unlist(lapply(seq_len(n), function(k) j_off[k] + lst[[k]]@j),
                  FALSE, FALSE)
    if(kind. != "n")
        r@x <- unlist(lapply(lst, slot, "x"), FALSE, FALSE)
    r
}

bdiag <- function(...)
{
    if((n <- ...length()) == 0L)
        new("dgCMatrix")
    else if(n > 1L)
        .M2C(.bdiag(list(...)))
    else if(!is.list(x <- ..1))
        as(x, "CsparseMatrix")
    else if(length(x) == 1L)
        as(x[[1L]], "CsparseMatrix")
    else .M2C(.bdiag(x))
}

bandSparse <- function(n, m = n, k, diagonals,
                       symmetric = FALSE,
                       repr = "C", giveCsparse = (repr == "C"))
{
    ## Purpose: Compute a band-matrix by speciyfying its (sub-)diagonal(s)
    ## ----------------------------------------------------------------------
    ## Arguments: (n,m) : Matrix dimension
    ##                k : integer vector of "diagonal numbers",  with identical
    ##                    meaning as in  band(*, k)
    ##         diagonals: (optional!) list of (sub/super)diagonals
    ##         symmetric: if TRUE, specify only upper or lower triangle;
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 20 Feb 2009, 22:42

    if(use.x <- !missing(diagonals)) # when specified, must be matrix or list
        diag.isMat <- is.matrix(diagonals)
    len.k <- length(k)
    stopifnot(!use.x || is.list(diagonals) || diag.isMat,
              k == as.integer(k), n == as.integer(n), m == as.integer(m))
    k <- as.integer(k)
    n <- as.integer(n)
    m <- as.integer(m)
    stopifnot(n >= 0, m >= 0, -n+1 <= (mik <- min(k)), (mak <- max(k)) <= m - 1)
    if(missing(repr) && !giveCsparse) {
        warning("'giveCsparse' has been deprecated; setting 'repr = \"T\"' for you")
        repr <- "T"
    } else if(!missing(repr) && !missing(giveCsparse))
        warning("'giveCsparse' has been deprecated; will use 'repr' instead")
    if(use.x) {
        if(diag.isMat) {
            if(ncol(diagonals) != len.k)
                stop(gettextf("'diagonals' matrix must have %d columns (= length(k) )",
                              len.k), domain=NA)
            getD <- function(j) diagonals[,j]

        } else { ## is.list(diagonals):
            if(length(diagonals) != len.k)
                stop(gettextf("'diagonals' must have the same length (%d) as 'k'",
                              len.k), domain=NA)
            getD <- function(j) diagonals[[j]]
        }
    }
    sqr <- n == m
    if(symmetric) {
        if(!sqr) stop("matrix can only be symmetric if square, but n != m")
        if(mik < 0 && mak > 0)
            stop("for symmetric band matrix, only specify upper or lower triangle\n hence, all k must have the same sign")
    } else
        tri <- sqr && sign(mik)*sign(mak) >= 0 # triangular result
    dims <- c(n,m)
    k.lengths <- ## This is a bit "ugly"; I got the cases "by inspection"
        if(n >= m) {
            ifelse(k >= m-n,  m - pmax(0,k), n+k)
        } else { ## n < m (?? k >= -n+1 always !!)
            ifelse(k >= -n+1, n + pmin(0,k), m-k)
        }
    i <- j <- integer(sum(k.lengths))
    if(use.x)
        x <- if(len.k > 0) # carefully getting correct type/mode
                 rep.int(getD(1)[1], length(i))
    off.i <- 0L
    for(s in seq_len(len.k)) {
        kk <- k[s] ## *is* integer
        l.kk <- k.lengths[s] ## == length of (sub-)diagonal kk
        ii1 <- seq_len(l.kk)
        ind <- ii1 + off.i
        if(kk >= 0) {
            i[ind] <- ii1
            j[ind] <- ii1 + kk
        } else { ## k < 0
            i[ind] <- ii1 - kk
            j[ind] <- ii1
        }
        if(use.x) {
            xx <- getD(s)
            if(length(xx) < l.kk)
                warning(gettextf("the %d-th (sub)-diagonal (k = %d) is too short; filling with NA's",
                                 s, kk), domain=NA)
            x[ind] <- xx[ii1]
        }
        off.i <- off.i + l.kk
    }
    if(symmetric) { ## we should have smarter sparseMatrix()
        UpLo <- if(min(k) >= 0) "U" else "L"
        T <-
            if(use.x) {
                if(is.integer(x))
                    x <- as.double(x)
                cc <- paste0(.M.kind(x), "sTMatrix")
                new(cc, i= i-1L, j= j-1L, x = x, Dim= dims, uplo=UpLo)
            } else new("nsTMatrix", i= i-1L, j= j-1L, Dim= dims, uplo=UpLo)
        switch(repr, "C" = .M2C(T), "T" = T, "R" = .M2R(T),
               stop("invalid 'repr'; must be \"C\", \"T\", or \"R\""))
    }
    else { ## not symmetric, possibly triangular
        if(use.x)
            sparseMatrix(i=i, j=j, x=x, dims=dims, triangular=tri, repr=repr)
        else
            sparseMatrix(i=i, j=j,	dims=dims, triangular=tri, repr=repr)
    }
}

rsparsematrix <- function(nrow, ncol, density,
                          nnz = round(density * maxE),
                          symmetric = FALSE,
                          rand.x = function(n) signif(rnorm(n), 2L), ...)
{
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

Hilbert <- function(n)
{
    n <- as.integer(n)
    i <- seq_len(n)
    new("dpoMatrix", Dim = c(n, n), x = c(1/outer(i - 1L, i, `+`)))
}

spV2M <- function(x, nrow, ncol, byrow = FALSE,
                  check = TRUE, symmetric = FALSE)
{
    if(check && !is(x, "sparseVector"))
	stop("'x' must inherit from \"sparseVector\"")
    if(!missing(ncol)) { ncol <- as.integer(ncol)
			 if(ncol < 0) stop("'ncol' must be >= 0") }
    if(!missing(nrow)) { nrow <- as.integer(nrow)
			 if(nrow < 0) stop("'nrow' must be >= 0") }
    n <- length(x)
    if(symmetric) {
	if(missing(nrow)) stop("Must specify 'nrow' when 'symmetric' is true")
	if(!missing(ncol) && nrow != ncol)
	    stop("'nrow' and 'ncol' must be the same when 'symmetric' is true")
	## otherwise  ncol will not used at all when (symmetric)
	if(check && as.double(nrow)^2 != n)
	    stop("'x' must have length nrow^2 when 'symmetric' is true")
	## x <- x[indTri(nrow, upper=TRUE, diag=TRUE)]
    } else if(missing(nrow)) {
	nrow <- as.integer(
	    if(missing(ncol)) { ## both missing: --> (n x 1)
		ncol <- 1L
		n
	    } else {
		if(n %% ncol != 0) warning("'ncol' is not a factor of length(x)")
		as.integer(ceiling(n / ncol))
	    })
    } else if(missing(ncol)) {
        ncol <- if(symmetric) nrow else {
            if(n %% nrow != 0) warning("'nrow' is not a factor of length(x)")
            as.integer(ceiling(n / nrow)) }
    } else {                          ## both nrow and ncol specified
        n.n <- as.double(ncol) * nrow # no integer overflow
        if(n.n <  n) stop("nrow * ncol < length(x)", domain = NA)
        if(n.n != n) warning("nrow * ncol != length(x)", domain = NA)
    }
    ## now nrow * ncol >= n  (or 'symmetric')
    ##	   ~~~~~~~~~~~~~~~~
    kind <- .M.kind(x) # "d", "n", "l", "i", "z", ...
    has.x <- kind != "n"
    clStem <- if(symmetric) "sTMatrix" else "gTMatrix"
    ## "careful_new()" :
    cNam <- paste0(kind, clStem)
    chngCl <- is.null(newCl <- getClassDef(cNam))
    if(chngCl) { ## e.g. "igTMatrix" is not yet implemented
	if(kind == "z")
	    stop(gettextf("Class %s is not yet implemented", dQuote(cNam)),
		 domain = NA)
	## coerce to "double":
	newCl <- getClassDef(paste0("d", clStem))
    }
    r <- new(newCl, Dim = c(nrow, ncol))
    ## now "compute"  the (i,j,x) slots given x@(i,x)
    i0 <- x@i - 1L
    if(byrow) { ## need as.integer(.) since <sparseVector> @ i can be double
	j <- as.integer(i0 %% ncol)
	i <- as.integer(i0 %/% ncol)
    } else { ## default{byrow = FALSE}
	i <- as.integer(i0 %% nrow)
	j <- as.integer(i0 %/% nrow)
    }
    if(has.x)
	x <- if(chngCl) as.numeric(x@x) else x@x
    if(symmetric) {  ## using  uplo = "U"
	i0 <- i <= j ## i.e., indTri(nrow, upper=TRUE, diag=TRUE)
	i <- i[i0]
	j <- j[i0]
	if(has.x) x <- x[i0]
    }
    r@j <- j
    r@i <- i
    if(has.x) r@x <- x
    r
}

.sparseV2Mat <- function(from)
    spV2M(from, nrow = from@length, ncol = 1L, check = FALSE)

sp2vec <- function(x, mode = .type.kind[.M.kind(x)])
{
    ## sparseVector  ->  vector
    has.x <- .hasSlot(x, "x")## has "x" slot
    m.any <- (mode == "any")
    if(m.any)
	mode <- if(has.x) mode(x@x) else "logical"
    else if(has.x) # is.<mode>() is much faster than inherits() | is():
        xxOk <- switch(mode,
		       "double" = is.double(x@x),
		       "logical" = is.logical(x@x),
		       "integer" = is.integer(x@x),
		       "complex" = is.complex(x@x),
		       ## otherwise (does not happen with default 'mode'):
		       inherits(x@x, mode))
    r <- vector(mode, x@length)
    r[x@i] <-
	if(has.x) {
	    if(m.any || xxOk) x@x else as(x@x, mode)
	} else TRUE
    r
}

newSpV <- function(class, x, i, length, drop0 = TRUE, checkSort = TRUE)
{
    if(has.x <- !missing(x)) {
	if(length(x) == 1 && (li <- length(i)) != 1) ## recycle x :
	    x <- rep.int(x, li)
	if(drop0 && isTRUE(any(x0 <- x == 0))) {
	    keep <- is.na(x) | !x0
	    x <- x[keep]
	    i <- i[keep]
	}
    }
    if(checkSort && is.unsorted(i)) {
	ii <- sort.list(i)
	if(has.x) x <- x[ii]
	i <- i[ii]
    }
    if(has.x)
	new(class, x = x, i = i, length = length)
    else
	new(class,        i = i, length = length)
}

newSpVec <- function(class, x, prev)
    newSpV(class = class, x = x, i = prev@i, length = prev@length)

sparseVector <- function(x, i, length)
    newSpV(class = paste0(if(missing(x)) "n" else .M.kind(x), "sparseVector"),
           x = x, i = i, length = length)
