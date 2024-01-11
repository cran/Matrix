#### "Namespace private" Auxiliaries  such as method functions
#### (called from more than one place --> need to be defined early)

is0  <- function(x) !(is.na(x) | x)
isN0 <- function(x)   is.na(x) | x
is1  <- function(x)  !is.na(x) & x == 1L
isN1 <- function(x)   is.na(x) | x != 1L

allTrue  <- function(x)
                                           !is.na(a <- all( x)) &&  a
allFalse <- function(x)
    if(is.atomic(x)) .Call(R_all0, x) else !is.na(a <- any( x)) && !a
all0     <- function(x)
    if(is.atomic(x)) .Call(R_all0, x) else !is.na(a <- all(!x)) &&  a
anyFalse <- function(x)
    if(is.atomic(x)) .Call(R_any0, x) else !is.na(a <- all( x)) && !a
any0     <- function(x)
    if(is.atomic(x)) .Call(R_any0, x) else !is.na(a <- any(!x)) &&  a

## NB: change to using 'typeof' when we define iMatrix
as1 <- function(x, mode. = mode(x))
    switch(mode.,
           "logical" = TRUE,
           "integer" =   1L,
           "double"  =     ,
           "numeric" =    1,
           "complex" = 1+0i,
           stop(gettextf("invalid mode \"%s\"", mode.), domain = NA))

## NB: change to using 'typeof' when we define iMatrix
as0 <- function(x, mode. = mode(x))
    switch(mode.,
           "logical" = FALSE,
           "integer" =    0L,
           "double"  =      ,
           "numeric" =     0,
           "complex" =  0+0i,
           stop(gettextf("invalid mode \"%s\"", mode.), domain = NA))

.bail.out.2 <- function(name, cl1, cl2)
    stop(gettextf("%s(<%s>, <%s>) is not yet implemented; ask maintainer(\"%s\") to implement the missing method",
                  name, cl1[1L], cl2[1L], "Matrix"),
         call. = FALSE, domain = NA)

Matrix.verbose <- function()
    getOption("Matrix.verbose", .MatrixEnv[["verbose"]])

Matrix.warn <- function()
    getOption("Matrix.warn", .MatrixEnv[["warn"]])

## MJ: or maybe a general one?
if(FALSE) {
Matrix.getOption <- function(x, default = .Matrix.Env[[x]])
    getOption(paste0("Matrix.", x), default)
} ## MJ

Matrix.message <- function(..., .M.level = 1, call. = FALSE, domain = NULL) {
    if(Matrix.verbose() >= .M.level) {
        m <-
            if((w <- Matrix.warn()) < 1)
                function(..., call., domain) message(..., domain = domain)
            else if(w < 2)
                warning
            else stop
        m(..., call. = call., domain = domain)
    }
}

## MJ: Implement forceTriangular() and export this and that?
## MJ: Notably this provides a model for (maybe, in the future) allowing
## forceSymmetric(<non-square>) ... by truncating the "too long" dimension.
forceDiagonal <- function(x, diag = NA_character_) {
    y <- diag(x, names = FALSE) # FIXME? don't allocate if diag == "U"
    cl <- switch(typeof(y),
                 logical =
                     {
                         if(is.na(diag))
                             diag <-
                             if(!is.na(a <- all(y        )) && a) "U" else "N"
                         if(.M.kind(x) == "n") "ndiMatrix" else "ldiMatrix"
                     },
                 integer =
                     {
                         if(is.na(diag))
                             diag <-
                             if(!is.na(a <- all(y == 1L  )) && a) "U" else "N"
                         storage.mode(y) <- "double"
                         "ddiMatrix"
                     },
              ## integer =
              ##     {
              ##         if(is.na(diag))
              ##             diag <-
              ##             if(!is.na(a <- all(y == 1L  )) && a) "U" else "N"
              ##         "idiMatrix"
              ##     },
                 double =
                     {
                         if(is.na(diag))
                             diag <-
                             if(!is.na(a <- all(y == 1   )) && a) "U" else "N"
                         "ddiMatrix"
                     },
                 complex =
                     stop(gettextf("complex %s not yet implemented",
                                   "diagonalMatrix"),
                          domain = NA),
              ## complex =
              ##     {
              ##         if(is.na(diag))
              ##             diag <-
              ##             if(!is.na(a <- all(y == 1+0i)) && a) "U" else "N"
              ##         "zdiMatrix"
              ##     },
                 stop(gettextf("cannot coerce matrix of type \"%s\" to %s",
                               typeof(y), "diagonalMatrix"),
                      domain = NA))
    n <- length(y)
    d <- dim(x)
    dn <- dimnames(x) %||% list(NULL, NULL)
    if(any(d > n)) {
        d <- c(n, n)
        w <- if(d[1L] > n) 1L else 2L
        if(!is.null(dnw <- dn[[w]]))
            dn[[w]] <- dnw[seq_len(n)]
    }
    new(cl, Dim = d, Dimnames = dn, diag = diag,
        x = if(diag == "N") y else y[FALSE])
}

.tCRT <- function(x, lazy = TRUE) .Call(R_sparse_transpose, x, lazy)

.drop0 <- function(x, tol = 0, isM = TRUE) {
    if(isM)
        return(.Call(R_sparse_drop0, x, tol))
    ## TODO: write sparseVector code in C and respecting 'tol'
    if(.M.kind(x) == "n")
        return(x)
    x.x <- x@x
    k <- which(is.na(x.x) | x.x)
    if(length(k)) {
        x@i <- x@i[k]
        x@x <- x.x[k]
    }
    x
}

drop0 <- function(x, tol = 0, is.Csparse = NA, give.Csparse = TRUE) {
    tryCoerce <-
        if(give.Csparse)
            is.na(is.Csparse) || !is.Csparse
        else if(.M.kind(x) != "n")
            !.isCRT(x)
        else if(all(.M.repr(x) != c("C", "R", "T", "i")))
            TRUE
        else return(x) # n[gst][CRT]Matrix or indMatrix
    if(tryCoerce)
        x <- if(isS4(x)) .M2C(x) else .m2sparse.checking(x, ".", "C")
    .Call(R_sparse_drop0, x, as.double(tol))
}

indDiag <- function(n, upper = TRUE, packed = FALSE)
    .Call(R_index_diagonal, n, packed, upper)

indTri <- function(n, upper = TRUE, diag = FALSE, packed = FALSE)
    .Call(R_index_triangle, n, packed, upper, diag)

## For sparseness handling, return a
## 2-column (i,j) matrix of 0-based indices of non-zero entries:

##' the workhorse for non0ind(.), but occasionally used directly
non0.i <- function(M, cM = class(M), uniqT = TRUE) {
    cld <- getClassDef(cM)
    if(extends(cld, "CsparseMatrix"))
        .Call(compressed_non_0_ij, M, TRUE)
    else if(extends(cld, "RsparseMatrix"))
        .Call(compressed_non_0_ij, M, FALSE)
    else if(extends(cld, "TsparseMatrix")) {
        if(uniqT && !isUniqueT(M))
	    .Call(compressed_non_0_ij, .M2C(M), TRUE)
	else cbind(M@i, M@j, deparse.level = 0L)
    } else if(extends(cld, "diagonalMatrix")) {
        i <- seq.int(from = 0L, length.out = M@Dim[1L])
        if(M@diag == "N")
	    i <- i[isN0(M@x)]
	cbind(i, i, deparse.level = 0L)
    } else if(extends(cld, "indMatrix")) {
        perm <- M@perm
        i <- seq.int(from = 0L, length.out = length(perm))
        if(M@margin == 1L)
            cbind(i, perm - 1L, deparse.level = 0L)
        else cbind(perm - 1L, i, deparse.level = 0L)
    } else
        stop(gettextf("non0.i() not yet implemented for class %s", dQuote(cM)),
             domain = NA)
}

##' the "more versatile / user" function (still not exported):
non0ind <- function(x, cld = getClassDef(class(x)),
		    uniqT = TRUE, xtendSymm = TRUE, check.Udiag = TRUE)
{
    if(is.numeric(x)) {
        n <- length(x)
        return(if(n == 0L)
                   integer(0L)
               else if(is.matrix(x))
                   arrayInd(seq_len(n)[isN0(x)], dim(x)) - 1L
               else (0:(n-1L))[isN0(x)])
    }
    stopifnot(extends(cld, "sparseMatrix"))

    ij <- non0.i(x, cld, uniqT = uniqT)
    if(xtendSymm && extends(cld, "symmetricMatrix")) {
        ## also get "other" triangle, but not the diagonal again
	notdiag <- ij[, 1L] != ij[, 2L]
	rbind(ij, ij[notdiag, 2:1], deparse.level = 0L)
    } else if(check.Udiag && extends(cld, "triangularMatrix") &&
              x@diag == "U") {
        i <- seq.int(from = 0L, length.out = x@Dim[1L])
        rbind(ij, cbind(i, i, deparse.level = 0L), deparse.level = 0L)
    } else ij
}

##' Decode "encoded" (i,j) indices back to  cbind(i,j)
##' This is the inverse of encodeInd(.)
##'
##' @title Decode "Encoded" (i,j) Indices
##' @param code integer in 0:((n x m - 1)  <==> encodeInd(.) result
##' @param nr the number of rows
##' @return
##' @author Martin Maechler
decodeInd <- function(code, nr)
    cbind(as.integer(code %% nr), as.integer(code %/% nr),
          deparse.level = 0L)

complementInd <- function(ij, dim, orig1=FALSE, checkBnds=FALSE) {
    ## Purpose: Compute the complement of the 2-column 0-based ij-matrix
    ##		but as 1-based indices
    n <- prod(dim)
    if(n == 0L)
        return(integer(0L))
    seq_len(n)[-(1L + .Call(m_encodeInd, ij, dim, orig1, checkBnds))]
}

WhichintersectInd <- function(ij1, ij2, di, orig1=FALSE, checkBnds=FALSE) {
    ## from 2-column (i,j) matrices where i \in {0,.., nrow-1},
    ## find *where*  common entries are in ij1 & ij2
    m1 <- match(.Call(m_encodeInd, ij1, di, orig1, checkBnds),
                .Call(m_encodeInd, ij2, di, orig1, checkBnds))
    ni <- !is.na(m1)
    list(which(ni), m1[ni])
}

anyDuplicatedT <- function(x, ...) {
    mn <- prod(d <- x@Dim)
    if(mn <= .Machine$integer.max)
        anyDuplicated.default(          x@j  * d[1L] + x@i, ...)
    else if(mn <= 0x1p+53)
        anyDuplicated.default(as.double(x@j) * d[1L] + x@i, ...)
    else anyDuplicated.default(.mapply(c, list(x@i, x@j), NULL), ...)
}

isUniqueT <- function(x, byrow = FALSE, isT = is(x, "TsparseMatrix"))
    isT && !is.unsorted(if(byrow) order(x@i, x@j) else order(x@j, x@i),
                        strictly = TRUE)

asUniqueT <- function(x, byrow = FALSE, isT = is(x, "TsparseMatrix"))
    if(isUniqueT(x, byrow, isT)) x else .M2T(if(byrow) .M2R(x) else .M2C(x))

aggregateT <- function(x) .Call(Tsparse_aggregate, x)

mat2triplet <- function(x, uniqT = FALSE) {
    T <- as(x, "TsparseMatrix")
    if(uniqT)
        T <- asUniqueT(T, isT = TRUE)
    if(is(T, "nsparseMatrix"))
         list(i = T@i + 1L, j = T@j + 1L)
    else list(i = T@i + 1L, j = T@j + 1L, x = T@x)
}

.validateCsparse <- function(x, sort.if.needed = FALSE) {
    if(sort.if.needed)
        .Call(CsparseMatrix_validate_maybe_sorting, x)
    else .Call(CsparseMatrix_validate, x)
}

dmperm <- function(x, nAns = 6L, seed = 0L) {
    ## NB: result is determined entirely by nonzero pattern of 'x'
    stopifnot(length(nAns <- as.integer(nAns)) == 1L, any(nAns == 2L * (1L:3L)),
              length(seed <- as.integer(seed)) == 1L, any(seed == -1L:1L))
    x <- if(isS4(x)) .M2gen(.M2C(x), "n") else .m2sparse(x, "ngC")
    .Call(Csparse_dmperm, x, nAns, seed)
}

## (matrix|denseMatrix)->denseMatrix as similar as possible to "target"
as_denseClass <- function(x, cl, cld = getClassDef(cl)) {
    kind <- .M.kind(x)
    cl <- .M.nonvirtual(new(cld))
    if(cl == "indMatrix")
        cl <- "ngeMatrix"
    symmetric <- substr(cl, 2L, 2L) == "s" && isSymmetric(x)
    triangular <- !symmetric &&
        substr(cl, 2L, 2L) == "t" && (it <- isTriangular(x))
    packed <- (symmetric || triangular) && substr(cl, 3L, 3L) == "p"
    if(isS4(x)) {
        r <- if(symmetric || triangular)
                 .M2kind(if(symmetric)
                             forceSymmetric(x)
                         else if(attr(it, "kind") == "U")
                             triu(x)
                         else tril(x),
                         kind)
             else .M2gen(x, kind)
        if(packed) pack(r) else r
    }
    else if(symmetric)
        .m2dense(x, paste0(kind, "s", if(packed) "p" else "y"), "U", NULL)
    else if(triangular)
        .m2dense(x, paste0(kind, "t", if(packed) "p" else "r"), attr(it, "kind"), "N")
    else .m2dense(x, paste0(kind, "ge"))
}

## (matrix|sparseMatrix)->CsparseMatrix as similar as possible to "target"
as_CspClass <- function(x, cl, cld = getClassDef(cl)) {
    kind <- .M.kind(x)
    cl <- .M.nonvirtual(new(cld))
    if(cl == "indMatrix")
        cl <- "ngeMatrix"
    symmetric <- substr(cl, 2L, 2L) == "s" && isSymmetric(x)
    triangular <- !symmetric &&
        substr(cl, 2L, 2L) == "t" && (it <- isTriangular(x))
    if(isS4(x)) {
        r <- if(symmetric || triangular)
                 .M2kind(if(symmetric)
                             forceSymmetric(x)
                         else if(attr(it, "kind") == "U")
                             triu(x)
                         else tril(x),
                         kind)
             else .M2gen(x, kind)
        .M2C(r)
    }
    else if(symmetric)
        .m2sparse(x, paste0(kind, "sC"), "U", NULL)
    else if(triangular)
        .m2sparse(x, paste0(kind, "tC"), attr(it, "kind"), "N")
    else .m2sparse(x, paste0(kind, "gC"))
}

## as(<[Mm]atrix>, <non-unit triangular CsparseMatrix>)
asCspN <- function(x) .Call(R_sparse_diag_U2N, as(x, "CsparseMatrix"))

diagU2N <- function (x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(extends(cl, "triangularMatrix") && x@diag == "U")
        .diagU2N(x, cl = cl, checkDense = checkDense)
    else x
}

.diagU2N <- function(x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(!checkDense && extends(cl, "denseMatrix"))
        x <- as(x, "CsparseMatrix")
    ..diagU2N(x)
}

..diagU2N <- function(x) {
    diag(x) <- TRUE
    x
}

diagN2U <- function(x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(extends(cl, "triangularMatrix") && x@diag == "N")
        .diagN2U(x, cl = cl, checkDense = checkDense)
    else x
}

.diagN2U <- function(x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(!checkDense & (isDense <- extends(cl, "denseMatrix")))
        ..diagN2U(as(x, "CsparseMatrix"), sparse = TRUE)
    else ..diagN2U(x, sparse = !isDense)
}

..diagN2U <- function(x, sparse) {
    if(sparse && x@Dim[1L] > 0L)
        x <- switch(x@uplo,
                    U = .Call(R_sparse_band, x, 1L, NULL),
                    L = .Call(R_sparse_band, x, NULL, -1L))
    x@diag <- "U"
    x
}

## Caches 'value' in the 'factors' slot of 'x' { NOT a copy of 'x' ! }
## and returns 'value'
.set.factor <- function(x, name, value, warn.no.slot = FALSE)
    .Call(R_set_factor, x, name, value, warn.no.slot)

##' Compute the three "parts" of two sets:
.setparts <- function(x, y) {
    n1 <- length(m1 <- match(x, y, 0L))
    n2 <- length(m2 <- match(y, x, 0L))
    ix <- which(m1 == 0L)
    iy <- which(m2 == 0L)
    list(x.only = x[ix], ix.only = ix, mx = m1,
         y.only = y[iy], iy.only = iy, my = m2,
         int = if(n1 < n2) y[m1] else x[m2])
}

##' *Only* to be used as function in
##'    setMethod.("Compare", ...., .Cmp.swap)  -->  ./Ops.R  & ./diagMatrix.R
.Cmp.swap <- function(e1, e2) {
    ## "swap RHS and LHS" and use the method below:
    switch(.Generic,
	   "==" =,
           "!=" = callGeneric(e2, e1),
	   "<"	= e2 >	e1,
	   "<=" = e2 >= e1,
	   ">"	= e2 <	e1,
	   ">=" = e2 <= e1)
}

.diag.dsC <- function(x, Chx = Cholesky(x, LDL = TRUE), res.kind = "diag") {
    if(!missing(Chx))
        stopifnot(.CHM.is.LDL(Chx), is.integer(Chx@p), is.double(Chx@x))
    .Call(tCsparse_diag, Chx, res.kind)
    ##    ^^^^^^^^^^^^^ in ../src/Csparse.c
}

dimScale <- function(x, d1 = sqrt(1/diag(x, names = FALSE)), d2 = d1) {
    dim.x <- dim(x)
    D1 <- Diagonal(n = dim.x[1L], x = d1)
    D2 <- if(missing(d2)) D1 else Diagonal(n = dim.x[2L], x = d2)
    y <- D1 %*% x %*% D2 # inefficient for symmetricMatrix 'x', but "general"
    if(isS4(x) && is(x, "symmetricMatrix") && identical(d1, d2))
        y <- forceSymmetric(y, x@uplo)
    if(is.list(dn <- dimnames(x)))
        y@Dimnames <- dn
    y
}

rowScale <- function(x, d) {
    y <- Diagonal(n = nrow(x), x = d) %*% x
    if(is.list(dn <- dimnames(x)))
        y@Dimnames <- dn
    y
}

colScale <- function(x, d) {
    y <- x %*% Diagonal(n = ncol(x), x = d)
    if(is.list(dn <- dimnames(x)))
        y@Dimnames <- dn
    y
}
