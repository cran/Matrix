####  Utilities  for  Sparse Model Matrices

## > from <- factor(sample(c(1:100, NA), size = 1e+04, replace = TRUE,
## +                       prob = rep.int(c(1, 20), c(100L, 1L))))
## > microbenchmark::microbenchmark(
## +                     Matrix:::fac2sparse1(from, drop.unused.levels = FALSE),
## +                     Matrix:::fac2sparse2(from, drop.unused.levels = FALSE),
## +                     Matrix:::fac2sparse (from, drop.unused.levels = FALSE),
## +                     times = 10000L)
## Unit: microseconds
##                                                    expr     min      lq     mean  median      uq       max neval
##  Matrix:::fac2sparse1(from, drop.unused.levels = FALSE)  94.997 105.001 118.6117 107.789 110.413  4439.111 10000
##  Matrix:::fac2sparse2(from, drop.unused.levels = FALSE) 538.289 586.341 643.7947 595.197 602.700 43128.556 10000
##  Matrix:::fac2sparse (from, drop.unused.levels = FALSE) 163.139 185.115 221.7241 190.117 193.889 44887.784 10000

if(FALSE) {
## The "first" version, no longer used:
fac2sparse1 <- function(from,
                        to = c("d", "i", "l", "n", "z"),
                        drop.unused.levels = FALSE)
{
    ## factor(-like) --> sparseMatrix {also works for integer, character}
    fact <- if (drop.unused.levels) factor(from) else as.factor(from)
    levs <- levels(fact)
    n <- length(fact)
    to <- match.arg(to)
    ## MM: using new() and then assigning slots has efficiency "advantage"
    ##     of *not* validity checking
    res <- new(paste0(to, "gCMatrix"))
    res@i <- as.integer(fact) - 1L # 0-based
    res@p <- 0:n
    res@Dim <- c(length(levs), n)
    res@Dimnames <- list(levs, NULL)
    if(to != "n")
	res@x <- rep.int(switch(to,
				"d" = 1., "i" = 1L, "l" = TRUE, "z" = 1+0i),
			 n)
    res
}

## The "second" version, no longer used:
## * this one handles NAs correctly but may be slightly less efficient
##   than fac2sparse1() above
fac2sparse2 <- function(from,
                        to = c("d", "i", "l", "n", "z"),
                        drop.unused.levels = TRUE,
                        repr = c("C", "T", "R"),
                        giveCsparse)
{
    ## factor(-like) --> sparseMatrix {also works for integer, character}
    fact <- if (drop.unused.levels) factor(from) else as.factor(from)
    levs <- levels(fact)
    n <- length(fact)
    to <- match.arg(to)
    i <- as.integer(fact) - 1L ## 0-based indices
    df <- data.frame(i = i, j = if(n) 0:(n-1L) else integer())[!is.na(i),]
    if(to != "n")
	df$x <- rep.int(switch(to,
			       "d" = 1., "i" = 1L, "l" = TRUE, "z" = 1+0i),
			nrow(df))
    T <- do.call(new, c(list(Class = paste0(to, "gTMatrix"),
                             Dim = c(length(levs), n),
                             Dimnames = list(levs, names(fact))), df))
    ## silent, back compatible (not yet warning about 'giveCsparse' deprecation):
    repr <- if(missing(repr) && !missing(giveCsparse))
		if(giveCsparse) "C" else "T"
	    else match.arg(repr)
    switch(repr,
	   "C" = .T2C(T),
	   "T" =    T,# TsparseMatrix
	   "R" = .T2R(T))
}
} # if(FALSE)

## The "third" version, dealing with NAs _and_ fast:
fac2sparse <- function(from,
                       to = c("d", "l", "n"), # no "i" or "z" yet
		       drop.unused.levels = TRUE,
                       repr = c("C", "R", "T"),
                       giveCsparse)
{
    cl <- ".g.Matrix"
    substr(cl, 1L, 1L) <- to <- match.arg(to)
    substr(cl, 3L, 3L) <- repr <-
        if(!missing(repr) || missing(giveCsparse))
            match.arg(repr)
        else if(giveCsparse)
            "C"
        else "T"

    ## Arguments are valid: _now_ allocate
    from <-
        if(drop.unused.levels && is.factor(from))
            factor(from)
        else as.factor(from)
    n <- length(from)
    lv <- levels(from)
    nlv <- length(lv)

    res <- new(cl)
    res@Dim <- c(nlv, n)
    res@Dimnames <- list(lv, names(from))
    i. <-
        switch(repr,
               "C" =
                   {
                       not.na <- !is.na(from)
                       res@p <- c(0L, cumsum(not.na))
                       res@i <- as.integer(from)[not.na] - 1L
                   },
               "R" =
                   {
                       res@p <- c(0L, cumsum(tabulate(from, nlv)))
                       res@j <- order(from, na.last = NA) - 1L
                   },
               "T" =
                   {
                       which.not.na <- which(!is.na(from))
                       res@i <- as.integer(from)[which.not.na] - 1L
                       res@j <- which.not.na - 1L
                   }
               )
    if(to != "n")
        res@x <- rep.int(switch(to, "d" = 1, "l" = TRUE, "i" = 1L, "z" = 1+0i),
                         length(i.))
    res
}

setAs("factor", "sparseMatrix", function(from) fac2sparse(from, to = "d"))
## FIXME? support factor->[dlnCRT]sparseMatrix?

##' fac2Sparse() := fac2sparse w/ contrasts
##'
##' @param from factor of which we want the "contrasted" (indicator)
##'   design matrix
##' @param to character string specifying the response type
##' @param drop.unused.level logical indicating if non-present factor
##'   levels should be dropped, via  factor(from)
##' @param factorPatt12 logical vector fp[] of length 2
##'   fp[1] : give contrasted t(X);  fp[2] : give "dummy" t(X) [=fac2sparse()]
##' @param contrasts.arg character string or NULL or (coercible to)
##'		sparseMatrix, specifying the contrast
##'
##' @return a list of length two, each with the corresponding t(model matrix),
##'	when the corresponding factorPatt12 is true.
fac2Sparse <- function(from,
                       to = c("d", "l", "n"), # no "i" or "z" yet
		       drop.unused.levels = TRUE,
                       repr = c("C", "R", "T"),
                       giveCsparse,
		       factorPatt12,
                       contrasts.arg = NULL)
{
    stopifnot(is.logical(factorPatt12), length(factorPatt12) == 2L)
    if(!any(factorPatt12))
        return(list(NULL, NULL)) # nothing to do

    m <- fac2sparse(from, to = to, drop.unused.levels = drop.unused.levels,
                    repr = repr, giveCsparse = giveCsparse)
    list(if(factorPatt12[1L]) {
             ## contrasted coding, using 'contrasts.arg'
             if(is.null(contrasts.arg))
                 contrasts.arg <- getOption("contrasts")[[if(is.ordered(from))
                                                              "ordered"
                                                          else "unordered"]]
             crossprod(if(is.character(contrasts.arg)) {
                           ## calling  contr.*() with level names directly:
                           contr <- get(contrasts.arg, mode = "function")
                           contr(m@Dimnames[[1L]], sparse = TRUE)
                       } else as(contrasts.arg, "sparseMatrix"),
                       m)
         },
         if(factorPatt12[2L])
             ## uncontrasted ("dummy") coding
             m
         )
}

## Cut and paste from stats:::deparse2() in stats/R/models.R
deparse2 <- function(x)
    paste(deparse(x, width.cutoff = 500L,
                  backtick = !is.symbol(x) && is.language(x)),
          collapse = " ")

## Cut and paste from stats:::model.matrix.default() in stats/R/models.R,
## with some adaptation, most notably at the very end where we do _not_
## call the C-level utility of 'stats'
sparse.model.matrix <- function(object,
                                data = environment(object),
                                contrasts.arg = NULL,
                                xlev = NULL,
                                transpose = FALSE,
                                drop.unused.levels = FALSE,
                                row.names = TRUE,
                                sep = "",
                                verbose = FALSE,
                                ...)
{
    t <- if(missing(data)) terms(object) else terms(object, data=data)
    if (is.null(attr(data, "terms")))
	data <- model.frame(object, data, xlev=xlev)
    else {
	reorder <- match(vapply(attr(t, "variables"), deparse2, "")[-1L],
                         names(data))
	if (anyNA(reorder))
	    stop("model frame and formula mismatch in sparse.model.matrix()")
	if(!identical(reorder, seq_len(ncol(data))))
	    data <- data[,reorder, drop=FALSE]
    }
    int <- attr(t, "response")
    if(length(data)) {
        contr.funs <- as.character(getOption("contrasts"))
        namD <- names(data)
        ## turn any character columns into factors
        for(i in namD)
            if(is.character(data[[i]]))
                data[[i]] <- factor(data[[i]])
        isF <- vapply(data, function(x) is.factor(x) || is.logical(x), NA)
        isF[int] <- FALSE
        isOF <- vapply(data, is.ordered, NA)
        for(nn in namD[isF])            # drop response
            if(is.null(attr(data[[nn]], "contrasts")))
                contrasts(data[[nn]]) <- contr.funs[1 + isOF[nn]]
        ## it might be safer to have numerical contrasts:
        ##	  get(contr.funs[1 + isOF[nn]])(nlevels(data[[nn]]))
        if (!is.null(contrasts.arg)) {
            if (!is.list(contrasts.arg))
                warning("non-list contrasts argument ignored")
            else {  ## contrasts.arg is a list
                if (is.null(namC <- names(contrasts.arg)))
                    stop("'contrasts.arg' argument must be named")
                for (nn in namC) {
                    if (is.na(ni <- match(nn, namD)))
                        warning(gettextf("variable '%s' is absent, its contrast will be ignored", nn),
                                domain = NA)
                    else {
                        ca <- contrasts.arg[[nn]]
                        ## contrasts(*, ncol(m)) <- m works also
                        ## for function||character 'm' in R >= 4.2,
                        ## which supports how.many=NULL
                        if(is.matrix(ca) || is(ca, "Matrix"))
                            contrasts(data[[ni]], ncol(ca)) <- ca
                        else contrasts(data[[ni]]) <- ca
                    }
                }
            }
        } ## non-null contrasts.arg
    } else { #  no rhs terms ('~1', or '~0'): internal model.matrix needs some variable
	isF <- FALSE
	data[["x"]] <- raw(nrow(data))
    }

    ## <stats>
    ## ans <- .External2(C_modelmatrix, t, data)
    ## if(any(isF))
    ##     attr(ans, "contrasts") <- lapply(data[isF], attr, "contrasts")
    ## </stats>
    ## <Matrix>
    if(verbose) {
	cat("model.spmatrix(t, data, ...) with t =\n")
        str(t, give.attr = FALSE)
    }
    ans <- model.spmatrix(trms = t,
                          mf = data,
                          transpose = transpose,
                          drop.unused.levels = drop.unused.levels,
			  row.names = row.names,
                          sep = sep,
                          verbose = verbose)
    ## MJ: hmm ... our tests require that this "slot" exists,
    ##     even in the empty case, i.e., !any(isF) ... why?
    attr(ans, "contrasts") <- lapply(data[isF], attr, "contrasts")
    ## </Matrix>
    ans
} # sparse.model.matrix

##' Produce the t(Z); Z = "design matrix" of (X : Y), where
##'             --- t(Z) : aka rowwise -version : "r"
##'
##' @title sparse model matrix for 2-way interaction
##' @param X and Y either are numeric matrices {maybe 1-column}
##' @param Y       or "as(<factor>, sparseM)" {dgCMatrix}
##' @param do.names logical
##' @param forceSparse logical
##' @return
##' @author Martin Maechler
.sparse.interaction.2 <- # formerly sparse2int()
function(X, Y, do.names = TRUE, forceSparse = FALSE, verbose = FALSE)
{
### FIXME -- the    X[rep(..), ] * Y[rep(..), ]   construct can become HUGE, even for sparse X[],Y[]
### ----- --> Matrix bug #1330 and  ~/R/MM/Pkg-ex/Matrix/sparse-matrix-fix.R

    ## MJ: Moreover, as(<factor>, "sparseMatrix") is a dgCMatrix,
    ##     for which row-indexing is already rather inefficient ... FIXME?

    sx <- isS4(X)
    sy <- isS4(Y)
    nx <- (dx <- dim(X))[1L]
    ny <- (dy <- dim(Y))[1L]
    if(verbose)
        cat(sprintf(".sparse.interaction.2(%s[%d], %s[%d])\n",
                    if(sx) "<sparse>" else "<N>", nx,
                    if(sy) "<sparse>" else "<N>", ny))
    if(do.names) {
	dnx <- dimnames(X)
	dny <- dimnames(Y)
    }
    dimnames(X) <- dimnames(Y) <- list(NULL, NULL)

    r <-
        if(sx && sy) {

            ## 'X' and 'Y' are dgCMatrix
            (if(ny == 1L) X else X[rep.int(seq_len(nx), times = ny), ]) *
            (if(nx == 1L) Y else Y[rep    (seq_len(ny),  each = nx), ])

        } else if (sx) {

            ## 'X' is a dgCMatrix, 'Y' is a numeric matrix
            if(ny <= 1L) {
                ## FIXME: a similar trick would be applicable for ny > 1
                r <- X
                dp <- X@p[-1L] - X@p[-(dx[2L]+1L)]
                ## stopifnot(all(dp %in% 0:1))
                r@x <- Y[dp == 1L] * X@x
                r
            } else {
                                     X[rep.int(seq_len(nx), times = ny), ] *
                (if(nx == 1L) Y else Y[rep    (seq_len(ny),  each = nx), ])
            }

        } else if(sy) {

            ## 'X' is a numeric matrix, 'Y' is a dgCMatrix
            if(nx <= 1L) {
                ## FIXME: a similar trick would be applicable for nx > 1
                r <- Y
                dp <- Y@p[-1L] - Y@p[-(dy[2L]+1L)]
                ## stopifnot(all(dp %in% 0:1))
                r@x <- X[dp == 1L] * Y@x
                r
            } else {
                (if(ny == 1L) X else X[rep.int(seq_len(nx), times = ny), ]) *
                                     Y[rep    (seq_len(ny),  each = nx), ]
            }

        } else {

            ## 'X' and 'Y' are numeric matrices
            r <- (if(ny == 1L) X else X[rep.int(seq_len(nx), times = ny), ]) *
                 (if(nx == 1L) Y else Y[rep    (seq_len(ny),  each = nx), ])
            if(forceSparse)
                .Call(R_dense_as_sparse, r, "dgC", NULL, NULL)
            else r

        }

    ## FIXME: This 'names' business needs a good solution ...
    ##        but maybe "up in the caller" ...
    if(do.names &&
       !is.null(dim(r)) &&
       !is.null(rnx <- dnx[[1L]]) &&
       !is.null(rny <- dny[[1L]]))
        dimnames(r)[[1L]] <- outer(rnx, rny, paste, sep = ":")
    r
} # .sparse.interaction.2

##' Sparse Model Matrix for a (high order) interaction term  A:B:x:C
##'
##' @param rList list(.) of (transposed) single-factor model matrices,
##'	belonging to, say, factors  a, b, c,...
##' @param do.names
##' @param forceSparse
##' @param verbose
##' @return the model matrix corresponding to a:b:...
.sparse.interaction.N <- # formerly sparseInt.r()
function(rList, do.names = TRUE, forceSparse = FALSE, verbose = FALSE)
{
    if((n <- length(rList)) == 0L)
        return(NULL) # caller beware
    if(verbose)
	cat(sprintf(".sparse.interaction.N(<list>[%d], fS=%s): is.mat=(%s)\n",
                    n, forceSparse, paste0(symnum(vapply(rList, is.matrix, NA)),
                                           collapse = "")),
            sep = "")
    r <- rList[[1L]]
    if(n > 1L)
        for(i in 2:n)
	    r <- .sparse.interaction.2(r, rList[[i]],
                                       forceSparse = forceSparse,
                                       do.names = do.names,
                                       verbose = verbose)
    if(forceSparse && (is.matrix(r) || is(r, "denseMatrix")))
        .Call(R_dense_as_sparse, r, "dgC", NULL, NULL)
    else r
} # .sparse.interaction.N

## MJ: unused
if(FALSE) {
is.model.frame <- function(x)
{
  ## Purpose: check if x is a "valid" model.frame
  ## ------------------------------------------------------------
  ## Author: Martin Maechler, Date: 3 Jul 2009
    is.data.frame(x) &&
    !is.null(tms <- attr(x, "terms")) &&
    inherits(tms, "terms") && ## is.terms() would be better
    inherits(tms, "formula") &&
    is.matrix(attr(tms, "factors")) &&
    is.language(vv <- attr(tms, "variables")) &&
    vv[[1]] == as.symbol("list") &&
    all(vapply(as.list(vv[-1]), as.character, "") %in% colnames(x))
    ## all((vars <- sapply(as.list(vv[-1]), as.character)) %in% colnames(x))
    ## and we could go on testing vars
}
} ## MJ

##' Create a sparse model matrix from a model frame.
##'
##' @title Sparse Model Matrix from Model Frame
##' @param trms a "terms" object
##' @param mf a data frame, typically resulting from  model.frame()
##' @param transpose logical indicating if  X' = t(X) {is faster!}
##'	or X should be returned
##' @param drop.unused.levels logical indicating if unused factor
##'	levels should be dropped
##' @param row.names
##' @return sparse matrix (class "dgCMatrix")
##' @author Martin Maechler
model.spmatrix <- function(trms, mf, transpose=FALSE,
			   drop.unused.levels = FALSE, row.names=TRUE, sep="",
                           verbose=FALSE)
{
    ## Author: Martin Maechler, Date:  7 Jul 2009

    ## mf is a model frame or a "simple" data.frame [after reorder !]
    stopifnot(is.data.frame(mf))
    n <- nrow(mf)
    if(row.names)
	rnames <- row.names(mf)
    ## mf:  make into list, dropping all attributes (but the names)
### FIXME: for poly(., 5)  mf has a 5-column matrix as "one column" => looses names here
    fnames <- names(mf <- unclass(mf))
    attributes(mf) <- list(names = fnames)

    if(length(factorPattern <- attr(trms, "factors"))) {
	d <- dim(factorPattern)
	nVar <- d[1]
	nTrm <- d[2]
	n.fP <- dimnames(factorPattern)
	fnames <- n.fP[[1]] # == names of variables {incl. "F(var)"} in the model
	Names  <- n.fP[[2]] # == colnames == names of terms:  "a", "b:c", ...
    } else { ## degenerate, e.g.  'Y ~ 1'
	nVar <- nTrm <- 0L
	fnames <- Names <- character(0)
    }
    ## all the "variables in the model" are also in "mf", including "sin(x)";
    ## actually, ..../src/main/model.c even assumes
    stopifnot((m <- length(mf)) >= nVar)
    if(verbose)
	cat(sprintf("model.spmatrix(): (n=%d, nVar=%d (m=%d), nTrm=%d)\n",
		    n, nVar,m, nTrm))
    if(m > nVar) mf <- mf[seq_len(nVar)]
    stopifnot(fnames == names(mf), allow.logical0 = TRUE)
    noVar <- nVar == 0
    ##>> this seems wrong; we use  1:nVar for indexing mf[] below ..
    ##>> if(noVar) nVar <- 1L # (as in ~/R/D/r-devel/R/src/main/model.c)
    ## Note: "character" variables have been changed to factor in the caller;
    ##     hence: both factor and *logical*  should be dealt as factor :
    is.f <- if(noVar) logical(0) else vapply(mf, function(.)
					     is.factor(.) | is.logical(.), NA)
    indF <- which(is.f)
    if(verbose) { cat(" --> indF =\n"); print(indF) }
    hasInt <- attr(trms, "intercept") == 1
    ## the degree of interaction:
    ## intOrder <- attr(trms, "order")
    ##
    if(!hasInt && length(indF)) {
	## change the '1' of the first factor into a '2' :
	if(any(i1 <- factorPattern[indF, ] == 1))
	    ## replace at the first '1' location:
	    factorPattern[indF,][which.max(i1)] <- 2L
	else {}
	## nothing to do
    }
    ## Convert "factors" to "Rowwise- sparseMatrix ("dummy"-matrix) -----------
    ## Result: a list of sparse model matrices for the "factor"s :
    f.matr <- structure(vector("list", length = length(indF)),
			names = fnames[indF])
    i.f <- 0
    ## ---- For each variable in the model -------------------
    for(i in seq_len(nVar)) {
	nam <- fnames[i]
	f <- mf[[i]]
	if(is.f[i]) {
	    fp <- factorPattern[i,] ## == factorPattern[nam,]
	    contr <- attr(f, "contrasts")
	    f.matr[[(i.f <- i.f + 1)]] <- # a list of 2
		lapply(fac2Sparse(f, to = "d",
				  drop.unused.levels=drop.unused.levels,
				  factorPatt12 = 1:2 %in% fp,
				  contrasts.arg = contr),
		       function(s) {
                           if(!is.null(s)) {
                               ## for some contr.*(), have lost rownames; hmm..
                               if(is.null(rn <- rownames(s)))
                                   rn <- seq_len(nrow(s))
                               rownames(s) <- paste(nam, rn, sep = sep)
			   }
                           s
		       })
	} else { ## continuous variable --> "matrix" - for all of them
	    if(any(iA <- (cl <- class(f)) == "AsIs")) # drop "AsIs" class
		class(f) <- if(length(cl) > 1L) cl[!iA]
	    nr <- if(is.matrix(f)) nrow(f <- t(f)) else (dim(f) <- c(1L, length(f)))[1]
	    if(is.null(rownames(f)))
		rownames(f) <- if(nr == 1) nam else paste(nam, seq_len(nr), sep=sep)
	    mf[[i]] <- f
	}
    }
    if(verbose) {
	cat(" ---> f.matr list :\n")
	str(f.matr, max.level = as.integer(verbose))
	fNms <- format(dQuote(Names))
	dim.string <- gsub('5', as.character(floor(1+log10(n))),
			   " -- concatenating (r, rj): dim = (%5d,%5d) | (%5d,%5d)\n")
    }

    ## FIXME: do all this in C --

    getR <- function(N)			# using 'nm'
	if(!is.null(r <- f.matr[[N]])) r[[factorPattern[N, nm]]] else mf[[N]]
    vNms <- "(Intercept)"[hasInt]
    counts <- integer(nTrm)
    r <-
	if(hasInt) ## column of 1's - as sparse
	    new("dgCMatrix", i = 0:(n-1L), p = c(0L, n),
		Dim = c(n, 1L), x = rep.int(1, n))
	else new("dgCMatrix", Dim = c(n, 0L))
    if(transpose) r <- t(r)
    iTrm <- seq_len(nTrm)
    for(j in iTrm) { ## j-th term
	nm <- Names[j]
	if(verbose) cat(sprintf("term[%2d] %s .. ", j, fNms[j]))
	nmSplits <- strsplit(nm, ":", fixed=TRUE)[[1]]
	## NOTA BENE: This can be very slow when many terms are involved
	## FIXME ??? why does it use *much* memory in those cases ??
	rj <- .sparse.interaction.N(lapply(nmSplits, getR),
                                    do.names = TRUE,
                                    forceSparse = TRUE,
                                    verbose = verbose) # or just (verbose >= 2)
	if(verbose) cat(sprintf(dim.string, nrow(r), ncol(r), nrow(rj),ncol(rj)))
	## fast version of cbind2() / rbind2(), w/o checks, dimnames, etc
	r <- if(transpose) .Call(Csparse_vertcat, r, rj)
		else	   .Call(Csparse_horzcat, r, t(rj))
	## if(verbose) cat(" [Ok]\n")
	vNms <- c(vNms, dimnames(rj)[[1]])
	counts[j] <- nrow(rj)
    }
    rns <- if(row.names) rnames
    dimnames(r) <- if(transpose) list(vNms, rns) else list(rns, vNms)
    attr(r, "assign") <- c(if(hasInt) 0L, rep(iTrm, counts))
    r
} ## model.spmatrix()
