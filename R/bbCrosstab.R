bbCrosstab <- function(flist) {
    ind <- function(i,j) ((i-1) * i)/2 + j # index in lower triangle of matrix
    flist <- lapply(as(flist, "list"), function (x) as(x, "factor")[drop = TRUE])
    nfac <- length(flist)
    if (nfac < 1) return(new("bbSparseSy", x = list(), uplo = 'L'))
    nlev <- unlist(lapply(flist, function(x) length(levels(x))))
    nobs <- length(flist[[1]])
    if (any(lapply(flist, length) != nobs))
        stop("all factors in flist must have the same length")
    ones <- rep(1, nobs)
    lst <- vector("list", choose(nfac + 1, 2))
    if (!is.null(nms <- names(flist))) {
        nmat <- outer(nms, nms, FUN = "paste", sep = ":")
        names(lst) <- nmat[lower.tri(nmat, diag = TRUE)]
    }
    zb <- lapply(flist, function(x) as.integer(x) - 1:1) # zero-based indices
    for (i in seq(along = flist)) {
        for (j in 1:i) {
            lst[[ind(i,j)]] <- as(as(new("tripletMatrix",
                                       i = zb[[i]], j = zb[[j]],
                                       x = ones, Dim = nlev[c(i,j)]),
                                   "cscMatrix"),
                                "cscBlocked")
        }
    }
    new("bbCrosstab", x = lst, uplo = "L")
}

setMethod("isNested", signature(object = "bbCrosstab"),
          function(object, ...)
          all(unlist(lapply(object@x, function(x) any(diff(x@p) > 1)))))
