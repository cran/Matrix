bCrosstab <- function(flist, reorder = TRUE) {
    ind <- function(i,j) ((i-1) * i)/2 + j # index in rowwise lower triangle
    ## Coerce flist to a list of factors with no unused levels
    flist <- lapply(as.list(flist), function (x) as.factor(x)[drop = TRUE])

    nfac <- length(flist)
    if (nfac < 1) stop("flist must be a non-empty list of factors")
                                        # Check for consistent lengths
    nobs <- length(flist[[1]])
    if (any(lapply(flist, length) != nobs))
        stop("all factors in flist must have the same length")
    ones <- rep(1, nobs)
    nlev <- unlist(lapply(flist, function(x) length(levels(x))))
    ord <- rev(order(nlev))
    if (reorder && any(ord != seq(along = ord))) {
        nlev <- nlev[ord]
        flist <- flist[ord]
    }
                                          # establish output list and names
    lst <- vector("list", choose(nfac + 1, 2))
    if (is.null(nms <- names(flist))) nms <- paste("V", 1:nfac, sep = "")
    nmat <- outer(nms, nms, FUN = "paste", sep = ":")
    nms <- vector("character", length(lst))
    zb <- lapply(flist, function(x) as.integer(x) - 1:1) # zero-based indices
    for (j in 1:nfac) {
        for (i in j:nfac) {
            IND <- ind(i,j)
            lst[[IND]] <- as(new("dgTMatrix",
                                 i = zb[[i]], j = zb[[j]],
                                 x = ones, Dim = nlev[c(i,j)]),
                             ifelse(j != i,"dgCMatrix","dsCMatrix"))
            nms[ind(i,j)] <- nmat[i, j]
        }
    }
    names(lst) <- nms
    list(flist = flist, ctab = lst)
}

ctab2bblist <- function(ctab, nf, nc)
{
    ind <- function(i,j) ((i-1) * i)/2 + j # index in compressed lower triangle
    ans <- vector("list", length(ctab))
    names(ans) <- names(ctab)
    for (j in 1:nf) {
        for (i in j:nf) {
            IND <- ind(i, j)
            ctij <- ctab[[IND]]
            ans[[IND]] <- new("dgBCMatrix", p = ctij@p, i = ctij@p,
                              x = array(0, dim = c(nc[i], nc[j], length(ctij@x))))
        }
    }
    ans
}

Lstruct <- function(bcr, nc = rep(1, nf)) {
    ind <- function(i,j) ((i-1) * i)/2 + j # index in compressed lower triangle
    ctab <- bcr$ctab
    fl <- bcr$flist
    nf <- length(fl)
    Linv <- vector("list", length(fl))
    names(Linv) <- names(fl)

    if (all(unlist(lapply(ctab, function(x) all(diff(x@p)== 1))))) { # nested
        ZtZ <- ctab2bblist(ctab, nf, nc)
        for (j in 1:nf) Linv[[j]] <- ZtZ[[ind(j, j)]]
        return(list(flist = fl, ZtZ = ZtZ, Lmat = ZtZ, Linv = Linv))
    }
    ## non-nested case - here things get interesting
    ct11 <- ctab[[1]]
    Linv[[1]] <- new("dgBCMatrix", p = ct11@p, i = ct11@i,
                     x = array(0, dim = c(nc[1], nc[1], length(ct11@x))))
    for (i in 1:(nf - 1)) {
        ip1 <- i + 1
        res <- .Call("bCrosstab_project", ctab, i)
        fac <- fl[[ip1]]
        fl[[ip1]] <- factor(as.character(fac), levels = levels(fac)[1+res$perm])
        ctab <- .Call("bCrosstab_permute", res$ctab, i, res$perm)
        rLi <- res$Linv
        Linv[[ip1]] <- new("dgBCMatrix", p = rLi@p, i = rLi@i,
                           x = array(0, dim = c(nc[ip1], nc[ip1], length(rLi@x))))
    }
    list(flist = fl, ZtZ = ctab2bblist(bCrosstab(fl)$ctab, nf, nc),
         Lmat = ctab2bblist(ctab, nf, nc), Linv = Linv)
    ## FIXME: You probably don't want to define Lmat from ctab because
    ## it has off-diagonal elements generated from the inverses of the diagonals
}
