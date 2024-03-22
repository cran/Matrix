## METHODS FOR GENERIC: all.equal
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.V.a.e <- function(target, current, ...) {
    if((l1 <- length(target)) != (l2 <- length(current)))
        return(paste0("length(target) is ", l1, ", length(current) is ", l2))
    if(is.integer(l1) || l1 <= .Machine$integer.max) {
        i1 <- as.integer( target@i)
        i2 <- as.integer(current@i)
    } else {
        i1 <- trunc( target@i)
        i2 <- trunc(current@i)
    }
    x1 <- if(.hasSlot( target, "x"))  target@x else rep.int(TRUE, length(i1))
    x2 <- if(.hasSlot(current, "x")) current@x else rep.int(TRUE, length(i2))
    if(!identical(i1, i2)) {
        i3 <- sort.int(unique.default(c(i1, i2)))
        x1 <- replace(vector(typeof(x1), length(i3)), match(i1, i3, 0L), x1)
        x2 <- replace(vector(typeof(x2), length(i3)), match(i2, i3, 0L), x2)
    }
    all.equal(x1, x2, ...)
}

.M.attributes <-
function(x, exclude.informal, exclude.factors) {
    a <- attributes(x)
    if(isS4(x) && exclude.informal)
        a <- a[.slotNames(a)]
    if(length(a) == 0L)
        return(NULL)
    exclude <-
    if(!isS4(x))
        c("class", "dim", "dimnames")
    else if(.isMatrix(x))
        c("class", "Dim", "Dimnames",
          switch(.M.repr(x),
                 "C" = c("p", "i", if(.M.kind(x) != "n") "x"),
                 "R" = c("p", "j", if(.M.kind(x) != "n") "x"),
                 "T" = c("i", "j", if(.M.kind(x) != "n") "x"),
                 "d" = , "u" = , "p" = "x",
                 "i" = "perm"),
          switch(.M.shape(x),
                 "g" = if(exclude.factors) "factors",
                 "s" = c("uplo", if(exclude.factors) "factors"),
                 "t" = c("uplo", "diag"),
                 "d" = "diag",
                 "i" = if(exclude.factors) "factors"))
    else "class"
    nms <- names(a)
    i <- match(nms, exclude, 0L) == 0L
    if(any(i)) a[sort.int(nms[i])] else NULL
}

.M.attr.all.equal <-
function(target, current,
         check.type, check.class, check.attributes, check.factors, ...) {
    msg <- msg. <- NULL
    if(check.type  && !identical(t1 <- typeof(target), t2 <- typeof(current)))
        msg <- c(msg, paste0("typeof(target) is ", deparse(t1), ", typeof(current) is ", deparse(t2)))
    if(check.class && !identical(c1 <-  class(target), c2 <-  class(current)))
        msg <- c(msg, paste0( "class(target) is ", deparse(c1),  ", class(current) is ", deparse(c2)))
    if(is.na(check.attributes) || check.attributes) {
        if(!isTRUE(ae <- all.equal.raw(dim(target), dim(current), ...)))
            msg <- c(msg, paste0("dim: < ", ae, " >"))
        if(!isTRUE(ae <- all.equal.list(dimnames( target) %||% list(NULL, NULL),
                                        dimnames(current) %||% list(NULL, NULL),
                                        ...)))
            msg <- c(msg, paste0("dimnames: < ", ae, " >"))
        a1 <- .M.attributes( target, is.na(check.attributes), !check.factors)
        a2 <- .M.attributes(current, is.na(check.attributes), !check.factors)
        if(!((is.null(a1) && is.null(a2)) ||
             isTRUE(ae <- all.equal.list(a1, a2, ...))))
            msg <- msg. <- c(msg, paste0("Attributes: < ", ae, " >"))
    }
    list(msg, is.null(msg) != is.null(msg.))
}

.M.all.equal <-
function(target, current,
         check.type = check.class,
         check.class = TRUE,
         check.attributes = TRUE,
         check.factors = FALSE, ...) {
    msg <- .M.attr.all.equal(target, current,
                             check.type = check.type,
                             check.class = check.class,
                             check.attributes = check.attributes,
                             check.factors = check.factors, ...)
    if(!msg[[2L]]) {
        ae <-
        if(.isVector( target) || .isSparse( target) ||
           .isVector(current) || .isSparse(current)) {
            v1 <- as( target, "sparseVector")
            v2 <- as(current, "sparseVector")
            ae <- .V.a.e(v1, v2, ...)
        } else {
            v1 <- as( target, "vector")
            v2 <- as(current, "vector")
            ae <- all.equal(v1, v2, ...)
        }
        if(!isTRUE(ae))
            return(c(msg[[1L]], ae))
    }
    if(is.null(msg[[1L]])) TRUE else msg[[1L]]
}

setMethod("all.equal", c(target = "Matrix", current = "vector"),
          .M.all.equal)

setMethod("all.equal", c(target = "vector", current = "Matrix"),
          .M.all.equal)

setMethod("all.equal", c(target = "Matrix", current = "Matrix"),
          .M.all.equal)

## And for completeness:

setMethod("all.equal", c(target = "Matrix", current = "sparseVector"),
          .M.all.equal)

.V.attributes <-
function(x, exclude.informal) {
    a <- attributes(x)
    if(isS4(x) && exclude.informal)
        a <- a[.slotNames(a)]
    if(length(a) == 0L)
        return(NULL)
    exclude <-
    if(.isVector(x))
        c("class", "length", "i", if(.M.kind(x) != "n") "x")
    else "class"
    nms <- names(a)
    i <- match(nms, exclude, 0L) == 0L
    if(any(i)) a[sort.int(nms[i])] else NULL
}

.V.attr.all.equal <-
function(target, current,
         check.type, check.class, check.attributes, ...) {
    msg <- msg. <- NULL
    if(check.type  && !identical(t1 <- typeof(target), t2 <- typeof(current)))
        msg <- c(msg, paste0("typeof(target) is ", deparse(t1), ", typeof(current) is ", deparse(t2)))
    if(check.class && !identical(c1 <-  class(target), c2 <-  class(current)))
        msg <- c(msg, paste0( "class(target) is ", deparse(c1),  ", class(current) is ", deparse(c2)))
    if(is.na(check.attributes) || check.attributes) {
        if((l1 <- length(target)) != (l2 <- length(current)))
            msg <- c(msg, paste0("length(target) is ", l1, ", length(current) is ", l2))
        a1 <- .V.attributes( target, is.na(check.attributes))
        a2 <- .V.attributes(current, is.na(check.attributes))
        if(!((is.null(a1) && is.null(a2)) ||
             isTRUE(ae <- all.equal.list(a1, a2, ...))))
            msg <- msg. <- c(msg, paste0("Attributes: < ", ae, " >"))
    }
    list(msg, is.null(msg) != is.null(msg.))
}

.V.all.equal <-
function(target, current,
         check.type = check.class,
         check.class = TRUE,
         check.attributes = TRUE, ...) {
    msg <- .V.attr.all.equal(target, current,
                             check.type = check.type,
                             check.class = check.class,
                             check.attributes = check.attributes, ...)
    if(!msg[[2L]]) {
        if(.isVector( target) || .isSparse( target) ||
           .isVector(current) || .isSparse(current)) {
            v1 <- as( target, "sparseVector")
            v2 <- as(current, "sparseVector")
            ae <- .V.a.e(v1, v2, ...)
        } else {
            v1 <- as( target, "vector")
            v2 <- as(current, "vector")
            ae <- all.equal(v1, v2, ...)
        }
        if(!isTRUE(ae))
            return(c(msg[[1L]], ae))
    }
    if(is.null(msg[[1L]])) TRUE else msg[[1L]]
}

setMethod("all.equal", c(target = "sparseVector", current = "vector"),
          .V.all.equal)

setMethod("all.equal", c(target = "vector", current = "sparseVector"),
          .V.all.equal)

setMethod("all.equal", c(target = "sparseVector", current = "sparseVector"),
          .V.all.equal)

## And for completeness:

setMethod("all.equal", c(target = "sparseVector", current = "Matrix"),
          .V.all.equal)
