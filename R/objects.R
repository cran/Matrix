## if strict=FALSE then gives "...Matrix" or ".sparseVector" or ""
## if strict= TRUE then may also give one of these:
## "dpoMatrix", "dppMatrix", "corMatrix", "copMatrix", "pMatrix"
.M.nonvirtual <- function(x, strict = FALSE)
    .Call(R_Matrix_nonvirtual, x, strict)

## "[nlidz]" for Matrix, sparseVector, logical, integer, double, complex 'x';
## otherwise ""
.M.kind  <- function(x) .Call(R_Matrix_kind , x)

## "[gshtdi]" for Matrix, sparseVector 'x'; otherwise ""
.M.shape <- function(x) .Call(R_Matrix_shape, x)

## "[upCRTdi]" for Matrix 'x'; otherwise ""
.M.repr  <- function(x) .Call(R_Matrix_repr , x)

.isMatrix   <- function(x)
    nzchar(cl <- .M.nonvirtual(x)) && substr(cl, 4L, 4L) == "M"
.isVector   <- function(x)
    nzchar(cl <- .M.nonvirtual(x)) && substr(cl, 4L, 4L) != "M"
.isUnpacked <- function(x) .M.repr(x) == "u"
.isPacked   <- function(x) .M.repr(x) == "p"
.isC        <- function(x) .M.repr(x) == "C"
.isR        <- function(x) .M.repr(x) == "R"
.isT        <- function(x) .M.repr(x) == "T"
.isDiagonal <- function(x) .M.repr(x) == "d"
.isInd      <- function(x) .M.repr(x) == "i"
.isDense    <- function(x) any(.M.repr(x) == c("u", "p"))
.isSparse   <- function(x) any(.M.repr(x) == c("C", "R", "T", "d", "i"))
.isCRT      <- function(x) any(.M.repr(x) == c("C", "R", "T"))

## for .type.kind[.M.kind(x)]:
.type.kind <- c("n" = "logical",
                "l" = "logical",
                "i" = "integer",
                "d" = "double",
                "z" = "complex")

## for .kind.type[ typeof(x)]:
.kind.type <- c("logical" = "l",
                "integer" = "i",
                "double"  = "d",
                "complex" = "z")

extends1of <- function(class, classes, ...) {
    if(is.character(class))
        class <- getClassDef(class[[1L]])
    for(cl in classes)
        if(extends(class, cl, ...))
            return(TRUE)
    FALSE
}

MatrixClass <- function(cl, cld = getClassDef(cl),
                        ...Matrix = TRUE, dropVirtual = TRUE, ...) {
    if(!is.character(cl) || length(cl) != 1L || is.na(cl))
        stop("'cl' is not a character string")
    if(is.null(pkg <- cld@package) && is.null(pkg <- attr(cl, "package")))
        return(character(0L))
    if(identical(pkg, "Matrix") && (!...Matrix ||
           grepl("^[nlidz](ge|sy|sp|tr|tp|di|[gst][CRT])Matrix$", cl)))
        return(cl)
    r <- .selectSuperClasses(cld@contains, dropVirtual = dropVirtual,
                             namesOnly = TRUE, ...)
    if(length(r) == 0L)
        return(character(0L))
    while({
        r1 <- Recall(r[1L], ...Matrix = ...Matrix, dropVirtual = dropVirtual, ...)
        length(r1) == 0L && length(r) > 1L
    })
        r <- r[-1L]
    r1
}

class2 <- function(cl, kind = "l")
    sub("^[nlidz]", kind, MatrixClass(cl))

copyClass <- function(from, to.class,
                      sNames = intersect(slotNames(to.class), slotNames(from)),
                      check = TRUE) {
    to <- new(to.class)
    if(check)
        for(nm in sNames)
            slot(to, nm) <- slot(from, nm)
    else
        for(nm in sNames)
            attr(to, nm) <- attr(from, nm)
    to
}
