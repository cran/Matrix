### Note that "in theory" even base::as.vector() should be overloaded.
### In practice that could be too much of a performance penalty in some cases.

.onLoad <- function(libname, pkgname)
{
    require(methods)
    require(utils) # -> assignInNamespace {but "anyway"}

    ## GOAL: all the functions in 'base' that start with something like
    ##	"x <- as.matrix(x)" or	"X <- as.array(X)"
    ## will work for 'Matrix'-matrices :

    ## works around namespace-protection on purpose:
    assignInNamespace("..Old..as.matrix", base::as.matrix, ns = "base")
    assignInNamespace("..Old..as.array",  base::as.array, ns = "base")

    ##  hack because base::as.matrix() is an S3 generic :
    tmp <- function(x, ...) if(isS4(x)) Matrix::as.matrix(x) else UseMethod("as.matrix")
    environment(tmp) <- baseenv()
    assignInNamespace("as.matrix", tmp,      ns = "base")
    assignInNamespace("as.array",  as.array, ns = "base")


    ## kronecker() / %x% -- in principle should re-assign base::kronecker
    ## -----------> ?? performance hit ?? in mantelhaen.test() ??
    ##
    ## This is formally identical to the base definition, but should use the
    ## generic kronecker
    assignInNamespace("%x%", function (X, Y) kronecker(X, Y), ns = "base")

}

## Instead, simply re-assign the [cr]bind()s which are recursively
## based on [cr]bind2 :
##
## save to cBind / rBind  ("rename")
cBind <- methods:::cbind
rBind <- methods:::rbind


.onUnload <- function(libpath)
{
    assignInNamespace("as.matrix", base::..Old..as.matrix, ns = "base")
    assignInNamespace("as.array",  base::..Old..as.array,  ns = "base")
    library.dynam.unload("Matrix", libpath)
}

if(R.version$`svn rev` < 48201) { ## later: if(getRversion() < "2.9.0")
    ## also   export() it in  ../NAMESPACE
    .selectSuperClasses <- function(ext, dropVirtual = FALSE, namesOnly = TRUE,
                                    directOnly = TRUE, simpleOnly = directOnly)
    {
        ## No argument checking here
        addCond <- function(xpr, prev)
            if(length(prev)) substitute(P && N, list(P = prev, N = xpr)) else xpr
        C <- if(dropVirtual) {
            isVirtualExt <- function(x) getClass(x@superClass)@virtual
            quote(!isVirtualExt(exti))
        } else expression()
        if(directOnly) C <- addCond(quote(length(exti@by) == 0), C)
        if(simpleOnly) C <- addCond(quote(exti@simple), C)
        if(length(C)) {
            F <- function(exti){}; body(F) <- C
            ext <- ext[unlist(lapply(ext, F), use.names=FALSE)]
        }
        if(namesOnly) names(ext) else ext
    }
}

## A wrapper for now [as long as  'methods' has no *exported* version]:
.M.classEnv <- function (Class) methods:::.classEnv(Class)
