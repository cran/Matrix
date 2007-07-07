### Note that "in theory" even base::as.vector() should be overloaded.
### In practice that could be too much of a performance penalty in some cases.

.onLoad <- function(libname, pkgname)
{
    require(methods)
    require(utils) # -> assignInNamespace {but "anyway"}

    Rv <- paste(R.version$major, R.version$minor, sep=".")

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
