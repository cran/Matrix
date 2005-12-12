### Note that "in theory" even base::as.vector() should be overloaded.
### In practice that could be too much of a performance penalty in some cases.

.onLoad <- function(libname, pkgname)
{
    require(methods)

    ## The following works around namespace-protection on purpose:
    assignInNamespace("..Old..as.matrix", base::as.matrix, ns = "base")
    assignInNamespace("..Old..as.array", base::as.array, ns = "base")
    assignInNamespace("as.matrix", as.matrix, ns = "base")
    assignInNamespace("as.array", as.array, ns = "base")
    ## Now all the functions in 'base' that start with something like
    ##  "x <- as.matrix(x)" or  "X <- as.array(X)"
    ## will work for 'Matrix'-matrices

    ## kronecker() / %x% -- in principle should re-assigne base::kronecker
    ## -----------> ?? performance hit ?? in mantelhaen.test() ??
    ##
    ## This is formally identical to the base definition, but should use the
    ## generic kronecker
    assignInNamespace("%x%", function (X, Y) kronecker(X, Y), ns = "base")

    if(paste(R.version$major, R.version$minor, sep=".") >= "2.2")
        methods:::bind_activation(TRUE)
}

.onUnload <- function(libpath)
{
    assignInNamespace("as.matrix", base::..Old..as.matrix, ns = "base")
    assignInNamespace("as.array",  base::..Old..as.array,  ns = "base")
    library.dynam.unload("Matrix", libpath)

    if(paste(R.version$major, R.version$minor, sep=".") >= "2.2")
        methods:::bind_activation(FALSE)
}
