.onLoad <- function(libname, pkgname)
{
    ## The following works around namespace-protection on purpose:
    assignInNamespace("..Old..as.matrix", base::as.matrix, ns = "base")
    assignInNamespace("as.matrix", as.matrix, ns = "base")
    ## Now all the functions in 'base' that start with something like
    ##  "x <- as.matrix(x)"
    ## will work for 'Matrix'-matrices

    if(paste(R.version$major, R.version$minor, sep=".") >= "2.2")
        methods:::bind_activation(TRUE)
}

.onUnload <- function(libpath)
{
    assignInNamespace("as.matrix", base::..Old..as.matrix, ns = "base")
    library.dynam.unload("Matrix", libpath)

    if(paste(R.version$major, R.version$minor, sep=".") >= "2.2")
        methods:::bind_activation(FALSE)
}
