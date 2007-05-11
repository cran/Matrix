### Note that "in theory" even base::as.vector() should be overloaded.
### In practice that could be too much of a performance penalty in some cases.

.onLoad <- function(libname, pkgname)
{
    require(methods)
    require(utils) # -> assignInNamespace {but "anyway"}

    Rv <- paste(R.version$major, R.version$minor, sep=".")

    ## The following works around namespace-protection on purpose:
    assignInNamespace("..Old..as.matrix", base::as.matrix, ns = "base")
    assignInNamespace("..Old..as.array", base::as.array, ns = "base")

    if(Rv >= "2.4") {
	## For R 2.4.0 and newer, need to also set the baseenv() --
	##  the following being really a hack:
	tmp <- function(x) {
	    if(methods:::seemsS4Object(x)) Matrix::as.matrix(x)
	    else UseMethod("as.matrix")
	}
	if(Rv >= "2.5")# change arglist to ' (x, ...) '
	    formals(tmp) <- alist(x=, ...=)
	environment(tmp) <- baseenv()
	assignInNamespace("as.matrix", tmp, ns = "base")
    } else {
	assignInNamespace("as.matrix", as.matrix, ns = "base")
    }
    ## does not (yet) need special treatment, since it's not S3 generic:
    assignInNamespace("as.array",  as.array, ns = "base")

    ## Now all the functions in 'base' that start with something like
    ##	"x <- as.matrix(x)" or	"X <- as.array(X)"
    ## will work for 'Matrix'-matrices

    ## kronecker() / %x% -- in principle should re-assign base::kronecker
    ## -----------> ?? performance hit ?? in mantelhaen.test() ??
    ##
    ## This is formally identical to the base definition, but should use the
    ## generic kronecker
    assignInNamespace("%x%", function (X, Y) kronecker(X, Y), ns = "base")

    if(Rv < "2.5") {
	## For	R  versions prior to 2.5.0 -- replace "diag<-" in base :
	tmp <- function(x, value) {
	    dx <- dim(x)
	    if(length(dx) != 2)
		## no further check, to also work with 'Matrix'
		stop("only matrix diagonals can be replaced")
	    len.i <- min(dx)
	    i <- seq_len(len.i)
	    len.v <- length(value)
	    if(len.v != 1 && len.v != len.i)
		stop("replacement diagonal has wrong length")
	    if(len.i > 0) x[cbind(i, i)] <- value
	    x
	}
	environment(tmp) <- baseenv()
	assignInNamespace("diag<-", tmp, ns = "base")
    }

##No more ## Activate the [cr]bind()s which are recursively based on [cr]bind2
##No more     methods:::bind_activation(TRUE)
}

## Instead, simply re-assign the [cr]bind()s which are recursively
## based on [cr]bind2 :
##
## save to cBind / rBind  ("rename")
cBind <- methods:::cbind
rBind <- methods:::rbind
## TODO? -- and export these {but users may need to use  base::cbind() ..?!}
## cbind <- cBind
## rbind <- rBind


.onUnload <- function(libpath)
{
    assignInNamespace("as.matrix", base::..Old..as.matrix, ns = "base")
    assignInNamespace("as.array",  base::..Old..as.array,  ns = "base")
    library.dynam.unload("Matrix", libpath)

##No more  ## deactivate the S4-aware [cr]bind()
##No more     methods:::bind_activation(FALSE)
}
