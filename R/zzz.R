### Note that "in theory" even base::as.vector() should be overloaded.
### In practice that could be too much of a performance penalty in some cases.

.MatrixEnv <- new.env(parent = emptyenv(), hash = FALSE)#  e.g., for  once-per-session warnings

.chm_common <- new.env(parent = emptyenv())
## environment in which to store some settings from cholmod_common

.onLoad <- function(libname, pkgname)
{
    .Call(CHM_set_common_env, .chm_common)
    ## S4 method dispatch ambiguity warnings
    if(is.null(getOption("ambiguousMethodSelection"))) {
	if((!is.null(v <- getOption("Matrix.verbose")) && v >= 1) ||
	   isTRUE(getOption("Matrix.ambiguityNotes")) ||
	   interactive() && identical(Sys.info()[["user"]], "maechler")) {
	    ## nothing
	} else { # ambiguity notices are trashed
	    options(ambiguousMethodSelection = function(cond) {})
	    assign("no.methods.ambiguityNotes", TRUE, envir=.MatrixEnv)
	}
    }
}

## Instead, simply re-assign the [cr]bind()s which are recursively
## based on [cr]bind2 :
##
## save to cBind / rBind  ("rename")
if(getRversion() >= "3.2.0") {
    ## New (2015-02)  base :: cbind(), rbind() which dispatch on S4 "when needed":
    cBind <- function (..., deparse.level = 1) {
	## Once per session warning (or if "Matrix.(warn|verbose)"):
	if(is.null(wrn <- get0("warned.cBind", .MatrixEnv)) ||
	   isTRUE(getOption("Matrix.warn")) ||
	   isTRUE(getOption("Matrix.verbose"))) {
	    if(is.null(wrn))
		assign("warned.cBind", TRUE, envir=.MatrixEnv)
	    .Deprecated(msg = "'cBind' is deprecated.
 Since R version 3.2.0, base's cbind() should work fine with S4 objects")
	}
	base::cbind(..., deparse.level=deparse.level)
    }
    rBind <- function (..., deparse.level = 1) {
	## Once per session warning (or if "Matrix.(warn|verbose)"):
	if(is.null(wrn <- get0("warned.rBind", .MatrixEnv)) ||
           isTRUE(getOption("Matrix.warn")) ||
	   isTRUE(getOption("Matrix.verbose"))) {
	    if(is.null(wrn))
		assign("warned.rBind", TRUE, envir=.MatrixEnv)
	    .Deprecated(msg = "'rBind' is deprecated.
 Since R version 3.2.0, base's rbind() should work fine with S4 objects")
	}
	base::rbind(..., deparse.level=deparse.level)
    }

} else {
    cBind <- methods:::cbind
    rBind <- methods:::rbind
    lengths <- function (x, use.names = TRUE) vapply(x, length, 1L, USE.NAMES = use.names)
}

.onUnload <- function(libpath)
{
    library.dynam.unload("Matrix", libpath)
    if(isTRUE(.MatrixEnv $ no.methods.ambiguityNotes))# revert
	options(ambiguousMethodSelection = NULL)
}

.SuiteSparse_version <- function() {
    ssv <- .Call(get_SuiteSparse_version)
    package_version(list(major = ssv[1], minor = paste(ssv[2:3], collapse=".")))
}

if(getRversion() < "3.1.0") {
    if(getRversion() < "3.0.0") {
        rep_len <- function(x, length.out) rep(x, length.out=length.out)
    }
    anyNA <- function(x) any(is.na(x))
}
