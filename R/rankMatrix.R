#### Determine *the* rank of a matrix
#### --------------------------------
##
## As this is not such a well-defined problem as people think,
## we provide *some* possibilities here, including the Matlab one.
##
## Ideas by Martin Maechler (April 2007) and Ravi Varadhan (October 2007)

qr2rankMatrix <- function(qr, tol = NULL, isBqr = is.qr(qr), do.warn=TRUE) {
    ## NB: 1) base::qr(*, LAPACK = TRUE/FALSE)  differ via attr(.,"useLAPACK")
    ##     2) if LAPACK=TRUE, .$rank is useless (always = full rank)
    ##
    ## return ( . ) :
    if(isBqr && !isTRUE(attr(qr, "useLAPACK")))
        qr$rank
    else {
        diagR <- if(isBqr) # hence "useLAPACK" here
                     diag(qr$qr) # faster than, but equivalent to   diag(qr.R(q.r))
                 else ## ==> assume Matrix::qr() i.e., currently "sparseQR"
                     ## FIXME: Here, we could be quite a bit faster,
                     ## by not returning the full sparseQR, but just
                     ## doing the following in C, and return the rank.
                     diag(qr@R)

        if(anyNA(diagR) || !all(is.finite(diagR))) {
            if(do.warn) {
                ifi <- is.finite(diagR)
                warning(gettextf(
                    "qr2rankMatrix(.): QR with only %d out of %d finite diag(R) entries",
                    sum(ifi), length(ifi)))
            }
            ## return
            NA_integer_
            ## alternative: gives *too* small rank in typical cases
            ## reduce the maximal rank by omitting all non-finite entries:
            ## diagR <- diagR[is.finite(diagR)]
            ## if(length(diagR) == 0)
            ##     return(NA_integer_)
        } else {
            if(isBqr) diagR <- abs(diagR) # in base qr(), sign( diag(R) ) are *not* coerced to positive
            else if(do.warn && any(diagR < 0))
                warning(gettextf("qr2rankMatrix(.): QR has negative diag(R) entries"))
            ## declare those entries to be zero that are < tol*max(.)
            if((mdi <- max(diagR, na.rm=TRUE)) > 0) {
                if(!is.numeric(tol)) {
                    ## d := dim(x) extracted from qr, in both (dense and sparse) qr() cases
                    d <- dim(if(isBqr) qr$qr else qr)
                    tol <- max(d) * .Machine$double.eps
                }
                sum(diagR >= tol * mdi)
                ## was sum(diag(q.r@R) != 0)
            }
            else 0L # for 0-matrix or all NaN or negative diagR[]
        }
    } ## else {Lapack or sparseQR}
}

rankMatrix <- function(x, tol = NULL,
                       method = c("tolNorm2", "qr.R", "qrLINPACK", "qr",
                                  "useGrad", "maybeGrad"),
                       sval = svd(x, 0,0)$d, warn.t = TRUE, warn.qr = TRUE)
{
    ## Purpose: rank of a matrix ``as Matlab'' or "according to Ravi V"
    ## ----------------------------------------------------------------------
    ## Arguments: x: a numerical matrix, maybe non-square
    ##          tol: numerical tolerance (compared to singular values)
    ##         sval: vector of non-increasing singular values of  x
    ##               (pass as argument if already known)
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 7 Apr 2007, 16:16
    ## ----------------------------------------------------------------------
    ##
    ## maybeGrad (Ravi V.): This algorithm determines the rank based on the
    ##	"gradient" of the
    ## absolute, singular values, rather than enforcing a rigid
    ## tolerance criterion,
    ##
    ## Author: Ravi Varadhan, Date: 22 October 2007 // Tweaks: MM, Oct.23
    ## ----------------------------------------------------------------------

    stopifnot(length(d <- dim(x)) == 2)
    p <- min(d)
    ## miss.meth <- missing(method)
    method <- match.arg(method)

    if(useGrad <- (method %in% c("useGrad", "maybeGrad"))) {
	stopifnot(length(sval) == p)
	if(p > 1) stopifnot(diff(sval) <= 0) # must be sorted non-increasingly: max = s..[1]
	if(sval[1] == 0) { ## <==> all singular values are zero  <==> Matrix = 0  <==> rank = 0
	    useGrad <- FALSE
	    method <- eval(formals()[["method"]])[[1]]
	} else {
	    ln.av <- log(abs(sval))
	    if(any(s0 <- sval == 0)) ln.av[s0] <- - .Machine$double.xmax # so we get diff() == 0
	    diff1 <- diff(ln.av)
	    if(method == "maybeGrad") {
		grad <- (min(ln.av) - max(ln.av)) / p
		useGrad <- !is.na(grad) && p > 1 && min(diff1) <= min(-3, 10 * grad)
	    }#  -------
	}
    }
    if(!useGrad) {
	x.dense <- is.numeric(x) || is(x,"denseMatrix")
        ## "qr" is allowed for backcompatibility [change @ 2013-11-24]
        if((Meth <- method) == "qr")
            method <- if(x.dense) "qrLINPACK" else "qr.R"
        else Meth <- substr(method, 1,2)

	if(Meth == "qr") {
	    if(is.null(tol)) tol <- max(d) * .Machine$double.eps
	} else { ## (Meth != "qr"), i.e. "tolNorm2"
	    if(is.null(tol)) {
		if(!x.dense && missing(sval) && prod(d) >= 100000L)
		    warning(gettextf(
 "rankMatrix(<large sparse Matrix>, method = '%s') coerces to dense matrix.
 Probably should rather use method = 'qr' !?",
				     method),
			    immediate.=TRUE, domain=NA)
                ## the "Matlab" default:
                if(p > 1) stopifnot(diff(sval) <= 0) #=> sval[1]= max(sval)
                tol <- max(d) * .Machine$double.eps
	    } else stopifnot((tol <- as.numeric(tol)[[1]]) >= 0)
	}
    }

    structure(## rank :
	      if(useGrad) which.min(diff1)
	      else if(Meth == "qr") {
		  if((do.t <- (d[1L] < d[2L])) && warn.t)
		      warning(gettextf(
			"rankMatrix(x, method='qr'): computing t(x) as nrow(x) < ncol(x)"))
		  q.r <- qr(if(do.t) t(x) else x, tol=tol, LAPACK = method != "qrLINPACK")
                  qr2rankMatrix(q.r, tol=tol, isBqr = x.dense, do.warn = warn.qr)
	      }
	      else if(sval[1] > 0) sum(sval >= tol * sval[1]) else 0L, ## "tolNorm2"
	      "method" = method,
	      "useGrad" = useGrad,
	      "tol" = if(useGrad) NA else tol)
}

## Ravi's plot of the absolute singular values:
if(FALSE) {
## if (plot.eigen) {
    plot(abs(sval), type = "b", xlab = "Index", xaxt = "n",
         log = "y", ylab = "|singular value|   [log scaled]")
    axis(1, at = unique(c(axTicks(1), rank, p)))
    abline(v = rank, lty = 3)
    mtext(sprintf("rank = %d (used %s (%g))", rank,
                  if(use.grad)"'gradient'" else "fixed tol.",
                  if(use.grad) min(diff1)  else tol))
}

