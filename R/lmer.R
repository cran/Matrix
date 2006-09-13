# Methods for lmer and for the objects that it produces

## To Do: Check if it would be worthwhile using a few ECME iterations
##   to stabilize the variance parameters at the beginning a Laplace
##   fit.

## To Do: Determine why the names of the components of the values of
##   the ranef and coef extractor methods are not printed.

## To Do: Change the output format for lmer objects to always print
##   the values of the log-likelihood and the
##   restricted-log-likelihood.  Base AIC and BIC on the log-likelihood.

## Some utilities

## Return the pairs of expressions separated by vertical bars
findbars <- function(term)
{
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    if (!is.call(term)) stop("term must be of class call")
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}

## Return the formula omitting the pairs of expressions
## that are separated by vertical bars
nobars <- function(term)
{
    if (!('|' %in% all.names(term))) return(term)
    if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
    if (length(term) == 2) {
	nb <- nobars(term[[2]])
	if (is.null(nb)) return(NULL)
	term[[2]] <- nb
	return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

## Substitute the '+' function for the '|' function
subbars <- function(term)
{
    if (is.name(term) || !is.language(term)) return(term)
    if (length(term) == 2) {
	term[[2]] <- subbars(term[[2]])
	return(term)
    }
    stopifnot(length(term) == 3)
    if (is.call(term) && term[[1]] == as.name('|'))
	term[[1]] <- as.name('+')
    term[[2]] <- subbars(term[[2]])
    term[[3]] <- subbars(term[[3]])
    term
}

## Return the list of '/'-separated terms in an expression that
## contains slashes
slashTerms <- function(x) {
    if (!("/" %in% all.names(x))) return(x)
    if (x[[1]] != as.name("/"))
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}

## from a list of length 2 return recursive interaction terms
makeInteraction <- function(x) {
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}


factorNames2char <- function(nms, collapse = ", ") {
    ## utility in messages / print etc:
    nms <- sQuote(nms)
    if(length(nms) == 1) paste("factor", nms)
    else paste("factors", paste(nms, collapse = collapse))
}

## expand any slashes in the grouping factors returned by findbars
expandSlash <- function(bb) {
    if (!is.list(bb)) return(expandSlash(list(bb)))
    ## I really do mean lapply(unlist(... - unlist returns a
    ## flattened list in this case
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
            return(lapply(unlist(makeInteraction(trms)),
                          function(trm) substitute(foo|bar,
                                                   list(foo = x[[2]],
                                                        bar = trm))))
        x
    }))
}

abbrvNms <- function(gnm, cnms)
{
    ans <- paste(abbreviate(gnm), abbreviate(cnms), sep = '.')
    if (length(cnms) > 1) {
	anms <- lapply(cnms, abbreviate, minlength = 3)
	nmmat <- outer(anms, anms, paste, sep = '.')
	ans <- c(ans, paste(abbreviate(gnm, minlength = 3),
			    nmmat[upper.tri(nmmat)], sep = '.'))
    }
    ans
}

## Control parameters for lmer
lmerControl <-
  function(maxIter = 200, # used in ../src/lmer.c only
	   tolerance = sqrt(.Machine$double.eps),# ditto
	   msMaxIter = 200,
	   ## msTol = sqrt(.Machine$double.eps),
	   ## FIXME:  should be able to pass tolerances to nlminb()
	   msVerbose = getOption("verbose"),
	   niterEM = 15,
	   EMverbose = getOption("verbose"),
	   PQLmaxIt = 30,# FIXME: unused; PQL currently uses 'maxIter' instead
	   usePQL = TRUE,
	   gradient = TRUE,
	   Hessian = FALSE # unused _FIXME_
	   )
{
    list(maxIter = as.integer(maxIter),
	 tolerance = as.double(tolerance),
	 msMaxIter = as.integer(msMaxIter),
	 ## msTol = as.double(msTol),
	 msVerbose = as.integer(msVerbose),# "integer" on purpose
	 niterEM = as.integer(niterEM),
	 EMverbose = as.logical(EMverbose),
	 PQLmaxIt = as.integer(PQLmaxIt),
	 usePQL = as.logical(usePQL),
	 gradient = as.logical(gradient),
	 Hessian = as.logical(Hessian))
}

rWishart <- function(n, df, invScal)
    .Call(Matrix_rWishart, n, df, invScal)

setMethod("coef", signature(object = "mer"),
	  function(object, ...)
      {
          if (length(list(...)))
              warning(paste('arguments named "',
                            paste(names(list(...)), collapse = ", "),
                                  '" ignored', sep = ''))
          fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
          ref <- ranef(object)
          val <- lapply(ref, function(x) fef[rep(1, nrow(x)),,drop = FALSE])
          for (i in seq(a = val)) {
              refi <- ref[[i]]
              row.names(val[[i]]) <- row.names(refi)
              nmsi <- colnames(refi)
              if (!all(nmsi %in% names(fef)))
                  stop("unable to align random and fixed effects")
              for (nm in nmsi) val[[i]][[nm]] <- val[[i]][[nm]] + refi[,nm]
          }
          new("coef.lmer", val)
       })

setMethod("plot", signature(x = "coef.lmer"),
          function(x, y, ...)
      {
          varying <- unique(do.call("c",
                                    lapply(x, function(el)
                                           names(el)[sapply(el,
                                                            function(col)
                                                            any(col != col[1]))])))
          gf <- do.call("rbind", lapply(x, "[", j = varying))
          gf$.grp <- factor(rep(names(x), sapply(x, nrow)))
          switch(min(length(varying), 3),
                 qqmath(eval(substitute(~ x | .grp,
                                        list(x = as.name(varying[1])))), gf, ...),
                 xyplot(eval(substitute(y ~ x | .grp,
                                        list(y = as.name(varying[1]),
                                             x = as.name(varying[2])))), gf, ...),
                 splom(~ gf | .grp, ...))
      })

setMethod("plot", signature(x = "ranef.lmer"),
	  function(x, y, ...)
      {
	  lapply(x, function(x) {
	      cn <- lapply(colnames(x), as.name)
	      switch(min(ncol(x), 3),
		     qqmath(eval(substitute(~ x, list(x = cn[[1]]))), x, ...),
		     xyplot(eval(substitute(y ~ x,
					    list(y = cn[[1]],
						 x = cn[[2]]))), x, ...),
		     splom(~ x, ...))
	  })
      })

setMethod("with", signature(data = "lmer"),
	  function(data, expr, ...) {
	      dat <- eval(data@call$data)
	      if (!is.null(na.act <- attr(data@frame, "na.action")))
		  dat <- dat[-na.act, ]
	      lst <- c(list(. = data), data@flist, data@frame, dat)
	      eval(substitute(expr), lst[unique(names(lst))])
	  })

setMethod("terms", signature(x = "lmer"),
	  function(x, ...) x@terms)


setMethod("lmer", signature(formula = "formula"),
	  function(formula, data, family = gaussian,
		   method = c("REML", "ML", "PQL", "Laplace", "AGQ"),
		   control = list(), start,
		   subset, weights, na.action, offset, contrasts = NULL,
		   model = TRUE, ...)
      {
	  ## match and check parameters
	  if (length(formula) < 3) stop("formula must be a two-sided formula")
	  cv <- do.call("lmerControl", control)

	  ## Must evaluate the model frame first and then fit the glm using
	  ## that frame.  Otherwise missing values in the grouping factors
	  ## cause inconsistent numbers of observations.
	  mf <- match.call()
	  m <- match(c("data", "subset", "weights",
		       "na.action", "offset"), names(mf), 0)
	  mf <- mf[c(1, m)]
	  frame.form <- subbars(formula) # substitute `+' for `|'
	  fixed.form <- nobars(formula)	 # remove any terms with `|'
	  if (inherits(fixed.form, "name")) # RHS is empty - use a constant
	      fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))
	  environment(fixed.form) <- environment(frame.form) <- environment(formula)

	  ## evaluate a model frame for fixed and random effects
	  mf$formula <- frame.form
	  mf$drop.unused.levels <- TRUE
	  mf[[1]] <- as.name("model.frame")
	  fe <- mf
	  mf <- eval(mf, parent.frame())

	  ## get the terms for the fixed-effects only
	  fe$formula <- fixed.form
	  fe <- eval(fe, parent.frame())
	  mt <- attr(fe, "terms")   # allow model.frame to update them
	  ## response vector
	  Y <- model.response(mf, "numeric")
	  ## avoid problems with 1D arrays, but keep names
	  if(length(dim(Y)) == 1) {
	      nm <- rownames(Y)
	      dim(Y) <- NULL
	      if(!is.null(nm)) names(Y) <- nm
	  }
	  ## null model support
	  X <- if (!is.empty.model(mt))
	      model.matrix(mt, mf, contrasts) else matrix(,NROW(Y),0)

	  weights <- model.weights(mf)
	  offset <- model.offset(mf)
	  ## check weights and offset
	  if( !is.null(weights) && any(weights < 0) )
	      stop("negative weights not allowed")
	  if(!is.null(offset) && length(offset) != NROW(Y))
	      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
			    length(offset), NROW(Y)), domain = NA)
	  if (is.null(weights)) weights <- rep.int(1, NROW(Y))
	  if (is.null(offset)) offset <- numeric(NROW(Y))

	  if(is.character(family))
	      family <- get(family, mode = "function", envir = parent.frame())
	  if(is.function(family)) family <- family()
	  if(is.null(family$family)) {
	      print(family)
	      stop("'family' not recognized")
	  }
	  ## check for a linear mixed model
	  lmm <- family$family == "gaussian" && family$link == "identity"
	  if (lmm) { # linear mixed model
	      method <- match.arg(method)
	      if (method %in% c("PQL", "Laplace", "AGQ")) {
		  warning(paste('Argument method = "', method,
				'" is not meaningful for a linear mixed model.\n',
				'Using method = "REML".\n', sep = ''))
		  method <- "REML"
	      }
	  } else { # generalized linear mixed model
	      if (missing(method)) method <- "PQL"
	      else {
		  method <- match.arg(method)
		  if (method == "ML") method <- "PQL"
		  if (method == "REML")
		      warning('Argument method = "REML" is not meaningful ',
			      'for a generalized linear mixed model.',
			      '\nUsing method = "PQL".\n')
	      }
	  }
	  if (method == "AGQ")
	      stop('method = "AGQ" not yet implemented for supernodal representation')
	  ## create factor list for the random effects
	  bars <- expandSlash(findbars(formula[[3]]))
	  names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
	  fl <- lapply(bars,
		       function(x)
		       eval(substitute(as.factor(fac)[,drop = TRUE],
				       list(fac = x[[3]])), mf))
	  ## order factor list by decreasing number of levels
	  nlev <- sapply(fl, function(x) length(levels(x)))
	  if(any(nlev == 0))
	      stop("resulting factor(s) with 0 levels in random effects part:\n ",
		   paste(sQuote(names(nlev[nlev == 0])), collapse=", "))
	  if (any(diff(nlev) > 0)) {
	      ord <- rev(order(nlev))
	      bars <- bars[ord]
	      fl <- fl[ord]
	  }
	  ## create list of transposed model matrices for random effects
	  Ztl <- lapply(bars, function(x)
			t(model.matrix(eval(substitute(~ expr,
						       list(expr = x[[2]]))),
				       mf)))
	  if (lmm) {
	      ## Create the mixed-effects representation (mer) object
	      mer <- .Call(mer_create, fl,
			   .Call(Zt_create, fl, Ztl),
			   X, Y, method, sapply(Ztl, nrow),
			   c(lapply(Ztl, rownames), list(.fixed = colnames(X))),
			   !(family$family %in% c("binomial", "poisson")),
			   match.call(), family)
	      .Call(mer_ECMEsteps, mer, cv$niterEM, cv$EMverbose)
	      LMEoptimize(mer) <- cv
	      return(new("lmer", mer,
			 frame = if (model) mf else data.frame(),
			 terms = mt))
	  }

	  ## The rest of the function applies to generalized linear mixed models
	  gVerb <- getOption("verbose")
	  glmFit <- glm.fit(X, Y, weights = weights, offset = offset, family = family,
			    intercept = attr(mt, "intercept") > 0)
	  weights <- glmFit$prior.weights
	  eta <- glmFit$linear.predictors
	  Y <- as.double(glmFit$y)
	  wtssqr <- weights * weights
	  linkinv <- quote(family$linkinv(eta))
	  mu.eta <- quote(family$mu.eta(eta))
	  mu <- family$linkinv(eta)
	  variance <- quote(family$variance(mu))
	  dev.resids <- quote(family$dev.resids(Y, mu, wtssqr))
	  LMEopt <- get("LMEoptimize<-")
	  doLMEopt <- quote(LMEopt(x = mer, value = cv))
	  mer <- .Call(mer_create, fl,
		       .Call(Zt_create, fl, Ztl),
		       X, Y, method, sapply(Ztl, nrow),
		       c(lapply(Ztl, rownames), list(.fixed = colnames(X))),
		       !(family$family %in% c("binomial", "poisson")),
		       match.call(), family)

	  GSpt <- .Call(glmer_init, environment())
	  if (cv$usePQL) {
	      .Call(glmer_PQL, GSpt)  # obtain PQL estimates
	      PQLpars <- c(fixef(mer),
			   .Call(mer_coef, mer, 2))
	  } else {
	      PQLpars <- c(coef(glmFit),
			   .Call(mer_coef, mer, 2))
	  }
	  if (method == "PQL") {
	      .Call(glmer_devLaplace, PQLpars, GSpt)
	      .Call(glmer_finalize, GSpt)
	      return(new("glmer", mer,
			 frame = if (model) mf else data.frame(),
			 terms = mt, weights = weights))
	  }

	  fixInd <- seq(ncol(X))
	  ## pars[fixInd] == beta, pars[-fixInd] == theta
	  ## indicator of constrained parameters
	  const <- c(rep(FALSE, length(fixInd)),
		     unlist(lapply(mer@nc[seq(along = fl)],
				   function(k) 1:((k*(k+1))/2) <= k)
			    ))
	  devLaplace <- function(pars)
	      .Call(glmer_devLaplace, pars, GSpt)

	  optimRes <-
	      nlminb(PQLpars, devLaplace,
		     lower = ifelse(const, 5e-10, -Inf),
		     control = list(trace = cv$msVerbose,
		     iter.max = cv$msMaxIter))
	  .Call(glmer_finalize, GSpt)
	  return(new("glmer", mer,
		     frame = if (model) mf else data.frame(),
		     terms = mt,
                     weights = weights))

      })

## Extract the L matrix
setAs("mer", "dtCMatrix", function(from)
      .Call(mer_dtCMatrix, from))

## Extract the fixed effects
setMethod("fixef", signature(object = "mer"),
	  function(object, ...)
	  .Call(mer_fixef, object))

## Extract the random effects
setMethod("ranef", signature(object = "mer"),
	  function(object, postVar = FALSE, ...) {
	      ans <- new("ranef.lmer",
                         lapply(.Call(mer_ranef, object),
                                data.frame, check.names = FALSE))
              names(ans) <- names(object@flist)
              if (postVar) {
                  pV <- .Call(mer_postVar, object)
                  for (i in seq(along = ans))
                      attr(ans[[i]], "postVar") <- pV[[i]]
              }
              ans
	  })
## Optimization for mer objects
setReplaceMethod("LMEoptimize", signature(x="mer", value="list"),
		 function(x, value)
	     {
		 if (value$msMaxIter < 1) return(x)
		 nc <- x@nc
		 constr <- unlist(lapply(nc, function(k) 1:((k*(k+1))/2) <= k))
		 fn <- function(pars)
		     deviance(.Call(mer_coefGets, x, pars, 2))
		 gr <- if (value$gradient)
		     function(pars) {
			 if (!isTRUE(all.equal(pars,
					       .Call(mer_coef, x,
						     2))))
			     .Call(mer_coefGets, x, pars, 2)
			 .Call(mer_gradient, x, 2)
		     }
		 else NULL
		 optimRes <- nlminb(.Call(mer_coef, x, 2),
				    fn, gr,
				    lower = ifelse(constr, 5e-10, -Inf),
				    control = list(iter.max = value$msMaxIter,
				    trace = as.integer(value$msVerbose)))
                 estPar <- optimRes$par
		 .Call(mer_coefGets, x, estPar, 2)

                 ## check for convergence on boundary
		 if (any(bd <- (estPar[constr] < 1e-9))) {
		     bpar <- rep.int(FALSE, length(estPar))
		     bpar[constr] <- bd
		     bgrp <- split(bpar,
				   rep(seq(along = nc),
				       unlist(lapply(nc,
						     function(k) (k*(k+1))/2))))
		     bdd <- unlist(lapply(bgrp, any))
		     lens <- unlist(lapply(bgrp, length))
		     if (all(lens[bdd] == 1)) { # variance components only
			 warning("Estimated variance for ",
				 factorNames2char(names(x@flist)[bdd]),
				 " is effectively zero\n")
		     } else {
			 warning("Estimated variance-covariance for ",
				 factorNames2char(names(x@flist)[bdd]),
				 " is singular\n")
		     }
		 }
		 if (optimRes$convergence != 0) {
		     warning("nlminb returned message ", optimRes$message,"\n")
		 }
		 return(x)
	     })

setMethod("qqmath", signature(x = "ranef.lmer"),
          function(x, data, ...) {
              prepanel.ci <- function(x, y, se, subscripts, ...) {
                  y <- as.numeric(y)
                  se <- as.numeric(se[subscripts])
                  hw <- 1.96 * se
                  list(ylim = range(y - hw, y + hw, finite = TRUE))
              }
              panel.ci <- function(x, y, se, subscripts, pch = 16, ...)  {
                  panel.grid(h = -1,v = -1)
                  panel.abline(h = 0)
                  x <- as.numeric(x)
                  y <- as.numeric(y)
                  se <- as.numeric(se[subscripts])
                  ly <- y - 1.96 * se
                  uy <- y + 1.96 * se
                  panel.segments(x, y - 1.96*se, x, y + 1.96 * se,
                                 col = 'black')
                  panel.xyplot(x, y, pch = pch, ...)
              }
              f <- function(x) {
                  if (!is.null(pv <- attr(x, "postVar"))) {
                      cols <- 1:(dim(pv)[1])
                      se <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
                      nr <- nrow(x)
                      nc <- ncol(x)
                      ord <- unlist(lapply(x, order)) +
                          rep((0:(nc - 1)) * nr, each = nr)
                      rr <- 1:nr
                      ind <- gl(ncol(x), nrow(x), labels = names(x))
                      xyplot(unlist(x)[ord] ~
                             rep(qnorm((rr - 0.5)/nr), ncol(x)) | ind[ord],
                             se = se[ord], prepanel = prepanel.ci, panel = panel.ci,
                             scales = list(y = list(relation = "free")),
                             xlab = "Standard normal quantiles",
                             ylab = NULL, aspect = 1, ...)
                  } else {
                      qqmath(~values|ind, stack(x),
                             scales = list(y = list(relation = "free")),
                             xlab = "Standard normal quantiles",
                             ylab = NULL, ...)
                  }
              }
              lapply(x, f)
          })

setMethod("deviance", signature(object = "mer"),
	  function(object, REML = NULL, ...) {
              if (is.null(REML)) REML <- object@method == "REML"
	      .Call(mer_factor, object)
	      object@deviance[[ifelse(REML, "REML", "ML")]]
	  })

## Mangle the names of the columns of the mcmcsamp result ans
## This operation is common to the methods for "lmer" and "glmer"
mcmccompnames <- function(ans, object, saveb, trans, glmer)
{
    gnms <- names(object@flist)
    cnms <- object@cnames
    ff <- fixef(object)
    colnms <- c(names(ff), if (glmer) character(0) else "sigma^2",
                unlist(lapply(seq(along = gnms),
                              function(i)
                              abbrvNms(gnms[i],cnms[[i]]))))
    if (trans) {
        ## parameter type: 0 => fixed effect, 1 => variance,
        ##		 2 => covariance
        ptyp <- c(integer(length(ff)), if (glmer) integer(0) else 1:1,
                  unlist(lapply(seq(along = gnms),
                                function(i)
                            {
                                k <- length(cnms[[i]])
                                rep(1:2, c(k, (k*(k-1))/2))
                            })))
        colnms[ptyp == 1] <-
            paste("log(", colnms[ptyp == 1], ")", sep = "")
        colnms[ptyp == 2] <-
            paste("atanh(", colnms[ptyp == 2], ")", sep = "")
    }
    colnms <- c(colnms, "deviance")
    if(saveb) {## maybe better colnames, "RE.1","RE.2", ... ?
        rZy <- object@rZy
        colnms <- c(colnms,
                    paste("b", sprintf(paste("%0",
                                             1+floor(log(length(rZy),10)),
                                             "d", sep = ''),
                                       seq(along = rZy)),
                          sep = '.'))
    }
    colnames(ans) <- colnms
    ans
}

setMethod("mcmcsamp", signature(object = "lmer"),
	  function(object, n = 1, verbose = FALSE, saveb = FALSE,
		   trans = TRUE, ...)
      {
          ans <- t(.Call(mer_MCMCsamp, object, saveb, n, trans, verbose))
	  attr(ans, "mcpar") <- as.integer(c(1, n, 1))
	  class(ans) <- "mcmc"
          mcmccompnames(ans, object, saveb, trans, FALSE)
      })

setMethod("mcmcsamp", signature(object = "glmer"),
	  function(object, n = 1, verbose = FALSE, saveb = FALSE,
		   trans = TRUE, ...)
      {
          family <- object@family
          mer <- as(object, "mer")
          weights <- object@weights
          cv <- lmerControl()
          eta <- .Call(mer_fitted, mer)
          offset <- numeric(length(eta)) ## change this, save the offset in mer
          Y <- object@y
          wtssqr <- weights * weights
          linkinv <- quote(family$linkinv(eta))
          mu.eta <- quote(family$mu.eta(eta))
          mu <- family$linkinv(eta)
          variance <- quote(family$variance(mu))
          dev.resids <- quote(family$dev.resids(Y, mu, wtssqr))
          LMEopt <- get("LMEoptimize<-")
          doLMEopt <- quote(LMEopt(x = mer, value = cv))
          GSpt <- .Call(glmer_init, environment())
          ans <- t(.Call(glmer_MCMCsamp, GSpt, saveb, n, trans, verbose))
          .Call(glmer_finalize, GSpt)
	  attr(ans, "mcpar") <- as.integer(c(1, n, 1))
	  class(ans) <- "mcmc"
          mcmccompnames(ans, object, saveb, trans, TRUE)
      })

setMethod("simulate", signature(object = "mer"),
	  function(object, nsim = 1, seed = NULL, ...)
      {
	  if(!exists(".Random.seed", envir = .GlobalEnv))
	      runif(1)		     # initialize the RNG if necessary
	  if(is.null(seed))
	      RNGstate <- .Random.seed
	  else {
	      R.seed <- .Random.seed
	      set.seed(seed)
	      RNGstate <- structure(seed, kind = as.list(RNGkind()))
	      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
	  }

	  family <- object@family
	  if (family$family != "gaussian" ||
	      family$link != "identity")
	      stop("simulation of generalized linear mixed models not yet implemented")
	  ## similate the linear predictors
	  lpred <- .Call(mer_simulate, object, nsim)
	  sc <-
              if (object@useScale)
                  .Call(mer_sigma, object, object@method == "REML")
              else 1

	  ## add fixed-effects contribution and per-observation noise term
	  lpred <- as.data.frame(lpred + drop(object@X %*% fixef(object)) +
				 rnorm(prod(dim(lpred)), sd = sc))
	  ## save the seed
	  attr(lpred, "seed") <- RNGstate
	  lpred
      })

simulestimate <- function(x, FUN, nsim = 1, seed = NULL, control = list())
{
    FUN <- match.fun(FUN)
    nsim <- as.integer(nsim[1])
    stopifnot(nsim > 0)
    if (!is.null(seed)) set.seed(seed)
    stopifnot(inherits(x, "lmer"))
    ## similate the linear predictors
    lpred <- .Call(mer_simulate, x, nsim)
    sc <- 1
    if (x@useScale)
        sc <- .Call(mer_sigma, x, x@method == "REML")
    ## add fixed-effects contribution and per-observation noise term
    lpred <- lpred + drop(x@X %*% fixef(x)) + rnorm(prod(dim(lpred)), sd = sc)

    cv <- do.call(lmerControl, control)
    Omega <- x@Omega
    x@wrkres <- x@y <- lpred[,1]
    .Call(mer_update_ZXy, x)
    LMEoptimize(x) <- cv
    template <- FUN(x)
    if (!is.numeric(template))
        stop("simulestimate currently only handles functions that return numeric vectors")
    ans <- matrix(template, nr = nsim, nc = length(template), byrow = TRUE)
    colnames(ans) <- names(template)
    for (i in 1:nsim) {
        x@wrkres <- x@y <- lpred[,i]
        x@Omega <- Omega
        .Call(mer_update_ZXy, x)
        LMEoptimize(x) <- cv
        foo <- try(FUN(x))
        ans[i,] <- if (inherits(foo, "try-error")) NA else foo
    }
    ans
}


formatVC <- function(varc, digits = max(3, getOption("digits") - 2))
{  ## "format()" the 'VarCorr'	matrix of the random effects -- for show()ing
    sc <- attr(varc, "sc")
    recorr <- lapply(varc, function(el) el@factors$correlation)
    reStdDev <- c(lapply(recorr, slot, "sd"), list(Residual = sc))
    reLens <- unlist(c(lapply(reStdDev, length)))
    nr <- sum(reLens)
    reMat <- array('', c(nr, 4),
		   list(rep.int('', nr),
			c("Groups", "Name", "Variance", "Std.Dev.")))
    reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
    reMat[,2] <- c(unlist(lapply(reStdDev, names)), "")
    reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
    reMat[,4] <- format(unlist(reStdDev), digits = digits)
    if (any(reLens > 1)) {
	maxlen <- max(reLens)
	corr <-
	    do.call("rbind",
		    lapply(recorr,
			   function(x, maxlen) {
			       x <- as(x, "matrix")
			       cc <- format(round(x, 3), nsmall = 3)
			       cc[!lower.tri(cc)] <- ""
			       nr <- dim(cc)[1]
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }, maxlen))
	colnames(corr) <- c("Corr", rep.int("", maxlen - 1))
	cbind(reMat, rbind(corr, rep.int("", ncol(corr))))
    } else reMat
}

## use S3 print method -- which can have non-trivial arguments
##     ^^ (such that print() remains S3 generic rather than S4)

## This is modeled a bit after  print.summary.lm :
print.mer <- function(x, digits = max(3, getOption("digits") - 3),
                      correlation = TRUE, symbolic.cor = x$symbolic.cor,
                      signif.stars = getOption("show.signif.stars"), ...)
{
    so <- summary(x)
    useScale <- so@useScale
    REML <- so@method == "REML"
    llik <- so@logLik
    dev <- so@deviance
    devc <- so@devComp
    glz <- so@isG

    cat(so@methTitle, "\n")
    if (!is.null(so@call$formula))
        cat("Formula:", deparse(so@call$formula),"\n")
    if (!is.null(so@call$data))
        cat("   Data:", deparse(so@call$data), "\n")
    if (!is.null(so@call$subset))
        cat(" Subset:",
            deparse(asOneSidedFormula(so@call$subset)[[2]]),"\n")
    if (glz)
        cat(" Family: ", so@family$family, "(",
            so@family$link, " link)\n", sep = "")
    print(so@AICtab, digits = digits)

    cat("Random effects:\n")
    print(so@REmat, quote = FALSE, digits = digits, ...)

    ngrps <- so@ngrps
    cat(sprintf("number of obs: %d, groups: ", devc[1]))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
    if (!useScale)
	cat("\nEstimated scale (compare to 1) ", so@sigma, "\n")
    if (nrow(so@coefs) > 0) {
	cat("\nFixed effects:\n")
	printCoefmat(so@coefs, zap.ind = 3, #, tst.ind = 4
		     digits = digits, signif.stars = signif.stars)
	if(correlation) {
	    rn <- rownames(so@coefs)
	    corF <- so@vcov@factors$correlation
	    if (!is.null(corF)) {
		p <- ncol(corF)
		if (p > 1) {
		    cat("\nCorrelation of Fixed Effects:\n")
		    if (is.logical(symbolic.cor) && symbolic.cor) {
			print(symnum(as(corF, "matrix"), abbr.col = NULL))
		    }
		    else {
			corF <- matrix(format(round(corF@x, 3), nsmall = 3),
				       nc = p)
			dimnames(corF) <- list(abbreviate(rn, minlen=11),
					       abbreviate(rn, minlen=6))
			corF[!lower.tri(corF)] <- ""
			print(corF[-1, -p, drop=FALSE], quote = FALSE)
		    }
		}
	    }
	}
    }
    invisible(x)
}

setMethod("show", "mer", function(object) print.mer(object))


setMethod("vcov", signature(object = "mer"),
	  function(object, REML = object@method == "REML",
		   useScale = object@useScale,...) {
	      sc <- if (object@useScale) {
		  .Call(mer_sigma, object, REML)
	      } else { 1 }
	      rr <- as(sc^2 * tcrossprod(solve(object@RXX)), "dpoMatrix")
	      rr@factors$correlation <- as(rr, "corMatrix")
	      rr
	  })


## calculates degrees of freedom for fixed effects Wald tests
## This is a placeholder.  The answers are generally wrong.  It will
## be very tricky to decide what a 'right' answer should be with
## crossed random effects.

setMethod("getFixDF", signature(object="mer"),
	  function(object, ...) {
	      devc <- object@devComp
	      rep(as.integer(devc[1]- devc[2]), devc[2])
	  })

setMethod("logLik", signature(object="mer"),
	  function(object, REML = object@method == "REML", ...) {
	      val <- -deviance(object, REML = REML)/2
	      devc <- as.integer(object@devComp[1:2])
	      attr(val, "nall") <- attr(val, "nobs") <- devc[1]
	      attr(val, "df") <- abs(devc[2]) +
		  length(.Call(mer_coef, object, 0))
	      attr(val, "REML") <- REML
	      class(val) <- "logLik"
	      val
	  })

setMethod("VarCorr", signature(x = "mer"),
	  function(x, REML = x@method == "REML", useScale = x@useScale, ...)
      {
	  sc <- if (useScale)
	      .Call(mer_sigma, x, REML) else 1
	  sc2 <- sc * sc
	  cnames <- x@cnames
	  ans <- x@Omega
	  for (i in seq(a = ans)) {
	      el <- as(sc2 * solve(ans[[i]]), "dpoMatrix")
	      el@Dimnames <- list(cnames[[i]], cnames[[i]])
	      el@factors$correlation <- as(el, "corMatrix")
	      ans[[i]] <- el
	  }
	  attr(ans, "sc") <- sc
	  ans
      })

setMethod("anova", signature(object = "mer"),
	  function(object, ...)
      {
	  mCall <- match.call(expand.dots = TRUE)
	  dots <- list(...)
	  modp <- if (length(dots))
	      sapply(dots, is, "mer") | sapply(dots, is, "lm") else logical(0)
	  if (any(modp)) {		# multiple models - form table
	      opts <- dots[!modp]
	      mods <- c(list(object), dots[modp])
	      names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)],
				    as.character)
	      mods <- mods[order(sapply(lapply(mods, logLik, REML = FALSE),
					attr, "df"))]
	      calls <- lapply(mods, slot, "call")
	      data <- lapply(calls, "[[", "data")
	      if (any(data != data[[1]]))
		  stop("all models must be fit to the same data object")
	      header <- paste("Data:", data[[1]])
	      subset <- lapply(calls, "[[", "subset")
	      if (any(subset != subset[[1]]))
		  stop("all models must use the same subset")
	      if (!is.null(subset[[1]]))
		  header <-
		      c(header, paste("Subset", deparse(subset[[1]]), sep = ": "))
	      llks <- lapply(mods, logLik, REML = FALSE)
	      Df <- sapply(llks, attr, "df")
	      llk <- unlist(llks)
	      chisq <- 2 * pmax(0, c(NA, diff(llk)))
	      dfChisq <- c(NA, diff(Df))
	      val <- data.frame(Df = Df,
				AIC = sapply(llks, AIC),
				BIC = sapply(llks, BIC),
				logLik = llk,
				"Chisq" = chisq,
				"Chi Df" = dfChisq,
				"Pr(>Chisq)" = pchisq(chisq, dfChisq, lower = FALSE),
				check.names = FALSE)
	      class(val) <- c("anova", class(val))
	      attr(val, "heading") <-
		  c(header, "Models:",
		    paste(names(mods),
			  unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
			  sep = ": "))
	      return(val)
	  }
	  else { ## ------ single model ---------------------
	      foo <- object
	      #foo@status["factored"] <- FALSE
	      #.Call(mer_factor, foo)
	      #dfr <- getFixDF(foo)
	      ss <- foo@rXy^2
	      ssr <- exp(foo@devComp["logryy2"])
	      names(ss) <- object@cnames[[".fixed"]]
	      asgn <- attr(foo@X, "assign")
	      terms <- foo@terms
	      nmeffects <- attr(terms, "term.labels")
	      if ("(Intercept)" %in% names(ss))
		  nmeffects <- c("(Intercept)", nmeffects)
	      ss <- unlist(lapply(split(ss, asgn), sum))
	      df <- unlist(lapply(split(asgn,  asgn), length))
	      #dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
	      ms <- ss/df
	      #f <- ms/(ssr/dfr)
	      #P <- pf(f, df, dfr, lower.tail = FALSE)
	      #table <- data.frame(df, ss, ms, dfr, f, P)
	      table <- data.frame(df, ss, ms)
	      dimnames(table) <-
		  list(nmeffects,
#			c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
		       c("Df", "Sum Sq", "Mean Sq"))
	      if ("(Intercept)" %in% nmeffects)
		  table <- table[-match("(Intercept)", nmeffects), ]
	      attr(table, "heading") <- "Analysis of Variance Table"
	      class(table) <- c("anova", "data.frame")
	      table
	  }
      })

setMethod("confint", signature(object = "mer"),
	  function(object, parm, level = 0.95, ...)
	  .NotYetImplemented()
	  )

setMethod("fitted", signature(object = "mer"),
	  function(object, ...)
	  .Call(mer_fitted, object)
	  )

setMethod("formula", signature(x = "mer"),
	  function(x, ...)
	  x@call$formula
	  )

setMethod("residuals", signature(object = "mer"),
	  function(object, ...) {
              fam <- object@family
              if (fam$family == "gaussian" && fam$link == "identity")
                  return(object@y - fitted(object))
              .NotYetImplemented()
	  })

## FIXME: There should not be two identical methods like this but I'm not
##        sure how to pass the ... argument to a method for another generic
##        cleanly.
setMethod("resid", signature(object = "mer"),
	  function(object, ...) {
              fam <- object@family
              if (fam$family == "gaussian" && fam$link == "identity")
                  return(object@y - fitted(object))
              .NotYetImplemented()
	  })

setMethod("summary", signature(object = "mer"),
	  function(object, ...) {

	      fcoef <- .Call(mer_fixef, object)
	      useScale <- object@useScale
	      vcov <- vcov(object)
	      corF <- vcov@factors$correlation
	      ## DF <- getFixDF(object)
	      coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd) #, DF = DF)
	      REML <- object@method == "REML"
	      llik <- logLik(object, REML)
	      dev <- object@deviance
	      devc <- object@devComp

	      glz <- !(object@method %in% c("REML", "ML"))
	      methTitle <-
		  if (glz)
		      paste("Generalized linear mixed model fit using",
			    object@method)
		  else paste("Linear mixed-effects model fit by",
			     if(REML) "REML" else "maximum likelihood")

	      AICframe <- {
		  if (glz)
		      data.frame(AIC = AIC(llik), BIC = BIC(llik),
				 logLik = c(llik),
				 deviance = -2*llik,
				 row.names = "")
		  else
		      data.frame(AIC = AIC(llik), BIC = BIC(llik),
				 logLik = c(llik),
				 MLdeviance = dev["ML"],
				 REMLdeviance = dev["REML"],
				 row.names = "")
	      }
	      REmat <- formatVC(VarCorr(object))
	      if (!useScale) REmat <- REmat[-nrow(REmat), , drop = FALSE]

	      if (nrow(coefs) > 0) {
		  if (useScale) {
		      stat <- coefs[,1]/coefs[,2]
		      ##pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
		      coefs <- cbind(coefs, "t value" = stat) #, "Pr(>|t|)" = pval)
		  } else {
		      coefs <- coefs[, 1:2, drop = FALSE]
		      stat <- coefs[,1]/coefs[,2]
		      pval <- 2*pnorm(abs(stat), lower = FALSE)
		      coefs <- cbind(coefs, "z value" = stat, "Pr(>|z|)" = pval)
		  }
	      } ## else : append columns to 0-row matrix ...

	      new(if(is(object, "lmer")) "summary.lmer" else "summary.mer",
		  object,
		  isG = glz,
		  methTitle = methTitle,
		  logLik = llik,
		  ngrps = sapply(object@flist, function(x) length(levels(x))),
		  sigma = .Call(mer_sigma, object, REML),
		  coefs = coefs,
		  vcov = vcov,
		  REmat = REmat,
		  AICtab= AICframe
		  )
	  })## summary()

## Methods for "summary.*" objects:
setMethod("vcov", signature(object = "summary.mer"),
	  function(object) object@vcov)
setMethod("logLik", signature(object = "summary.mer"),
	  function(object) object@logLik)
setMethod("deviance", signature(object = "summary.mer"),
 	  function(object) object@deviance)
setMethod("summary", signature(object = "summary.mer"), function(object) object)


setMethod("update", signature(object = "mer"),
	  function(object, formula., ..., evaluate = TRUE)
      {
	  call <- object@call
	  if (is.null(call))
	      stop("need an object with call slot")
	  extras <- match.call(expand.dots = FALSE)$...
	  if (!missing(formula.))
	      call$formula <- update.formula(formula(object), formula.)
	  if (length(extras) > 0) {
	      existing <- !is.na(match(names(extras), names(call)))
	      for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
	      if (any(!existing)) {
		  call <- c(as.list(call), extras[!existing])
		  call <- as.call(call)
	      }
	  }
	  if (evaluate)
	      eval(call, parent.frame())
	  else call
      })

simss <- function(fm0, fma, nsim)
{
    ysim <- simulate(fm0, nsim)
    cv <- list(gradient = FALSE, msMaxIter = 200:200,
	       msVerbose = 0:0)
    sapply(ysim, function(yy) {
	.Call(mer_update_y, fm0, yy)
	LMEoptimize(fm0) <- cv
	.Call(mer_update_y, fma, yy)
	LMEoptimize(fma) <- cv
	exp(c(H0 = fm0@devComp[["logryy2"]],
	      Ha = fma@devComp[["logryy2"]]))
    })
}

## Some leftover code from the old AGQ method in lmer.
if (FALSE) {
### FIXME: For nf == 1 change this to an AGQ evaluation.  Needs
### AGQ for nc > 1 first.
    fxd <- PQLpars[fixInd]
    loglik <- logLik(mer)

    if (method %in% c("Laplace", "AGQ")) {
	nAGQ <- 1
	if (method == "AGQ") {	  # determine nAGQ at PQL estimates
	    dev11 <- devAGQ(PQLpars, 11)
	    ## FIXME: Should this be an absolute or a relative tolerance?
	    devTol <- sqrt(.Machine$double.eps) * abs(dev11)
	    for (nAGQ in c(9, 7, 5, 3, 1))
		if (abs(dev11 - devAGQ(PQLpars, nAGQ - 2)) > devTol) break
	    nAGQ <- nAGQ + 2
	    if (gVerb)
		cat(paste("Using", nAGQ, "quadrature points per column\n"))
	}
	obj <- function(pars)
	    .Call(glmer_devAGQ, pars, GSpt, nAGQ)
	optimRes <-
	    nlminb(PQLpars, obj,
		   lower = ifelse(const, 5e-10, -Inf),
		   control = list(trace = getOption("verbose"),
		   iter.max = cv$msMaxIter))
	optpars <- optimRes$par
	if (optimRes$convergence != 0)
	    warning("nlminb failed to converge")
	deviance <- optimRes$objective
	if (gVerb)
	    cat(paste("convergence message", optimRes$message, "\n"))
	fxd[] <- optpars[fixInd]  ## preserve the names
	.Call(lmer_coefGets, mer, optpars[-fixInd], 2)
    }

    .Call(glmer_finalize, GSpt)
    loglik[] <- -deviance/2
}## end{leftover}

setMethod("isNested", "mer",
          function(x, ...) !(x@L@type[1]),
          valueClass = "logical")

setMethod("denomDF", "mer",
          function(x, ...)
      {
          mm <- x@X
          aa <- attr(mm, "assign")
          tt <- x@terms
          if (!isNested(x))
              return(list(coef = as.numeric(rep(NA, length(x@fixef))),
                          terms = as.numeric(rep(NA,
                          length(attr(tt, "order"))))))
          hasintercept <- attr(tt, "intercept") > 0
          ## check which variables vary within levels of grouping factors
          vars <- eval(attr(tt, "variables"), x@frame)
          fl <- x@flist
          vv <- matrix(0:0, nrow = length(vars), ncol = length(fl),
                        dimnames = list(NULL, names(fl)))
          ## replace this loop by C code.
          for (i in 1:nrow(ans))        # check if variables vary within factors
              for (j in 1:ncol(ans))
                  ans[i,j] <- all(tapply(vars[[i]], fl[[j]],
                                         function(x) length(unique(x)) == 1))
          ## which terms vary within levels of which grouping factors?
          tv <- crossprod(attr(tt, "factors"), !ans)
          ## maximum level at which the term is constant
          ml <- apply(tv, 1, function(rr) max(0, which(as.logical(rr))))
          ## unravel assignment applied to terms
          ll <- attr(tt, "term.labels")
          if (hasintercept)
              ll <- c("(Intercept)", ll)
          aaa <- factor(aa, labels = ll)
          asgn <- split(order(aa), aaa)
          nco <- lapply(asgn, length)   # number of coefficients per term
          nlev <- lapply(fl, function(x) length(levels(x)))
          if (hasintercept) asgn$"(Intercept)" <- NULL
          list(ml = ml, nco = nco, nlev = nlev)
      })

hatTrace <- function(x)
{
    stopifnot(is(x, "mer"))
    .Call(mer_hat_trace2, x)
}
