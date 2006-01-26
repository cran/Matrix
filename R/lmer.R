# Methods for lmer and for the objects that it produces

## Some utilities

## Return the pairs of expressions separated by vertical bars
findbars <- function(term)
{
    if (is.name(term) || is.numeric(term)) return(NULL)
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
    if (is.name(term) || is.numeric(term)) return(term)
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

## Expand an expression with colons to the sum of the lhs
## and the current expression
colExpand <- function(term)
{
    if (is.name(term) || is.numeric(term)) return(term)
    if (length(term) == 2) {
        term[[2]] <- colExpand(term[[2]])
        return(term)
    }
    stopifnot(length(term) == 3)
    if (is.call(term) && term[[1]] == as.name(':')) {
        return(substitute(A+B, list(A = term, B = colExpand(term[[2]]))))
    }
    term[[2]] <- colExpand(term[[2]])
    term[[3]] <- colExpand(term[[3]])
    term
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
         gradient = as.logical(gradient),
         Hessian = as.logical(Hessian))
}

rWishart <- function(n, df, invScal)
    .Call("Matrix_rWishart", n, df, invScal, PACKAGE = "Matrix")

setMethod("coef", signature(object = "mer"),
          function(object, ...)
      {
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
           val
       })

## setMethod("plot", signature(x = "lmer.coef"),
##           function(x, y, ...)
##       {
##           varying <- unique(do.call("c",
##                                     lapply(x, function(el)
##                                            names(el)[sapply(el,
##                                                             function(col)
##                                                             any(col != col[1]))])))
##           gf <- do.call("rbind", lapply(x, "[", j = varying))
##           gf$.grp <- factor(rep(names(x), sapply(x, nrow)))
##           switch(min(length(varying), 3),
##                  qqmath(eval(substitute(~ x | .grp,
##                                         list(x = as.name(varying[1])))), gf, ...),
##                  xyplot(eval(substitute(y ~ x | .grp,
##                                         list(y = as.name(varying[1]),
##                                              x = as.name(varying[2])))), gf, ...),
##                  splom(~ gf | .grp, ...))
##       })

setMethod("plot", signature(x = "lmer.ranef"),
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
                   model = TRUE, x = TRUE, y = TRUE , ...)
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
          fixed.form <- nobars(formula)  # remove any terms with `|'
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
          X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(,NROW(Y),0)

          weights <- model.weights(mf)
          offset <- model.offset(mf)
          ## check weights and offset
          if( !is.null(weights) && any(weights < 0) )
              stop("negative weights not allowed")
          if(!is.null(offset) && length(offset) != NROW(Y))
              stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(Y)), domain = NA)
          if (is.null(weights)) weights <- rep.int(1, length(Y))
          if (is.null(offset)) offset <- numeric(length(Y))

          ## fit a glm model to the fixed formula
##           fe$formula <- fixed.form
##           fe$subset <- NULL             # subset has already been created in call to data.frame
##           fe$data <- frm
##           fe$x <- fe$model <- fe$y <- TRUE
##           fe[[1]] <- as.name("glm")
##           glmFit <- eval(fe, parent.frame())
##           x <- glmFit$x
##           y <- as.double(glmFit$y)

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
          bars <- findbars(formula[[3]])
          names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
          fl <- lapply(bars,
                       function(x)
                       eval(substitute(as.factor(fac)[,drop = TRUE],
                                       list(fac = x[[3]])), mf))
          ## order factor list by decreasing number of levels
          nlev <- sapply(fl, function(x) length(levels(x)))
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
              mer <- .Call("mer_create", fl,
                           .Call("Zt_create", fl, Ztl, PACKAGE = "Matrix"),
                           X, Y, method, sapply(Ztl, nrow),
                           c(lapply(Ztl, rownames), list(.fixed = colnames(X))),
                           !(family$family %in% c("binomial", "poisson")),
                           match.call(), family,
                           PACKAGE = "Matrix")
              .Call("mer_ECMEsteps", mer, cv$niterEM, cv$EMverbose, PACKAGE = "Matrix")
              LMEoptimize(mer) <- cv
              return(new("lmer", mer,
                         frame = if (model) mf else data.frame(),
                         terms = mt))
          }

          ## The rest of the function applies to generalized linear mixed models
          gVerb <- getOption("verbose")
          glmFit <- glm.fit(X, Y, weights = weights, offset = offset, family = family,
                            intercept = attr(mt, "intercept") > 0)
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
          mer <- .Call("mer_create", fl,
                       .Call("Zt_create", fl, Ztl, PACKAGE = "Matrix"),
                       X, Y, method, sapply(Ztl, nrow),
                       c(lapply(Ztl, rownames), list(.fixed = colnames(X))),
                       !(family$family %in% c("binomial", "poisson")),
                       match.call(), family,
                       PACKAGE = "Matrix")

          GSpt <- .Call("glmer_init", environment(), PACKAGE = "Matrix")
          .Call("glmer_PQL", GSpt, PACKAGE = "Matrix")  # obtain PQL estimates
          PQLpars <- c(fixef(mer),
                       .Call("mer_coef", mer, 2, PACKAGE = "Matrix"))
          if (method == "PQL") {
              .Call("glmer_devLaplace", PQLpars, GSpt, PACKAGE = "Matrix")
              .Call("glmer_finalize", GSpt, PACKAGE = "Matrix")
              return(new("lmer", mer,
                         frame = if (model) mf else data.frame(),
                         terms = mt))
          }

          fixInd <- seq(ncol(X))
          ## pars[fixInd] == beta, pars[-fixInd] == theta
          ## indicator of constrained parameters
          const <- c(rep(FALSE, length(fixInd)),
                     unlist(lapply(mer@nc[seq(along = fl)],
                                   function(k) 1:((k*(k+1))/2) <= k)
                            ))
          devLaplace <- function(pars)
              .Call("glmer_devLaplace", pars, GSpt, PACKAGE = "Matrix")

          optimRes <-
              nlminb(PQLpars, devLaplace,
                     lower = ifelse(const, 5e-10, -Inf),
                     control = list(trace = cv$msVerbose,
                     iter.max = cv$msMaxIter))
          .Call("glmer_finalize", GSpt, PACKAGE = "Matrix")
          return(new("lmer", mer,
                     frame = if (model) mf else data.frame(),
                     terms = mt))

      })

## Extract the L matrix
setAs("mer", "dtCMatrix", function(from)
      .Call("mer_dtCMatrix", from, PACKAGE = "Matrix"))

## Extract the fixed effects
setMethod("fixef", signature(object = "mer"),
          function(object, ...)
          .Call("mer_fixef", object, PACKAGE = "Matrix"))

## Extract the random effects
setMethod("ranef", signature(object = "mer"),
          function(object, ...)
              new("lmer.ranef", .Call("mer_ranef", object, PACKAGE = "Matrix"))
          )

## Optimization for mer objects
setReplaceMethod("LMEoptimize", signature(x="mer", value="list"),
                 function(x, value)
             {
                 if (value$msMaxIter < 1) return(x)
                 nc <- x@nc
                 constr <- unlist(lapply(nc, function(k) 1:((k*(k+1))/2) <= k))
                 fn <- function(pars)
                     deviance(.Call("mer_coefGets", x, pars, 2, PACKAGE = "Matrix"))
                 gr <- if (value$gradient)
                     function(pars) {
                         if (!isTRUE(all.equal(pars,
                                               .Call("mer_coef", x,
                                                     2, PACKAGE = "Matrix"))))
                             .Call("mer_coefGets", x, pars, 2, PACKAGE = "Matrix")
                         .Call("mer_gradient", x, 2, PACKAGE = "Matrix")
                     }
                 else NULL
		 optimRes <- nlminb(.Call("mer_coef", x, 2, PACKAGE = "Matrix"),
                                    fn, gr,
                                    lower = ifelse(constr, 5e-10, -Inf),
                                    control = list(iter.max = value$msMaxIter,
                                    trace = as.integer(value$msVerbose)))
                 .Call("mer_coefGets", x, optimRes$par, 2, PACKAGE = "Matrix")
                 if (optimRes$convergence != 0) {
                     warning(paste("nlminb returned message",
                                   optimRes$message,"\n"))
                 }
                 return(x)
             })

setMethod("deviance", signature(object = "mer"),
          function(object, ...) {
              .Call("mer_factor", object, PACKAGE = "Matrix")
              object@deviance[[ifelse(object@method == "REML", "REML", "ML")]]
          })

setMethod("mcmcsamp", signature(object = "mer"),
          function(object, n = 1, verbose = FALSE, saveb = FALSE,
                   trans = TRUE, ...)
      {
          family <- object@family
          lmm <- family$family == "gaussian" && family$link == "identity"
          if (!lmm)
              stop("mcmcsamp for GLMMs not yet implemented in supernodal representation")
          ans <- t(.Call("mer_MCMCsamp", object, saveb, n,
                         trans, PACKAGE = "Matrix"))
          attr(ans, "mcpar") <- as.integer(c(1, n, 1))
          class(ans) <- "mcmc"
          glmer <- FALSE
          gnms <- names(object@flist)
          cnms <- object@cnames
          ff <- fixef(object)
          colnms <- c(names(ff), if (glmer) character(0) else "sigma^2",
                      unlist(lapply(seq(along = gnms),
                                    function(i)
                                    abbrvNms(gnms[i],cnms[[i]]))))
          if (trans) {
              ## parameter type: 0 => fixed effect, 1 => variance,
              ##                 2 => covariance
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
	  if(saveb) ## maybe better colnames, "RE.1","RE.2", ... ?
	      colnms <- c(colnms, rep.int("", length(object@rZy)))
          colnames(ans) <- colnms
          ans
      })

setMethod("simulate", signature(object = "mer"),
          function(object, nsim = 1, seed = NULL, ...)
      {
          if(!exists(".Random.seed", envir = .GlobalEnv))
              runif(1)               # initialize the RNG if necessary
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
          lpred <- .Call("mer_simulate", object, nsim, PACKAGE = "Matrix")
          sc <- 1
          if (object@useScale)
              sc <- .Call("mer_sigma", object, object@method == "REML",
                          PACKAGE = "Matrix")
          ## add fixed-effects contribution and per-observation noise term
          lpred <- as.data.frame(lpred + drop(object@X %*% fixef(object)) +
                                 rnorm(prod(dim(lpred)), sd = sc))
          ## save the seed
          attr(lpred, "seed") <- RNGstate
          lpred
      })


setMethod("show", "mer",
          function(object) {
              vcShow <- function(varc, useScale)
              {
                  digits <- max(3, getOption("digits") - 2)
                  sc <- attr(varc, "sc")
                  recorr <- lapply(varc, function(el) el@factors$correlation)
                  reStdDev <- c(lapply(recorr, slot, "sd"), list(Residual = sc))
                  reLens <- unlist(c(lapply(reStdDev, length)))
                  reMat <- array('', c(sum(reLens), 4),
                                 list(rep('', sum(reLens)),
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
                      colnames(corr) <- c("Corr", rep("", maxlen - 1))
                      reMat <- cbind(reMat, rbind(corr, rep("", ncol(corr))))
                  }
                  if (!useScale) reMat <- reMat[-nrow(reMat),]
                  print(reMat, quote = FALSE)
              }

              fcoef <- .Call("mer_fixef", object, PACKAGE = "Matrix")
              useScale <- object@useScale
              corF <- vcov(object)@factors$correlation
              #DF <- getFixDF(object)
              coefs <- cbind(fcoef, corF@sd)#, DF)
              dimnames(coefs) <-
                  list(names(fcoef), c("Estimate", "Std. Error"))#, "DF"))
              digits <- max(3, getOption("digits") - 2)
              REML <- object@method == "REML"
              llik <- logLik(object, REML)
              dev <- object@deviance
              devc <- object@devComp

              rdig <- 5
              if (glz <- !(object@method %in% c("REML", "ML"))) {
                  cat(paste("Generalized linear mixed model fit using",
                            object@method, "\n"))
              } else {
                  cat("Linear mixed-effects model fit by ")
                  cat(if(REML) "REML\n" else "maximum likelihood\n")
              }
              if (!is.null(object@call$formula)) {
                  cat("Formula:", deparse(object@call$formula),"\n")
              }
              if (!is.null(object@call$data)) {
                  cat("   Data:", deparse(object@call$data), "\n")
              }
              if (!is.null(object@call$subset)) {
                  cat(" Subset:",
                      deparse(asOneSidedFormula(object@call$subset)[[2]]),"\n")
              }
              if (glz) {
                  cat(" Family: ", object@family$family, "(",
                      object@family$link, " link)\n", sep = "")
                  print(data.frame(AIC = AIC(llik), BIC = BIC(llik),
                                   logLik = c(llik),
                                   deviance = -2*llik,
                                   row.names = ""))
              } else {
                  print(data.frame(AIC = AIC(llik), BIC = BIC(llik),
                                   logLik = c(llik),
                                   MLdeviance = dev["ML"],
                                   REMLdeviance = dev["REML"],
                                   row.names = ""))
              }
              cat("Random effects:\n")
              vcShow(VarCorr(object), useScale)
              ngrps <- lapply(object@flist, function(x) length(levels(x)))
              cat(sprintf("# of obs: %d, groups: ", devc[1]))
              cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
              cat("\n")
              if (!useScale)
                  cat("\nEstimated scale (compare to 1) ",
                      .Call("mer_sigma", object, FALSE, PACKAGE = "Matrix"),
                      "\n")
              if (nrow(coefs) > 0) {
                  if (useScale) {
                      stat <- coefs[,1]/coefs[,2]
                      #pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
                      nms <- colnames(coefs)
                      coefs <- cbind(coefs, stat) #, pval)
                      colnames(coefs) <- c(nms, "t value")#, "Pr(>|t|)")
                  } else {
                      coefs <- coefs[, 1:2, drop = FALSE]
                      stat <- coefs[,1]/coefs[,2]
                      pval <- 2*pnorm(abs(stat), lower = FALSE)
                      nms <- colnames(coefs)
                      coefs <- cbind(coefs, stat, pval)
                      colnames(coefs) <- c(nms, "z value", "Pr(>|z|)")
                  }
                  cat("\nFixed effects:\n")
                  printCoefmat(coefs, zap.ind = 3)#, tst.ind = 4)
                  rn <- rownames(coefs)
                  if (!is.null(corF)) {
                      p <- ncol(corF)
                      if (p > 1) {
                          cat("\nCorrelation of Fixed Effects:\n")
                          corF <- matrix(format(round(corF@x, 3), nsmall = 3),
                                         nc = p)
                          dimnames(corF) <- list(
                                                 abbreviate(rn, minlen=11),
                                                 abbreviate(rn, minlen=6))
                          corF[!lower.tri(corF)] <- ""
                          print(corF[-1, -p, drop=FALSE], quote = FALSE)
                      }
                  }
              }
              invisible(object)
          })

setMethod("vcov", signature(object = "mer"),
          function(object, REML = object@method == "REML",
                   useScale = object@useScale,...) {
              sc <- if (object@useScale) {
                  .Call("mer_sigma", object, REML, PACKAGE = "Matrix")
              } else { 1 }
              rr <- as(sc^2 * tcrossprod(solve(object@RXX)), "dpoMatrix")
              rr@factors$correlation <- as(rr, "correlation")
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
                  length(.Call("mer_coef", object, 0, PACKAGE = "Matrix"))
              attr(val, "REML") <- REML
              class(val) <- "logLik"
              val
          })

setMethod("VarCorr", signature(x = "mer"),
          function(x, REML = x@method == "REML", useScale = x@useScale, ...)
      {
          sc <- 1
          if (useScale)
              sc <- .Call("mer_sigma", x, REML, PACKAGE = "Matrix")
          sc2 <- sc * sc
          cnames <- x@cnames
          ans <- x@Omega
          for (i in seq(a=ans)) {
              el <- as(sc2 * solve(ans[[i]]), "dpoMatrix")
              el@Dimnames <- list(cnames[[i]], cnames[[i]])
              el@factors$correlation <- as(el, "correlation")
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
          modp <- logical(0)
          if (length(dots))
              modp <- sapply(dots, is, "mer") | sapply(dots, is, "lm")
          if (any(modp)) {              # multiple models - form table
              opts <- dots[!modp]
              mods <- c(list(object), dots[modp])
              names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)], as.character)
              mods <- mods[order(sapply(lapply(mods, logLik, REML = FALSE), attr, "df"))]
              calls <- lapply(mods, slot, "call")
              data <- lapply(calls, "[[", "data")
              if (any(data != data[[1]])) stop("all models must be fit to the same data object")
              header <- paste("Data:", data[[1]])
              subset <- lapply(calls, "[[", "subset")
              if (any(subset != subset[[1]])) stop("all models must use the same subset")
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
          } else {
              foo <- object
              #foo@status["factored"] <- FALSE
              #.Call("mer_factor", foo, PACKAGE="Matrix")
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
#                       c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
                       c("Df", "Sum Sq", "Mean Sq"))
              if ("(Intercept)" %in% nmeffects) table <- table[-match("(Intercept)", nmeffects), ]
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
          .Call("mer_fitted", object, PACKAGE = "Matrix")
          )

setMethod("formula", signature(x = "mer"),
          function(x, ...)
          x@call$formula
          )

setMethod("residuals", signature(object = "mer"),
          function(object, ...)
          .NotYetImplemented()
          )

setMethod("resid", signature(object = "mer"),
          function(object, ...)
          .NotYetImplemented()
          )

setMethod("summary", signature(object = "mer"),
          function(object, ...) object
          )

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
        .Call("mer_update_y", fm0, yy, PACKAGE = "Matrix")
        LMEoptimize(fm0) <- cv
        .Call("mer_update_y", fma, yy, PACKAGE = "Matrix")
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
        if (method == "AGQ") {    # determine nAGQ at PQL estimates
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
            .Call("glmer_devAGQ", pars, GSpt, nAGQ, PACKAGE = "Matrix")
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
        .Call("lmer_coefGets", mer, optpars[-fixInd], 2, PACKAGE = "Matrix")
    }

    .Call("glmer_finalize", GSpt, PACKAGE = "Matrix")
    loglik[] <- -deviance/2
}
