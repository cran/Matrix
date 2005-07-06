# Methods for lmer and for the objects that it produces

## Some utilities

## Return the index into the packed lower triangle
Lind <- function(i,j) {
    if (i < j) stop(paste("Index i=", i,"must be >= index j=", j))
    ((i - 1) * i)/2 + j
}

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
    # FIXME: is the is.name in the condition redundant?
    #   A name won't satisfy the first condition.
    if (!('|' %in% all.names(term)) || is.name(term)) return(term)
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
    if (is.call(term) && term[[1]] == as.name('|')) term[[1]] <- as.name('+')
    term[[2]] <- subbars(term[[2]])
    term[[3]] <- subbars(term[[3]])
    term
}
    
## Control parameters for lmer
lmerControl <-
  function(maxIter = 200,
           msMaxIter = 200,
           tolerance = sqrt((.Machine$double.eps)),
           niterEM = 15,
           msTol = sqrt(.Machine$double.eps),
           msVerbose = getOption("verbose"),
           PQLmaxIt = 30,
           EMverbose = getOption("verbose"),
           analyticGradient = TRUE,
           analyticHessian = FALSE)
{
    list(maxIter = as.integer(maxIter),
         msMaxIter = as.integer(msMaxIter),
         tolerance = as.double(tolerance),
         niterEM = as.integer(niterEM),
         msTol = as.double(msTol),
         msVerbose = as.logical(msVerbose),
         PQLmaxIt = as.integer(PQLmaxIt),
         EMverbose = as.logical(EMverbose),
         analyticGradient = as.logical(analyticGradient),
         analyticHessian = as.logical(analyticHessian))
}

setMethod("lmer", signature(formula = "formula"),
          function(formula, data, family,
                   method = c("REML", "ML", "PQL", "Laplace", "AGQ"),
                   control = list(),
                   subset, weights, na.action, offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {                                 ## match and check parameters
          if (length(formula) < 3) stop("formula must be a two-sided formula")
          cv <- do.call("lmerControl", control)
          if (lmm <- missing(family)) { # linear mixed model
              method <- match.arg(method)
              if (method %in% c("PQL", "Laplace", "AGQ")) {
                  warning(paste('Argument method = "', method,
                                '" is not meaningful for a linear mixed model.\n',
                                'Using method = "REML".\n', sep = ''))
                  method <- "REML"
              }
          } else {                      # generalized linear mixed model
              method <- if (missing(method)) "PQL" else match.arg(method)
              if (method == "ML") method <- "PQL"
              if (method == "REML") 
                  warning(paste('Argument method = "REML" is not meaningful',
                                'for a generalized linear mixed model.',
                                '\nUsing method = "PQL".\n'))
          }

          ## evaluate glm.fit, a generalized linear fit of fixed effects only
          mf <- match.call()           
          m <- match(c("family", "data", "subset", "weights",
                       "na.action", "offset"), names(mf), 0)
          mf <- mf[c(1, m)]
          frame.form <- subbars(formula) # substitute `+' for `|'
          fixed.form <- nobars(formula)  # remove any terms with `|'
          if (inherits(fixed.form, "name")) # RHS is empty - use a constant
              fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))
          environment(fixed.form) <- environment(frame.form) <- environment(formula)
          mf$formula <- fixed.form
          mf$x <- mf$model <- mf$y <- TRUE
          mf[[1]] <- as.name("glm")
          glm.fit <- eval(mf, parent.frame())
          x <- glm.fit$x
          y <- as.double(glm.fit$y)
          family <- glm.fit$family

          ## evaluate a model frame for fixed and random effects
          mf$formula <- frame.form
          mf$x <- mf$model <- mf$y <- mf$family <- NULL
          mf$drop.unused.levels <- TRUE
          mf[[1]] <- as.name("model.frame")
          frm <- eval(mf, parent.frame())

          ## grouping factors and model matrices for random effects
          bars <- findbars(formula[[3]])
          random <-
              lapply(bars,
                     function(x) list(model.matrix(eval(substitute(~term,
                                                                   list(term=x[[2]]))),
                                                   frm),
                                      eval(substitute(as.factor(fac)[,drop = TRUE],
                                                      list(fac = x[[3]])), frm)))
          names(random) <- unlist(lapply(bars, function(x) deparse(x[[3]])))

          ## order factor list by decreasing number of levels
          nlev <- sapply(random, function(x) length(levels(x[[2]])))
          if (any(diff(nlev) > 0)) {
              random <- random[rev(order(nlev))]
          }

          ## Create the model matrices and a mixed-effects representation (mer)
          mmats <- c(lapply(random, "[[", 1),
                     .fixed = list(cbind(glm.fit$x, .response = glm.fit$y)))
          mer <- .Call("lmer_create", lapply(random, "[[", 2),
                       mmats, method, PACKAGE = "Matrix")
          if (lmm) {                    ## linear mixed model
              .Call("lmer_initial", mer, PACKAGE="Matrix")
              .Call("lmer_ECMEsteps", mer, cv$niterEM, cv$EMverbose, PACKAGE = "Matrix")
              LMEoptimize(mer) <- cv
              fits <- .Call("lmer_fitted", mer, mmats, TRUE, PACKAGE = "Matrix")
              return(new("lmer",
                         mer,
                         assign = attr(x, "assign"),
                         call = match.call(),
                         family = family, fitted = fits,
                         fixed = fixef(mer),
                         frame = if (model) frm else data.frame(),
                         logLik = logLik(mer),
                         residuals = unname(model.response(frm) - fits),
                         terms = glm.fit$terms))
          }

          ## The rest of the function applies to generalized linear mixed models
          gVerb <- getOption("verbose")
          eta <- glm.fit$linear.predictors
          wts <- glm.fit$prior.weights
          wtssqr <- wts * wts
          offset <- glm.fit$offset
          if (is.null(offset)) offset <- numeric(length(eta))
          mu <- numeric(length(eta))

          dev.resids <- quote(family$dev.resids(y, mu, wtssqr))
          linkinv <- quote(family$linkinv(eta))
          mu.eta <- quote(family$mu.eta(eta))
          variance <- quote(family$variance(mu))
          LMEopt <- get("LMEoptimize<-")
          doLMEopt <- quote(LMEopt(x = mer, value = cv))

          GSpt <- .Call("glmer_init", environment())
          .Call("glmer_PQL", GSpt)  # obtain PQL estimates

          fixInd <- seq(ncol(x))
          ## pars[fixInd] == beta, pars[-fixInd] == theta
          PQLpars <- c(fixef(mer),
                       .Call("lmer_coef", mer, 2, PACKAGE = "Matrix"))
          ## set flag to skip fixed-effects in subsequent calls
          mer@nc[length(mmats)] <- -mer@nc[length(mmats)]
          ## indicator of constrained parameters
          const <- c(rep(FALSE, length(fixInd)),
                     unlist(lapply(mer@nc[seq(along = random)],
                                   function(k) 1:((k*(k+1))/2) <= k)
                            ))
          devAGQ <- function(pars, n)
              .Call("glmer_devAGQ", pars, GSpt, n, PACKAGE = "Matrix")
          
          deviance <- devAGQ(PQLpars, 1)
### FIXME: Change this to an AGQ evaluation once when nf == 1.  Needs
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
              if (exists("nlminb", mode = "function")) {
                  optimRes <-
                      nlminb(PQLpars, obj,
                             lower = ifelse(const, 5e-10, -Inf),
                             control = list(trace = getOption("verbose"),
                             iter.max = cv$msMaxIter))
                  optpars <- optimRes$par
                  if (optimRes$convergence != 0)
                      warning("nlminb failed to converge")
                  deviance <- optimRes$objective
              } else {
                  optimRes <-
                      optim(PQLpars, obj, method = "L-BFGS-B",
                            lower = ifelse(const, 5e-10, -Inf),
                            control = list(trace = getOption("verbose"),
                                           maxit = cv$msMaxIter))
                  optpars <- optimRes$par
                  if (optimRes$convergence != 0)
                      warning("optim failed to converge")
                  deviance <- optimRes$value
              }
              if (gVerb) {
                  cat(paste("convergence message", optimRes$message, "\n"))
              }
              fxd[] <- optpars[fixInd]  ## preserve the names
              .Call("lmer_coefGets", mer, optpars[-fixInd], 2,
                    PACKAGE = "Matrix")
          }

          .Call("glmer_finalize", GSpt, PACKAGE = "Matrix")
          loglik[] <- -deviance/2
          new("lmer", mer, frame = frm, terms = glm.fit$terms,
              assign = attr(glm.fit$x, "assign"),
              call = match.call(), family = family,
              logLik = loglik, fixed = fxd)
      })

setReplaceMethod("LMEoptimize", signature(x="mer", value="list"),
                 function(x, value)
             {
                 if (value$msMaxIter < 1) return(x)
                 nc <- x@nc
                 constr <- unlist(lapply(nc[1:(length(nc) - 2)],
                                         function(k) 1:((k*(k+1))/2) <= k))
                 fn <- function(pars)
                     deviance(.Call("lmer_coefGets", x, pars, 2, PACKAGE = "Matrix"))
                 gr <- NULL
                 if (value$analyticGradient)
                     gr <- 
                         function(pars) {
                             if (!isTRUE(all.equal(pars,
                                                   .Call("lmer_coef", x,
                                                         2, PACKAGE = "Matrix"))))
                                 .Call("lmer_coefGets", x, pars, 2, PACKAGE = "Matrix")
                             .Call("lmer_gradient", x, 2, PACKAGE = "Matrix")
                         }
                 optimRes <- 
                     if (exists("nlminb", mode = "function"))
                         nlminb(.Call("lmer_coef", x, 2, PACKAGE = "Matrix"),
                                fn, gr,
                                lower = ifelse(constr, 5e-10, -Inf),
                                control = list(iter.max = value$msMaxIter,
                                trace = as.integer(value$msVerbose)))
                     else
                         optim(.Call("lmer_coef", x, 2, PACKAGE = "Matrix"),
                               fn, gr, method = "L-BFGS-B",
                               lower = ifelse(constr, 5e-10, -Inf),
                               control = list(maxit = value$msMaxIter,
                               trace = as.integer(value$msVerbose)))
                 .Call("lmer_coefGets", x, optimRes$par, 2, PACKAGE = "Matrix")
                 if (optimRes$convergence != 0) {
                     warning(paste("optim or nlminb returned message",
                                   optimRes$message,"\n"))
                 }
                 return(x)
             })

setMethod("ranef", signature(object = "lmer"),
          function(object, accumulate = FALSE, ...) {
              val <- new("lmer.ranef",
                         lapply(.Call("lmer_ranef", object, PACKAGE = "Matrix"),
                                data.frame, check.names = FALSE),
                         varFac = object@bVar,
                         stdErr = .Call("lmer_sigma", object,
                         object@method == "REML", PACKAGE = "Matrix"))
              if (!accumulate || length(val@varFac) == 1) return(val)
              ## check for nested factors
              L <- object@L
              if (any(sapply(seq(a = val), function(i) length(L[[Lind(i,i)]]@i))))
                  error("Require nested grouping factors to accumulate random effects")
              val
          })

setMethod("fixef", signature(object = "mer"),
          function(object, ...)
              .Call("lmer_fixef", object, PACKAGE = "Matrix"))

setMethod("fixef", signature(object = "lmer"),
          function(object, ...) object@fixed)

setMethod("VarCorr", signature(x = "lmer"),
          function(x, REML = TRUE, useScale = TRUE, ...) {
              val <- .Call("lmer_variances", x, PACKAGE = "Matrix")
              for (i in seq(along = val)) {
                  dimnames(val[[i]]) = list(x@cnames[[i]], x@cnames[[i]])
                  val[[i]] = as(as(val[[i]], "pdmatrix"), "corrmatrix")
              }
              new("VarCorr",
                  scale = .Call("lmer_sigma", x, REML, PACKAGE = "Matrix"),
                  reSumry = val,
                  useScale = useScale)
          })

setMethod("gradient", signature(x = "lmer"),
          function(x, unconst, ...)
          .Call("lmer_gradient", x, unconst, PACKAGE = "Matrix"))

setMethod("summary", signature(object = "lmer"),
          function(object, ...)
          new("summary.lmer", object,
              showCorrelation = TRUE,
              useScale = !((object@family)$family %in% c("binomial", "poisson"))))

setMethod("show", signature(object = "lmer"),
          function(object)
          show(new("summary.lmer", object,
                   showCorrelation = FALSE,
                   useScale = !((object@family)$family %in% c("binomial", "poisson")))))
          
setMethod("show", "summary.lmer",
          function(object) {
              fcoef <- object@fixed
              useScale <- object@useScale
              corF <- as(as(vcov(object, useScale = useScale), "pdmatrix"),
                         "corrmatrix")
              DF <- getFixDF(object)
              coefs <- cbind(fcoef, corF@stdDev, DF)
              nc <- object@nc
              dimnames(coefs) <-
                  list(names(fcoef), c("Estimate", "Std. Error", "DF"))
                            digits <- max(3, getOption("digits") - 2)
              REML <- object@method == "REML"
              llik <- object@logLik
              dev <- object@deviance
              
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
              show(VarCorr(object, useScale = useScale))
              ngrps <- lapply(object@flist, function(x) length(levels(x)))
              cat(sprintf("# of obs: %d, groups: ", object@nc[length(object@nc)]))
              cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
              cat("\n")
              if (!useScale)
                  cat("\nEstimated scale (compare to 1) ",
                      .Call("lmer_sigma", object, FALSE, PACKAGE = "Matrix"),
                      "\n")
              if (nrow(coefs) > 0) {
                  if (useScale) {
                      stat <- coefs[,1]/coefs[,2]
                      pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
                      nms <- colnames(coefs)
                      coefs <- cbind(coefs, stat, pval)
                      colnames(coefs) <- c(nms, "t value", "Pr(>|t|)")
                  } else {
                      coefs <- coefs[, 1:2, drop = FALSE]
                      stat <- coefs[,1]/coefs[,2]
                      pval <- 2*pnorm(abs(stat), lower = FALSE)
                      nms <- colnames(coefs)
                      coefs <- cbind(coefs, stat, pval)
                      colnames(coefs) <- c(nms, "z value", "Pr(>|z|)")
                  }
                  cat("\nFixed effects:\n")
                  printCoefmat(coefs, tst.ind = 4, zap.ind = 3)
                  if (length(object@showCorrelation) > 0 && object@showCorrelation[1]) {
                      rn <- rownames(coefs)
                      dimnames(corF) <- list(
                                               abbreviate(rn, minlen=11),
                                               abbreviate(rn, minlen=6))
                      if (!is.null(corF)) {
                          p <- NCOL(corF)
                          if (p > 1) {
                              cat("\nCorrelation of Fixed Effects:\n")
                              corF <- format(round(corF, 3), nsmall = 3)
                              corF[!lower.tri(corF)] <- ""
                              print(corF[-1, -p, drop=FALSE], quote = FALSE)
                          }
                      }
                  }
              }
              invisible(object)
          })




## calculates degrees of freedom for fixed effects Wald tests
## This is a placeholder.  The answers are generally wrong.  It will
## be very tricky to decide what a 'right' answer should be with
## crossed random effects.

setMethod("getFixDF", signature(object="lmer"),
          function(object, ...)
      {
          nc <- object@nc[-seq(along = object@Omega)]
          p <- abs(nc[1]) - 1
          n <- nc[2]
          rep(n-p, p)
      })

setMethod("logLik", signature(object="mer"),
          function(object, REML = object@method == "REML", ...) {
              val <- -deviance(object, REML = REML)/2
              nc <- object@nc[-seq(a = object@Omega)]
              attr(val, "nall") <- attr(val, "nobs") <- nc[2]
              attr(val, "df") <- abs(nc[1]) +
                  length(.Call("lmer_coef", object, 0, PACKAGE = "Matrix"))
              attr(val, "REML") <- REML 
              class(val) <- "logLik"
              val
          })

setMethod("logLik", signature(object="lmer"),
          function(object, ...) object@logLik)

setMethod("anova", signature(object = "lmer"),
          function(object, ...)
      {
          mCall <- match.call(expand.dots = TRUE)
          dots <- list(...)
          modp <- logical(0)
          if (length(dots))
              modp <- sapply(dots, inherits, "lmer") | sapply(dots, inherits, "lm")
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
              foo@status["factored"] <- FALSE
              .Call("lmer_factor", foo, PACKAGE="Matrix")
              dfr <- getFixDF(foo)
              rcol <- ncol(foo@RXX)
              ss <- foo@RXX[ , rcol]^2
              ssr <- ss[[rcol]]
              ss <- ss[seq(along = dfr)]
              names(ss) <- object@cnames[[".fixed"]][seq(along = dfr)]
              asgn <- foo@assign
              terms <- foo@terms
              nmeffects <- attr(terms, "term.labels")
              if ("(Intercept)" %in% names(ss))
                  nmeffects <- c("(Intercept)", nmeffects)
              ss <- unlist(lapply(split(ss, asgn), sum))
              df <- unlist(lapply(split(asgn,  asgn), length))
              dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
              ms <- ss/df
              f <- ms/(ssr/dfr)
              P <- pf(f, df, dfr, lower.tail = FALSE)
              table <- data.frame(df, ss, ms, dfr, f, P)
              dimnames(table) <-
                  list(nmeffects,
                       c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
              if ("(Intercept)" %in% nmeffects) table <- table[-1,]
              attr(table, "heading") <- "Analysis of Variance Table"
              class(table) <- c("anova", "data.frame")
              table
          }
      })

setMethod("update", signature(object = "lmer"),
          function(object, formula., ..., evaluate = TRUE)
      {
          call <- object@call
          if (is.null(call))
              stop("need an object with call component")
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


setMethod("confint", signature(object = "lmer"),
          function (object, parm, level = 0.95, ...) 
      {
          cf <- fixef(object)
          pnames <- names(cf)
          if (missing(parm)) 
              parm <- seq(along = pnames)
          else if (is.character(parm)) 
              parm <- match(parm, pnames, nomatch = 0)
          a <- (1 - level)/2
          a <- c(a, 1 - a)
          pct <- paste(round(100 * a, 1), "%")
          ci <- array(NA, dim = c(length(parm), 2),
                      dimnames = list(pnames[parm], pct))
          ses <- sqrt(diag(vcov(object)))[parm]
          ci[] <- cf[parm] + ses * t(outer(a, getFixDF(object)[parm], qt))
          ci
      })

setMethod("deviance", "mer",
          function(object, REML = NULL, ...) {
              .Call("lmer_factor", object, PACKAGE = "Matrix")
              if (is.null(REML))
                  REML <- object@method == "REML"
              object@deviance[[ifelse(REML, "REML", "ML")]]
          })


setMethod("deviance", "lmer",
          function(object, ...) -2 * c(object@logLik))


setMethod("chol", signature(x = "lmer"),
          function(x, pivot = FALSE, LINPACK = pivot) {
              x@status["factored"] <- FALSE # force a decomposition
              .Call("lmer_factor", x, PACKAGE = "Matrix")
          })

setMethod("solve", signature(a = "lmer", b = "missing"),
          function(a, b, ...)
          .Call("lmer_invert", a, PACKAGE = "Matrix")
          )

setMethod("formula", "lmer", function(x, ...) x@call$formula)

setMethod("vcov", signature(object = "lmer"),
          function(object, REML = object@method == "REML", useScale = TRUE,...) {
              sc <- .Call("lmer_sigma", object, REML, PACKAGE = "Matrix")
              rr <- object@RXX
              nms <- object@cnames[[".fixed"]]
              dimnames(rr) <- list(nms, nms)
              nr <- nrow(rr)
              rr <- rr[-nr, -nr, drop = FALSE]
              rr <- rr %*% t(rr)
              if (useScale) {
                  rr = sc^2 * rr
              }
              rr
          })

## Extract the L matrix 
setAs("lmer", "dtTMatrix",
      function(from)
  {
      ## force a refactorization if the factors have been inverted
      if (from@status["inverted"]) from@status["factored"] <- FALSE
      .Call("lmer_factor", from, PACKAGE = "Matrix")
      L <- lapply(from@L, as, "dgTMatrix")
      nf <- length(from@D)
      Gp <- from@Gp
      nL <- Gp[nf + 1]
      Li <- integer(0)
      Lj <- integer(0)
      Lx <- double(0)
      for (i in 1:nf) {
          for (j in 1:i) {
              Lij <- L[[Lind(i, j)]]
              Li <- c(Li, Lij@i + Gp[i])
              Lj <- c(Lj, Lij@j + Gp[j])
              Lx <- c(Lx, Lij@x)
          }
      }
      new("dtTMatrix", Dim = as.integer(c(nL, nL)), i = Li, j = Lj, x = Lx,
          uplo = "L", diag = "U")
  })

## Extract the ZZX matrix
setAs("lmer", "dsTMatrix",
      function(from)
  {
      .Call("lmer_inflate", from, PACKAGE = "Matrix")
      ZZpO <- lapply(from@ZZpO, as, "dgTMatrix")
      ZZ <- lapply(from@ZtZ, as, "dgTMatrix")
      nf <- length(ZZpO)
      Gp <- from@Gp
      nZ <- Gp[nf + 1]
      Zi <- integer(0)
      Zj <- integer(0)
      Zx <- double(0)
      for (i in 1:nf) {
          ZZpOi <- ZZpO[[i]]
          Zi <- c(Zi, ZZpOi@i + Gp[i])
          Zj <- c(Zj, ZZpOi@j + Gp[i])
          Zx <- c(Zx, ZZpOi@x)
          if (i > 1) {
              for (j in 1:(i-1)) {
                  ZZij <- ZZ[[Lind(i, j)]]
                  ## off-diagonal blocks are transposed
                  Zi <- c(Zi, ZZij@j + Gp[j])
                  Zj <- c(Zj, ZZij@i + Gp[i])
                  Zx <- c(Zx, ZZij@x)
              }
          }
      }
      new("dsTMatrix", Dim = as.integer(c(nZ, nZ)), i = Zi, j = Zj, x = Zx,
          uplo = "U")
  })

setMethod("fitted", signature(object = "lmer"),
          function(object, ...)
          napredict(attr(object@frame, "na.action"), object@fitted))

setMethod("residuals", signature(object = "lmer"),
          function(object, ...)
          naresid(attr(object@frame, "na.action"), object@residuals))

setMethod("resid", signature(object = "lmer"),
          function(object, ...) do.call("residuals", c(list(object), list(...))))

setMethod("coef", signature(object = "lmer"),
          function(object, ...)
      {
          fef <- data.frame(rbind(object@fixed), check.names = FALSE)
          ref <- as(ranef(object), "list")
          names(ref) <- names(object@flist)
          val <- lapply(ref, function(x) fef[rep(1, nrow(x)),])
          for (i in seq(a = val)) {
              refi <- ref[[i]]
              row.names(val[[i]]) <- row.names(refi)
              if (!all(names(refi) %in% names(fef)))
                  stop("unable to align random and fixed effects")
              val[[i]][ , names(refi)] <- val[[i]][ , names(refi)] + refi
          }
          new("lmer.coef", val)
      })

setMethod("plot", signature(x = "lmer.coef"),
          function(x, y, ...)
      {
          if (require("lattice", quietly = TRUE)) {
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
          }
      })

setMethod("plot", signature(x = "lmer.ranef"),
          function(x, y, ...)
      {
          if (require("lattice", quietly = TRUE)) 
              lapply(x, function(x) {
                  cn <- lapply(colnames(x), as.name)
                  switch(min(ncol(x), 3),
                         qqmath(eval(substitute(~ x,
                                                list(x = cn[[1]]))),
                                x, ...),
                         xyplot(eval(substitute(y ~ x,
                                                list(y = cn[[1]],
                                                     x = cn[[2]]))),
                                x, ...),
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

setMethod("show", signature(object="VarCorr"),
          function(object)
      {
          digits <- max(3, getOption("digits") - 2)
          useScale <- length(object@useScale) > 0 && object@useScale[1]
          sc <- ifelse(useScale, object@scale,  1.)
          reStdDev <- c(lapply(object@reSumry,
                               function(x, sc)
                               sc*x@stdDev,
                               sc = sc), list(Residual = sc))
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
                          lapply(object@reSumry,
                                 function(x, maxlen) {
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
      })

glmmMCMC <- function(obj, method = c("full"), nsamp = 1)
{
    if (!inherits(obj, "lmer")) stop("obj must be of class `lmer'")
    if (obj@family$family == "gaussian" && obj@family$link == "identity")
        warn("glmmMCMC not indended for Gaussian family with identity link")
    cv <- Matrix:::lmerControl()
    family <- obj@family
    frm <- obj@frame
    fixed.form <- Matrix:::nobars(obj@call$formula)
    if (inherits(fixed.form, "name")) # RHS is empty - use a constant
        fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))
    glm.fit <- glm(eval(fixed.form), family, frm, x = TRUE,
                   y = TRUE)
    x <- glm.fit$x
    y <- as.double(glm.fit$y)
    bars <- Matrix:::findbars(obj@call$formula[[3]])
    random <-
        lapply(bars,
               function(x) list(model.matrix(eval(substitute(~term,
                                                             list(term=x[[2]]))),
                                             frm),
                                eval(substitute(as.factor(fac)[,drop = TRUE],
                                                list(fac = x[[3]])), frm)))
    names(random) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    if (any(names(random) != names(obj@flist)))
        random <- random[names(obj@flist)]
    mmats <- c(lapply(random, "[[", 1),
               .fixed = list(cbind(glm.fit$x, .response = glm.fit$y)))
    mer <- as(obj, "mer")

    eta <- glm.fit$linear.predictors # perhaps later change this to obj@fitted?
    wts <- glm.fit$prior.weights
    wtssqr <- wts * wts
    offset <- glm.fit$offset
    if (is.null(offset)) offset <- numeric(length(eta))
    off <- numeric(length(eta))
    mu <- numeric(length(eta))

    dev.resids <- quote(family$dev.resids(y, mu, wtssqr))
    linkinv <- quote(family$linkinv(eta))
    mu.eta <- quote(family$mu.eta(eta))
    variance <- quote(family$variance(mu))
    LMEopt <- getAnywhere("LMEoptimize<-")
    doLMEopt <- quote(LMEopt(x = mer, value = cv))

    GSpt <- .Call("glmer_init", environment())
    nf <- length(obj@flist)
    fixed <- obj@fixed
    varc <- .Call("lmer_coef", mer, 2, PACKAGE = "Matrix")
    b <- .Call("lmer_ranef", mer, PACKAGE = "Matrix")
    ans <- list(fixed = matrix(0, nr = length(fixed), nc = nsamp),
                varc = matrix(0, nr = length(varc), nc = nsamp))
    for (i in 1:nsamp) {
        ## conditional means and variances of fixed effects
        fupd <- .Call("glmer_fixed_update", GSpt, b, fixed, PACKAGE = "Matrix")
        ans$fixed[,i] <- fixed <- fupd$fixed
        ## sample from the conditional distribution of beta given b and y
        ## conditional means and variances of random_effects
        .Call("glmer_bhat", GSpt, fixed, varc, PACKAGE = "Matrix")
        print(bhat <- .Call("lmer_ranef", mer, PACKAGE = "Matrix"))
        ## sample from the conditional distribution of b given beta and y
        ## sample from the conditional distribution of varc given b
        ans$varc[,i] <- varc
    }
    ans
}

