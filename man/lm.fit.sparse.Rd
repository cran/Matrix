\name{lm.fit.sparse}
\alias{lm.fit.sparse}
\title{Fitter Function for Sparse Linear Models}
\description{
  The basic computing engine for sparse linear least squares regression.

  Note that the exact interface (arguments, return value) currently is
  \bold{experimental}, and is bound to change.  Use at your own risk!
}
\usage{
lm.fit.sparse(x, y, w = NULL, offset = NULL,
              method = c("qr", "cholesky"),
              tol = 1e-7, singular.ok = TRUE, order = NULL,
              transpose = FALSE)
}
\arguments{
  \item{x}{\emph{sparse} design matrix of dimension \code{n * p}, i.e.,
    an \R object of a \code{\link{class}} extending
    \code{\linkS4class{dsparseMatrix}}; typically the result of
    \code{\link{sparse.model.matrix}}.}
  \item{y}{vector of observations of length \code{n}, or a matrix with
    \code{n} rows.}
  \item{w}{vector of weights (length \code{n}) to be used in the fitting
    process.  Weighted least squares is used with weights \code{w},
    i.e., \code{sum(w * e^2)} is minimized.

    \bold{Not yet implemented !}
  }
  \item{offset}{numeric of length \code{n}).  This can be used to
    specify an \emph{a priori} known component to be included in the
    linear predictor during fitting.}

  \item{method}{a character string specifying the (factorization)
    method.  Currently, \code{"qr"} or \code{"cholesky"}.}
  \item{tol}{[for back-compatibility only; unused:] tolerance for the
    \code{\link{qr}} decomposition.  Default is 1e-7.}
  \item{singular.ok}{[for back-compatibility only; unused:] logical. If
    \code{FALSE}, a singular model is an error.}
  \item{order}{integer or \code{NULL}, for \code{method == "qr"}, will
    determine how the fill-reducing ordering (aka permutation) for the
    \dQuote{symbolic} part is determined (in \code{cs_amd()}), with the
    options {0: natural}, {1: Chol}, {2: LU}, and {3: QR}, where
    \code{3} is the default.}
  \item{transpose}{
    logical; if true, use the transposed matrix \code{t(x)} instead of
    \code{x}.
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  Either a single numeric vector or a list of four numeric vectors.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\seealso{
  \code{\link{sparse.model.matrix}}; the non-sparse function in standard
    \R's package \pkg{stats}: \code{\link{lm.fit}()}.
}
\examples{
dd <- expand.grid(a = as.factor(1:3),
                  b = as.factor(1:4),
                  c = as.factor(1:2),
                  d= as.factor(1:8))
n <- nrow(dd <- dd[rep(seq_len(nrow(dd)), each = 10), ])
set.seed(17)
dM <- cbind(dd, x = round(rnorm(n), 1))
## randomly drop some
n <- nrow(dM <- dM[- sample(n, 50),])
dM <- within(dM, { A <- c(2,5,10)[a]
                   B <- c(-10,-1, 3:4)[b]
                   C <- c(-8,8)[c]
                   D <- c(10*(-5:-2), 20*c(0, 3:5))[d]
   Y <- A + B + A*B + C + D + A*D + C*x + rnorm(n)/10
   wts <- sample(1:10, n, replace=TRUE)
   rm(A,B,C,D)
})
str(dM) # 1870 x 7

X <- sparse.model.matrix( ~ (a+b+c+d)^2 + c*x, data = dM)
dim(X) # 1870 x 69
X[1:10, 1:20]

## For now, use  'Matrix:::'  --- TODO : export once interface is clear!

Xd <- as(X,"matrix")
system.time(fmDense <- lm.fit(Xd, y = dM[,"Y"]))
system.time( r1 <- Matrix:::lm.fit.sparse(X, y = dM[,"Y"]) ) # *is* faster
stopifnot(all.equal(r1, unname(fmDense$coeff), tol = 1e-12))
system.time(
     r2 <- Matrix:::lm.fit.sparse(X, y = dM[,"Y"], method = "chol") )
stopifnot(all.equal(r1, r2$coef, tol = 1e-12),
          all.equal(fmDense$residuals, r2$residuals, tol=1e-9)
         )
## with weights:
system.time(fmD.w <- with(dM, lm.wfit(Xd, Y, w = wts)))
system.time(fm.w1 <- with(dM, Matrix:::lm.fit.sparse(X, Y, w = wts)))
system.time(fm.w2 <- with(dM, Matrix:::lm.fit.sparse(X, Y, w = wts,
                                                     method = "chol") ))
stopifnot(all.equal(fm.w1, unname(fmD.w$coeff), tol = 1e-12),
          all.equal(fm.w2$coef, fm.w1, tol = 1e-12),
          all.equal(fmD.w$residuals, fm.w2$residuals, tol=1e-9)
          )
}
\keyword{regression}
\keyword{array}
