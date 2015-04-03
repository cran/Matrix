\name{rankMatrix}
\alias{rankMatrix}
\title{Rank of a Matrix}
\description{
  Compute \sQuote{the} matrix rank, a well-defined functional in theory,
  somewhat ambigous in practice.  We provide several methods, the
  default corresponding to Matlab's definition.
}
\usage{
rankMatrix(x, tol = NULL,
           method = c("tolNorm2", "qr.R", "qrLINPACK", "qr",
                      "useGrad", "maybeGrad"),
           sval = svd(x, 0, 0)$d, warn.t = TRUE)
}
\arguments{
  \item{x}{numeric matrix, of dimension \eqn{n \times m}{n x m}, say.}
  \item{tol}{nonnegative number specifying a (relative,
    \dQuote{scalefree}) tolerance for testing of
    \dQuote{practically zero} with specific meaning depending on
    \code{method}; by default, \code{max(dim(x)) * \link{.Machine}$double.eps}
    is according to Matlab's default (for its only method which is our
    \code{method="tolNorm2"}).}
  \item{method}{a character string specifying the computational method
    for the rank, can be abbreviated:
    \describe{
      \item{\code{"tolNorm2"}:}{the number of singular values
	\code{>= tol * max(sval)};}
      \item{\code{"qrLINPACK"}:}{for a dense matrix, this is the rank of
	\code{\link[base]{qr}(x, tol, LAPACK=FALSE)} (which is
	\code{qr(...)$rank});
	\cr
	This ("qr*", dense) version used to be \emph{the} recommended way to
	compute a matrix rank for a while in the past.

	For sparse \code{x}, this is equivalent to \code{"qr.R"}.
      }
      \item{\code{"qr.R"}:}{this is the rank of triangular matrix
	\eqn{R}, where \code{qr()} uses LAPACK or a "sparseQR" method
	(see \code{\link{qr-methods}}) to compute the decomposition
	\eqn{QR}.  The rank of \eqn{R} is then defined as the number of
	\dQuote{non-zero} diagonal entries \eqn{d_i} of \eqn{R}, and
	\dQuote{non-zero}s fulfill
	\eqn{|d_i| \ge \mathrm{tol}\cdot\max(|d_i|)}{|d_i| >= tol * max(|d_i|)}.
      }
      \item{\code{"qr"}:}{is for back compatibility; for dense \code{x},
	it corresponds to \code{"qrLINPACK"}, whereas for sparse
	\code{x}, it uses \code{"qr.R"}.

	For all the "qr*" methods, singular values \code{sval} are not
	used, which may be crucially important for a large sparse matrix
	\code{x}, as in that case, when \code{sval} is not specified,
	the default, computing \code{\link{svd}()} currently coerces
	\code{x} to a dense matrix.
      }

      \item{\code{"useGrad"}:}{considering the \dQuote{gradient} of the
	(decreasing) singular values, the index of the \emph{smallest} gap.}
      \item{\code{"maybeGrad"}:}{choosing method \code{"useGrad"} only when
	that seems \emph{reasonable}; otherwise using \code{"tolNorm2"}.}

%% FIXME say more
    }
  }
  \item{sval}{numeric vector of non-increasing singular values of
    \code{x}; typically unspecified and computed from \code{x} when
    needed, i.e., unless \code{method = "qr"}.}
  \item{warn.t}{logical indicating if \code{rankMatrix()} should warn
    when it needs \code{\link{t}(x)} instead of \code{x}.  Currently, for
    \code{method = "qr"} only, gives a warning by default because the
    caller often could have passed \code{t(x)} directly, more efficiently.}
}
% \details{
%   FIXME
% }
\note{For large sparse matrices \code{x}, unless you can specify
  \code{sval} yourself, currently \code{method = "qr"} may
  be the only feasible one, as the others need \code{sval} and call
  \code{\link{svd}()} which currently coerces \code{x} to a
  \code{\linkS4class{denseMatrix}} which may be very slow or impossible,
  depending on the matrix dimensions.

  Note that in the case of sparse \code{x}, \code{method = "qr"}, all
  non-strictly zero diagonal entries \eqn{d_i} where counted, up to
  including \pkg{Matrix} version 1.1-0, i.e., that method implicitly
  used \code{tol = 0}, see also the seed(42) example below.
}
\value{
  positive integer in \code{1:min(dim(x))}, with attributes detailing
  the method used.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{Martin Maechler; for the "*Grad" methods, building on
  suggestions by Ravi Varadhan.
}
\seealso{
 \code{\link{qr}}, \code{\link{svd}}.
}
\examples{
rankMatrix(cbind(1, 0, 1:3)) # 2

(meths <- eval(formals(rankMatrix)$method))

## a "border" case:
H12 <- Hilbert(12)
rankMatrix(H12, tol = 1e-20) # 12;  but  11  with default method & tol.
sapply(meths, function(.m.) rankMatrix(H12, method = .m.))
## tolNorm2    qr   qr.R  qrLINPACK  useGrad maybeGrad
##       11    12     11         12       11        11
## The meaning of 'tol' for method="qrLINPACK" and *dense* x is not entirely "scale free"
rMQL <- function(ex, M) rankMatrix(M, method="qrLINPACK",tol = 10^-ex)
rMQR <- function(ex, M) rankMatrix(M, method="qr.R",     tol = 10^-ex)
sapply(5:15, rMQL, M = H12) # result is platform dependent
##  7  7  8 10 10 11 11 11 12 12 12  {x86_64}
sapply(5:15, rMQL, M = 1000 * H12) # not identical unfortunately
##  7  7  8 10 11 11 12 12 12 12 12
sapply(5:15, rMQR, M = H12)
##  5  6  7  8  8  9  9 10 10 11 11
sapply(5:15, rMQR, M = 1000 * H12) # the *same*
\dontshow{
  (r12 <- sapply(5:15, rMQR, M = H12))
  stopifnot(identical(r12, sapply(5:15, rMQR, M = H12 / 100)),
            identical(r12, sapply(5:15, rMQR, M = H12 * 1e5)))

  rM1 <- function(ex, M) rankMatrix(M, tol = 10^-ex)
  (r12 <- sapply(5:15, rM1, M = H12))
  stopifnot(identical(r12, sapply(5:15, rM1, M = H12 / 100)),
            identical(r12, sapply(5:15, rM1, M = H12 * 1e5)))
}

## "sparse" case:
M15 <- kronecker(diag(x=c(100,1,10)), Hilbert(5))
sapply(meths, function(.m.) rankMatrix(M15, method = .m.))
#--> all 15, but 'useGrad' has 14.

## "large" sparse
n <- 250000; p <- 33; nnz <- 10000
L <- sparseMatrix(i = sample.int(n, nnz, replace=TRUE),
                  j = sample.int(p, nnz, replace=TRUE), x = rnorm(nnz))
(st1 <- system.time(r1 <- rankMatrix(L)))                #  warning +  ~1.5 sec
(st2 <- system.time(r2 <- rankMatrix(L, method = "qr"))) # considerably faster!
r1[[1]] == print(r2[[1]]) ## -->  ( 33  TRUE )
\dontshow{ stopifnot(r1[[1]] == 33, 33 == r2[[1]])
 if(doExtras <- interactive() ||
	identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras"))))
  stopifnot(st2[[1]] < 0.2) # seeing 0.03 (on ~ 2010-hardware; R 3.0.2)
}
## another sparse-"qr" one, which ``failed'' till 2013-11-23:
set.seed(42)
f1 <- factor(sample(50, 1000, replace=TRUE))
f2 <- factor(sample(50, 1000, replace=TRUE))
f3 <- factor(sample(50, 1000, replace=TRUE))
rbind. <- if(getRversion() < "3.2.0") rBind else rbind
D <- t(do.call(rbind., lapply(list(f1,f2,f3), as, 'sparseMatrix')))
dim(D); nnzero(D) ## 1000 x 150 // 3000 non-zeros (= 2\%)
stopifnot(rankMatrix(D,           method='qr') == 148,
	  rankMatrix(crossprod(D),method='qr') == 148)
}
\keyword{algebra}
\keyword{array}
