\name{sparseQR-class}
\docType{class}
\title{Sparse QR decomposition of a sparse matrix}
\alias{sparseQR-class}
\alias{qr.R,sparseQR-method}
\alias{qr.Q}% the generic
\alias{qr.Q,sparseQR-method}
\alias{qr.coef,sparseQR,ddenseMatrix-method}
\alias{qr.coef,sparseQR,matrix-method}
\alias{qr.coef,sparseQR,Matrix-method}
\alias{qr.coef,sparseQR,numeric-method}
\alias{qr.fitted,sparseQR,ddenseMatrix-method}
\alias{qr.fitted,sparseQR,matrix-method}
\alias{qr.fitted,sparseQR,Matrix-method}
\alias{qr.fitted,sparseQR,numeric-method}
\alias{qr.qty,sparseQR,ddenseMatrix-method}
\alias{qr.qty,sparseQR,matrix-method}
\alias{qr.qty,sparseQR,numeric-method}
\alias{qr.qty,sparseQR,Matrix-method}
\alias{qr.qy,sparseQR,ddenseMatrix-method}
\alias{qr.qy,sparseQR,matrix-method}
\alias{qr.qy,sparseQR,numeric-method}
\alias{qr.qy,sparseQR,Matrix-method}
\alias{qr.resid,sparseQR,ddenseMatrix-method}
\alias{qr.resid,sparseQR,matrix-method}
\alias{qr.resid,sparseQR,Matrix-method}
\alias{qr.resid,sparseQR,numeric-method}
%\alias{solve,sparseQR,ANY-method}--> solve-methods.Rd
%
\description{Objects class \code{"sparseQR"} represent a QR decomposition of a
  sparse \eqn{m \times n}{m x n} (\dQuote{long}: \eqn{m \ge n}{m >= n})
  rectangular matrix \eqn{A}, typically resulting from
  \code{\link{qr}()}, see \sQuote{Details} notably about row and column
  permutations for pivoting.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("sparseQR", ...)}
  but are more commonly created by function \code{\link[base]{qr}} applied
  to a sparse matrix such as a matrix of class
  \code{\linkS4class{dgCMatrix}}.
}
\section{Slots}{
  \describe{
    \item{\code{V}:}{Object of class \code{"dgCMatrix"}.  The columns of
      \code{V} are the vectors that generate the Householder
      transformations of which the matrix Q is composed.}
    \item{\code{beta}:}{Object of class \code{"numeric"}, the normalizing
      factors for the Householder transformations.}
    \item{\code{p}:}{Object of class \code{"integer"}: Permutation
      (of \code{0:(n-1)}) applied to the rows of the original matrix.}
    \item{\code{R}:}{Object of class \code{"dgCMatrix"}: An upper
      triangular matrix of the same dimension as \eqn{X}.}
    \item{\code{q}:}{Object of class \code{"integer"}: Permutation
      applied from the right, i.e., to the \emph{columns} of the
      original matrix.  Can be of length 0 which implies no permutation.}
  }
}
\section{Methods}{
  \describe{
    \item{qr.R}{\code{signature(qr = "sparseQR")}: compute the upper
      triangular \eqn{R} matrix of the QR decomposition.
      Note that this currently warns because of possible permutation
      mismatch with the classical \code{qr.R()} result, \emph{and} you
      can suppress these warnings by setting \code{\link{options}()}
      either \code{"Matrix.quiet.qr.R"} or (the more general)
      either \code{"Matrix.quiet"} to \code{\link{TRUE}}.}
    \item{qr.Q}{\code{signature(qr = "sparseQR")}: compute the
      orthogonal \eqn{Q} matrix of the QR decomposition.}
    \item{qr.coef}{\code{signature(qr = "sparseQR", y = "ddenseMatrix")}: ... }
    \item{qr.coef}{\code{signature(qr = "sparseQR", y = "matrix")}: ... }
    \item{qr.coef}{\code{signature(qr = "sparseQR", y = "numeric")}: ... }
    \item{qr.fitted}{\code{signature(qr = "sparseQR", y = "ddenseMatrix")}: ... }
    \item{qr.fitted}{\code{signature(qr = "sparseQR", y = "matrix")}: ... }
    \item{qr.fitted}{\code{signature(qr = "sparseQR", y = "numeric")}: ... }
    \item{qr.qty}{\code{signature(qr = "sparseQR", y = "ddenseMatrix")}: ... }
    \item{qr.qty}{\code{signature(qr = "sparseQR", y = "matrix")}: ... }
    \item{qr.qty}{\code{signature(qr = "sparseQR", y = "numeric")}: ... }
    \item{qr.qy}{\code{signature(qr = "sparseQR", y = "ddenseMatrix")}: ... }
    \item{qr.qy}{\code{signature(qr = "sparseQR", y = "matrix")}: ... }
    \item{qr.qy}{\code{signature(qr = "sparseQR", y = "numeric")}: ... }
    \item{qr.resid}{\code{signature(qr = "sparseQR", y = "ddenseMatrix")}: ... }
    \item{qr.resid}{\code{signature(qr = "sparseQR", y = "matrix")}: ... }
    \item{qr.resid}{\code{signature(qr = "sparseQR", y = "numeric")}: ... }

    \item{solve}{\code{signature(a = "sparseQR", b = "ANY")}: For
      \code{solve(a,b)}, simply uses \code{qr.coef(a,b)}.}
  }
}

%\references{}
%\author{}
\details{
  For a sparse \eqn{m \times n}{m x n} (\dQuote{long}: \eqn{m \ge n}{m >= n})
  rectangular matrix \eqn{A}, the sparse QR decomposition is either 
  \cr of the form \eqn{P A = Q R} with a (row)
  permutation matrix \eqn{P}, (encoded in the \code{p} slot of the
  result) if the \code{q} slot is of length 0,
  \cr or of the form \eqn{P A P* = Q R} with an extra (column) permutation
  matrix \eqn{P*} (encoded in the \code{q} slot).
  Note that the row permutation \eqn{P A} in \R is simply \code{A[p+1, ]}
  where \code{p} is the \code{p}-slot, a 0-based permutation of
  \code{1:m} applied to the rows of the original matrix.

  If the \code{q} slot has length \code{n} it is a 0-based permutation
  of \code{1:n} applied to the columns of the original matrix to reduce
  the amount of \dQuote{fill-in} in the matrix \eqn{R}, and
  \eqn{A P*} in \R is simply \code{A[ , q+1]}.

  \eqn{R} is an \eqn{m\times n}{m by n} matrix that is zero below the
  main diagonal, i.e., upper triangular (\eqn{m\times m}{m by n}) with
  \eqn{m-n} extra zero rows.

  The matrix \eqn{Q} is a "virtual matrix".  It is the product of
  \eqn{n} Householder transformations.  The information to generate
  these Householder transformations is stored in the \code{V} and
  \code{beta} slots.
  \cr
  Note however that \code{qr.Q()} returns the row permuted matrix
  \eqn{Q* := P^{-1}Q = P'Q}{Q* := P^(-1) Q = P'Q} as permutation matrices are
  orthogonal; and \eqn{Q*} is orthogonal itself because \eqn{Q} and \eqn{P} are.
  This is useful because then, as in the dense matrix and \pkg{base} \R
  matrix \code{\link{qr}} case, we have the mathematical identity
  \deqn{P A = Q* R,} in \R as \preformatted{            A[p+1,] == qr.Q(*) \%*\% R .}

  The \code{"sparseQR"} methods for the \code{qr.*} functions return
  objects of class \code{"dgeMatrix"} (see
  \code{\linkS4class{dgeMatrix}}).  Results from \code{qr.coef},
  \code{qr.resid} and \code{qr.fitted} (when \code{k == ncol(R)}) are
  well-defined and should match those from the corresponding dense matrix
  calculations.  However, because the matrix \code{Q} is not uniquely
  defined, the results of \code{qr.qy} and \code{qr.qty} do not
  necessarily match those from the corresponding dense matrix
  calculations.

  Also, the results of \code{qr.qy} and \code{qr.qty} apply to the
  permuted column order when the \code{q} slot has length \code{n}.
}
\seealso{
  \code{\link[base]{qr}}, \code{\link{qr.Q}},
  \code{\link{qr.R}},     \code{\link{qr.fitted}},
  \code{\link{qr.resid}}, \code{\link{qr.coef}},
  \code{\link{qr.qty}},   \code{\link{qr.qy}},

  Permutation matrices in the \pkg{Matrix} package: \code{\linkS4class{pMatrix}};
  \code{\linkS4class{dgCMatrix}}, \code{\linkS4class{dgeMatrix}}.

%%% not (yet) in CRAN-released version of Matrix:
%%   \code{\linkS4class{SPQR}} an alternative QR decomposition resulting
%%   from \code{\link{spqr}()}.
}
\examples{
data(KNex)
mm <- KNex $ mm
 y <- KNex $  y
 y. <- as(as.matrix(y), "dgCMatrix")
str(qrm <- qr(mm))
 qc  <- qr.coef  (qrm, y); qc. <- qr.coef  (qrm, y.) # 2nd failed in Matrix <= 1.1-0
 qf  <- qr.fitted(qrm, y); qf. <- qr.fitted(qrm, y.)
 qs  <- qr.resid (qrm, y); qs. <- qr.resid (qrm, y.)
stopifnot(all.equal(qc, as.numeric(qc.),  tolerance=1e-12),
          all.equal(qf, as.numeric(qf.),  tolerance=1e-12),
          all.equal(qs, as.numeric(qs.),  tolerance=1e-12),
          all.equal(qf+qs, y, tolerance=1e-12))
%% FIXME: add more similarly to ./lu.Rd, see also ./qr-methods.Rd
}
\keyword{classes}
\keyword{algebra}
\keyword{array}
