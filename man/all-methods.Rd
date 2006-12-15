\name{all-methods}
\docType{methods}
\alias{all-methods}
\alias{all,Matrix-method}
\alias{any,Matrix-method}
% \alias{all,lMatrix-method}
% \alias{any,lMatrix-method}
\alias{all,ldenseMatrix-method}
\alias{any,ldenseMatrix-method}
\alias{all,lsparseMatrix-method}
\alias{any,lsparseMatrix-method}
\alias{all,lsyMatrix-method}
%
\title{"Matrix" Methods for Functions all() and any()}
\description{
  The basic \R functions \code{\link{all}} and \code{\link{any}} now
  have methods for \code{\linkS4class{Matrix}} objects and should
  behave as for \code{\link{matrix}} ones.
}
\section{Methods}{
%% FIXME: write more
  \describe{
    \item{all}{\code{signature(x = "Matrix")}: ...}
    \item{any}{\code{signature(x = "Matrix")}: ...}

    \item{all}{\code{signature(x = "ldenseMatrix")}: ...}
    \item{all}{\code{signature(x = "lsparseMatrix")}: ...}

  }
}
\keyword{methods}
\examples{
M <- Matrix(1:12 +0, 3,4)
all(M >= 1) # TRUE
any(M < 0 ) # FALSE
% M[2,3] <- NA
% OOPS!  FIXME!! in 'Ops'
% all(M >= 0) # NA
% any(M >= 0) # NA
}
