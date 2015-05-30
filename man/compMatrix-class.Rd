\name{compMatrix-class}
\docType{class}
\title{Class "compMatrix" of Composite (Factorizable) Matrices}
\alias{compMatrix-class}
\alias{dimnames<-,compMatrix,list-method}
\alias{dimnames<-,compMatrix,NULL-method}
\description{
  Virtual class of \emph{composite} matrices; i.e., matrices
  that can be \emph{factorized}, typically as a product of simpler
  matrices.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{factors}:}{Object of class \code{"list"} - a list of
      factorizations of the matrix.  Note that this is typically empty,
      i.e., \code{list()}, initially and is \emph{updated
	\bold{automagically}} whenever a matrix factorization is computed.}
    \item{\code{Dim}, \code{Dimnames}:}{inherited from the
      \code{\linkS4class{Matrix}} class, see there.}
  }
}
\section{Extends}{
  Class \code{"Matrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{dimnames<-}{\code{signature(x = "compMatrix", value = "list")}: set
      the \code{dimnames} to a \code{\link{list}} of length 2, see
      \code{\link{dimnames<-}}.  The \code{factors} slot is currently reset to
      empty, as the factorization \code{dimnames} would have to be adapted, too.}
  }
}
\seealso{
  The matrix factorization classes
  \code{"\linkS4class{MatrixFactorization}"} and their generators,
  \code{\link{lu}()},
  \code{\link{qr}()},
  \code{\link{chol}()} and \code{\link{Cholesky}()},
  \code{\link{BunchKaufman}()}, \code{\link{Schur}()}.
}
%% FIXME:  add:
% \examples{
%
% }
\keyword{classes}
