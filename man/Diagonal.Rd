\name{Diagonal}
\alias{Diagonal}
\title{Create Diagonal Matrix Object}
\description{
  Create a diagonal matrix object, i.e., an object inheriting from
  \code{\link{diagonalMatrix-class}}.
}
\usage{
Diagonal(n, x = NULL)
}
\arguments{
  \item{n}{integer specifying the dimension of the (square) matrix.  If
    missing, \code{length(x)} is used.}
  \item{x}{numeric or logical; if missing, a \emph{unit} diagonal
    \eqn{n \times n}{n x n} matrix is created.}
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
\value{
  an object of class
  \code{\link{ddiMatrix-class}} or \code{\link{ldiMatrix-class}}
  (with \dQuote{superclass} \code{\link{diagonalMatrix-class}}).
}
%\author{Martin Maechler}
\seealso{the generic function \code{\link{diag}} for \emph{extraction}
  of the diagonal from a matrix;
  \code{\link{Matrix}}, \code{\link{diagonalMatrix-class}}
}
\examples{
Diagonal(3)
Diagonal(x = 10^(3:1))
Diagonal(x = (1:4) >= 2)#-> "ldiMatrix"
}
\keyword{array}
\keyword{algebra}
