\name{Matrix.class}
\alias{Matrix.class}
\title{
Subclass Determination for Matrices.
}
\description{
Determines the class of a Matrix relative to the \code{"Matrix"} class structure.
}
\usage{
Matrix.class(x, tol = 0, symmetry = T, unit.diagonal = T, 
             triangularity = c(T, T), orthogonality = c(T, T), 
             normality = c(T, T))
}
\arguments{
\item{x}{
a numeric or complex matrix.
}
\item{tol}{
tolerance for the various tests. The default is 0.
}
\item{symmetry}{
a logical value indicating whether or not to test for symmetry (conjugate symmetry
in the complex case). The default is to test for symmetry.
}
\item{unit.diagonal}{
a logical value indicating whether or not to test for a unit diagonal.
The default is to test for a unit diagonal.
}
\item{triangularity}{
a logical vector of length indicating whether or not to test lower and/or upper
triangularity, respectively.
The default is to test for both lower and upper triangularity.
}
\item{orthogonality}{
a logical vector of length indicating whether or not to test row and/or column
orthogonality, respectively.
The default is to test for both row and column orthogonality.
}
\item{normality}{
a logical vector of length indicating whether or not to see if the 
rows and/or columns have norm one, respectively.
The default is to test if both row and columns are normalized.
}}
\value{
Returns a class for the matrix, relative to the \code{"Matrix"} class structure.
The matrix will be of class \code{"Matrix"}, but there may also subclasses 
depending on the options chosen and their results.
}
\references{
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{is.Hermitian}}, \code{\link{is.LowerTriangular}}, \code{\link{is.Orthonormal}}.
}
\examples{
x <- Matrix( rnorm(9), 3, 3)
Matrix.class(x)
v <- rnorm(3)
v <- v /vecnorm(v)
x <- Diagonal(3) - 2 * (v \%*\% t(v))
Matrix.class(x)
Matrix.class(x, tol = sqrt(.Machine$double.eps))
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
