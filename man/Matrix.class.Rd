\name{Matrix.class}
\title{Determine the subclass of a Matrix}
\usage{
Matrix.class(x, tol = 0, symmetry = T, unit.diagonal = T, 
             triangularity = c(T, T), orthogonality = c(T, T), 
             normality = c(T, T))
}
\alias{Matrix.class}
\description{
    Determines the class of a Matrix in the \code{"Matrix"}
    class structure.
}
\arguments{
    \item{x}{
	a numeric or complex matrix.
    }
    \item{tol}{
	tolerance for the various tests. The default is 0.
    }
    \item{symmetry}{
	a logical value indicating whether or not to test for symmetry
	(conjugate symmetry in the complex case).
	The default is to test for symmetry.
    }
    \item{unit.diagonal}{
	a logical value indicating whether or not to test for a unit diagonal.
	The default is to test for a unit diagonal.
    }
    \item{triangularity}{
	a logical vector of length indicating whether or not to test
	lower and/or upper triangularity, respectively.
	The default is to test for both lower and upper triangularity.
    }
    \item{orthogonality}{
	a logical vector of length indicating whether or not to test row
	and/or column orthogonality, respectively.
	The default is to test for both row and column orthogonality.
    }
    \item{normality}{
	a logical vector of length indicating whether or not to see if the 
	rows and/or columns have norm one, respectively.
	The default is to test if both row and columns are normalized.
    }
}
\value{
    Returns a class for the matrix, relative to the \code{"Matrix"}
    class structure.
    The matrix will be of class \code{"Matrix"}, but there may also subclasses 
    depending on the options chosen and their results.
}
\author{Douglas Bates \email{bates@stat.wisc.edu},
    Saikat DebRoy \email{saikat@stat.wisc.edu}}
\seealso{
    \code{\link{is.Hermitian}}, \code{\link{is.LowerTriangular}},
    \code{\link{is.Orthonormal}}.
}
\examples{
x <- hilbert(9)
Matrix.class(x) # should be "Hermitian", "Matrix"
x <- Matrix( rnorm(9), 3, 3)
Matrix.class(x) # should be "Matrix"
x[row(x) < col(x)] <- 0
Matrix.class(x) # should be "LowerTriangular" "Matrix"
}
\keyword{array}
\keyword{algebra}
