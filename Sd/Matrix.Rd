\name{Matrix}
\alias{Matrix}
\title{
    Construct a Classed Matrix
}
\description{
    Construct a Matrix of class \code{Matrix}.
}
\usage{
Matrix(data=NA, nrow, ncol, byrow=F, dimnames=NULL, tol)
}
\arguments{
    \item{data}{
	vector of data values for the matrix.
    }
    \item{nrow,ncol}{
	The desired number of rows / columns for the matrix. If both
	dimensions are left unspecified, the default is a column vector
	if \code{byrow = F} and a row vector if \code{byrow = T}. 
	If only one of \code{nrow} or \code{ncol} is missing,
	the other is chosen according to the length of the data.
    }
    \item{byrow}{
	Data values are assigned to the matrix by row if \code{byrow =
	    T} and otherwise by column.
    }
    \item{dimnames}{
	A list of length 2 giving a dimnames attribute for the
	matrix. Each element of the list must be either \code{NULL} or
	else a vector of character strings of length equal to the
	corresponding element of the \code{dim} attribute of the result.
    }
}
\value{
    Returns an \code{nrow} by \code{ncol} matrix of class \code{"Matrix"}. 
}
\details{
    The input data must be a vector.  If you specify both \code{nrow}
    and \code{ncol}, you can supply a single value for \code{data} and
    that value will be used to fill the matrix.

    Although it is sometime possible to mix unclassed matrices (created
    with \code{matrix}) with ones of class \code{"Matrix"}, it is much
    safer to always use carefully constructed ones of class
    \code{"Matrix"}.
}
\seealso{
    \code{\link{as.Matrix}}, \code{\link{Matrix.class}}, \code{\link{matrix}}.
}
\examples{
Matrix(0, 3, 2)             # a 3 by 2 matrix of zeros

Matrix(1:6, 3, 2)           # a 3 by 2 matrix
Matrix(1:6, nrow=3)   
Matrix(1:6, ncol=2)   
Matrix(1:9, nrow=3, dimnames =list(c("a", "b", "c"), c("A", "B", "C")))

}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
