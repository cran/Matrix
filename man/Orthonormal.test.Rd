\name{Orthonormal.test}
\title{
    Check for Orthogonality or Orthonormality
}
\usage{
is.Orthonormal(x, tol = sqrt(.Machine$double.eps), byrow = FALSE)
is.ColOrthonormal(x, tol = sqrt(.Machine$double.eps))
is.RowOrthonormal(x, tol = sqrt(.Machine$double.eps))
Orthogonal.test(x, byrow = FALSE, normal = TRUE)
Orthonormal.test(x, byrow = FALSE)
}
\alias{is.Orthonormal}
\alias{is.ColOrthonormal}
\alias{is.RowOrthonormal}
\alias{Orthonormal.test}
\alias{Orthogonal.test}
\description{
    Check for column (or row) orthogonality (or orthonormality) in a
    matrix or return the numerical value used to perform this check.
}
\arguments{
    \item{x}{
	an numeric or complex matrix.
    }
    \item{tol}{
	tolerance in the test for orthonormality or orthogonality.
    }
    \item{byrow}{logical. If \code{FALSE} (the default) the columns are
	tested, otherwise the rows are tested.
    }
    \item{normal}{logical. If \code{TRUE} (the default) determine the
	maximum modulus of the deviation of the squared lengths of the
	columns (rows) from unity.  Otherwise, determine only the
	maximum modulus of the inner products of distinct columns (rows).
    }
}
\value{
    \code{Orthonormal.test} returns a numeric measure of the
    deviation of the columns (rows) of the matrix from orthogonality,
    when \code{normal} is \code{FALSE}, or orthonormality, when
    \code{normal} is \code{TRUE}.  This value is always at least the
    maximum modulus of the inner products of distinct columns (rows).
    When \code{normal} it \code{TRUE}, it is the maximum of this value
    and the maximum modulus of the deviation of the squared lengths of the
    columns (rows) from unity.  That is, when \code{normal} is
    \code{TRUE} the value the maximum modulus
    of the deviation of \code{crossprod(x)} from an identity matrix.
	
    \code{is.Orthonormal}, \code{is.ColOrthonormal} and
    \code{is.RowOrthonormal} return a logical value according to whether
    the columns (rows) of the matrix are orthonormal.  These simply
    compare a result from \code{Orthonormal.test} with appropriate
    options to \code{tol}.
%    If \code{tol} is missing from \code{is.Orthonormal}, this value is
%    equivalent to \code{inherits(x, "Orthonormal")}.
}
\seealso{
    \code{\link{.Machine}}.
}
\examples{
v <- rnorm(3)
v <- v / norm(v, 'F')
x <- diag(3) - 2 * (v \%*\% t(v))
is.Orthonormal(x)
# should be TRUE
is.Orthonormal(x, tol = .Machine$double.eps)
# should be FALSE
is.Orthonormal(x, tol = sqrt(.Machine$double.eps))
# should be TRUE
Orthonormal.test(x)
}
\keyword{algebra}
