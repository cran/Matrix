\name{norm}
\title{Matrix Norms}
\alias{norm}
\alias{norm,ANY,missing-method}
\alias{norm,matrix,character-method}
\alias{norm,Matrix,character-method}
\usage{
norm(x, type, \dots)
}
\description{
  Computes a matrix norm of \code{x}, using Lapack for dense matrices.
  The norm can be the one  (\code{"O"}, or \code{"1"}) norm, the
  infinity (\code{"I"}) norm, the Frobenius  (\code{"F"}) norm,
  the maximum modulus (\code{"M"}) among elements of a matrix, or the
  spectral norm or 2-norm (\code{"2"}), as determined by the value of
  \code{type}.
}
\arguments{
  \item{x}{
    a real or complex matrix.
  }
  \item{type}{
    A character indicating the type of norm desired.
    \describe{
      \item{\code{"O"}, \code{"o"} or \code{"1"}}{specifies the one norm,
	(maximum absolute column sum);}
      \item{\code{"I"} or \code{"i"}}{specifies the infinity norm (maximum
	absolute row sum);}
      \item{\code{"F"} or \code{"f"}}{specifies the Frobenius norm (the
	Euclidean norm of \code{x} treated as if it were a vector);}
      \item{\code{"M"} or \code{"m"}}{specifies the maximum modulus of
	all the elements in \code{x}; and}
      \item{\code{"2"}}{specifies the \dQuote{spectral norm} or {2-norm}, which
	is the largest singular value (\code{\link{svd}}) of \code{x}.}
    }
    The default is \code{"O"}.  Only the first character of
    \code{type[1]} is used.
  }
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
    A numeric value of class \code{"norm"}, representing the quantity
    chosen according to \code{type}.
}
\details{
  For dense matrices, the methods eventually call the Lapack functions
  \code{dlange}, \code{dlansy}, \code{dlantr}, \code{zlange},
  \code{zlansy}, and \code{zlantr}.
}
\seealso{
  \code{\link{onenormest}()}, an \emph{approximate} randomized estimate
  of the 1-norm condition number, efficient for large sparse matrices.

  The \code{\link[base]{norm}()} function from \R's \pkg{base} package.
}
\references{
  Anderson, E., et al. (1994).
  \emph{LAPACK User's Guide,}
  2nd edition, SIAM, Philadelphia.
}
\examples{
x <- Hilbert(9)
norm(x)# = "O" = "1"
stopifnot(identical(norm(x), norm(x, "1")))
norm(x, "I")# the same, because 'x' is symmetric

allnorms <- function(d) vapply(c("1","I","F","M","2"),
                               norm, x = d, double(1))
allnorms(x)
allnorms(Hilbert(10))

i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
A <- sparseMatrix(i, j, x = x)                      ##  8 x 10 "dgCMatrix"
(sA <- sparseMatrix(i, j, x = x, symmetric = TRUE)) ## 10 x 10 "dsCMatrix"
(tA <- sparseMatrix(i, j, x = x, triangular= TRUE)) ## 10 x 10 "dtCMatrix"
(allnorms(A) -> nA)
allnorms(sA)
allnorms(tA)
stopifnot(all.equal(nA, allnorms(as.matrix(A))),
	  all.equal(nA, allnorms(tA))) # because tA == rbind(A, 0, 0)
A. <- A; A.[1,3] <- NA
stopifnot(is.na(allnorms(A.))) # gave error
}
\keyword{algebra}
