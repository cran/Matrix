\name{la.env}
\alias{la.env}
\title{
Set Tuning Parameters for Linear Algebra Computations
}
\description{
Sets values of tuning parameters within LAPACK.
}
\usage{
la.env(NB=.laenv$NB, NBMIN=.laenv$NBMIN, NX=.laenv$NX, NS=.laenv$NS, 
       NXSVD=.laenv$NXSVD, MAXB=.laenv$MAXB)
}
\arguments{
\item{NB}{
optimal blocksize. 
}
\item{NBMIN}{
minimum block size.
If the usable block size is less than \code{NBMIN}, an unblocked method is used.
The value of \code{NBMIN} is irrelevant if \code{NB = 1}.
}
\item{NX}{
crossover point.
In a block method, when blocks reach dimension less than \code{NX}, 
an unblocked method is used. 
The value of \code{NX} is irrelevant if \code{NB = 1}.
}
\item{NS}{
number of shifts for nonsymmetric eigenvalue computations.
}
\item{NXSVD}{
crossover point for the singular-value computation. When reducing an 
m by n matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds this value,
the matrix is first reduced to triangular form by QR factorization.
}
\item{MAXB}{
crossover point within nonsymmetric eigenvalue computations.
}}
\value{
Returns the list of parameters use din the call to \code{la.env}, with missing
parameters replaced by their default values.
}
\section{SIDE EFFECTS}{
The corresponding parameters within Lapack are set to the values supplied 
in \code{la.env}.
}
\section{DETAILS}{
The following is from Anderson et al. (1994) :


The range of problem sizes needed to determine the optimal block size or
crossover point is machine dependent. For algorithms that require a
crossover point, it is best to start by finding the best block size with the
crossover point set to 0, and then locate the point at which blocked algorithm
begins to outperform then unblocked algorithm using this block size. The best 
crossover point will be somewhat smaller than the latter value.


By experimenting with small values, it should be straightforward to choose 
\code{NBMIN}, the smallest block size that gives improvement over an unblocked
method.
}
\source{
Anderson, E., et al. (1994),
\emph{LAPACK Users' Guide.}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{.laenv}}
}
\examples{
la.env( NB = 4, NBMIN = 2, NX = 2) 
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.










