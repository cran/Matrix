\name{writeHB}
\docType{genericFunction}
\alias{writeHB}
\alias{writeMM}
\alias{writeHB,dgCMatrix-method}
\alias{writeHB,dgTMatrix-method}
\alias{writeHB,dsCMatrix-method}
\alias{writeHB,dsTMatrix-method}
\alias{writeMM,dgCMatrix-method}
\alias{writeMM,dgTMatrix-method}
\alias{writeMM,dsCMatrix-method}
\alias{writeMM,dsTMatrix-method}
\title{Write Harwell-Boeing or MatrixMarket formats}
\description{
  Write a matrix to a file in the Harwell-Boeing format or in the
  MatrixMarket format.
}
\usage{
writeHB(obj, file, ...)
writeMM(obj, file, ...)
}
\arguments{
 \item{obj}{a real sparse matrix}
 \item{file}{character - the name of the file to be written.}
 \item{\dots}{optional additional arguments. Currently none are used in
   any methods.}
}
\value{
  NULL.  As a side effect, the matrix \code{obj} is written to
  \code{file} in the Harwell-Boeing format (writeHB)
  or the MatrixMarket format (writeMM).
}
\examples{
writeHB(mm, "mmHB.rua")
writeMM(mm, "mmMM.mtx")
}
\keyword{IO}
\keyword{array}
\keyword{algebra}
