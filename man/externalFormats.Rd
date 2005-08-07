\name{externalFormats}
\alias{readHB}
\alias{readMM}
\title{Read and write external matrix formats}
\description{
    Read matrices stored in the Harwell-Boeing or MatrixMarket formats.
}
\usage{
readHB(filename)
readMM(filename)
}
\arguments{
  \item{filename}{the name of the file to read, as a character scalar.
    The names of files storing matrices in the Harwell-Boeing format
    usually end in \code{".rua"} or \code{".rsa"}.  Those storing
    matrices in the MatrixMarket format usually end in {".mtx"}.}
}
\value{
  An object that inherits from the \code{"Matrix"} class.
}
\references{
  \url{http://math.nist.gov/MatrixMarket}

  \url{http://www.cise.ufl.edu/research/sparse/matrices}
}
\examples{
str(pores <- readMM(system.file("external/pores_1.mtx",
                                package = "Matrix")))
str(utm <- readHB(system.file("external/utm300.rua",
                               package = "Matrix")))
str(lundA <- readMM(system.file("external/lund_a.mtx",
                                package = "Matrix")))
str(lundA <- readHB(system.file("external/lund_a.rsa",
                                package = "Matrix")))
}
\keyword{IO}
\keyword{algebra}

