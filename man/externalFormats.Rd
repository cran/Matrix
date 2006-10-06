\name{externalFormats}
\docType{genericFunction}
\alias{readHB}
\alias{readMM}
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
\title{Read and write external matrix formats}
\description{
    Read matrices stored in the Harwell-Boeing or MatrixMarket formats
    or write sparseMatrix objects to one of these formats.
}
\usage{
readHB(file)
readMM(file)
writeHB(obj, file, ...)  # deprecated
writeMM(obj, file, ...)
}
\arguments{
 \item{obj}{a real sparse matrix}
 \item{file}{for \code{writeMM} - the name of the file to be written.
   For \code{readHB} and \code{readMM} the name of the file to read, as
   a character scalar.  The names of files storing matrices in the
   Harwell-Boeing format usually end in \code{".rua"} or \code{".rsa"}.
   Those storing matrices in the MatrixMarket format usually end in
   {".mtx"}.

   Alternatively, \code{readHB} and \code{readMM} accept connection objects.}
 \item{\dots}{optional additional arguments. Currently none are used in
   any methods.}
}
\value{
  The \code{readHB} and \code{readMM} functions return an object that
  inherits from the \code{"Matrix"} class.  Methods for the
  \code{writeMM} generic functions usually return
  \code{NULL} and, as a side effect, the matrix \code{obj} is written to
  \code{file} in the MatrixMarket format (writeMM).
}
\note{
  The Harwell-Boeing format is older and less flexible than the
 MatrixMarket format.  The function \code{writeHB} is now deprecated.
 We recommend the use of \code{writeMM} instead.
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
\dontrun{
## NOTE: The following examples take quite some time
## ----  even on a fast internet connection:
str(sm <-
 readHB(gzcon(url("http://www.cise.ufl.edu/research/sparse/HBformat/Boeing/msc00726.rsa.gz"))))
str(jgl009 <-
 readMM(gzcon(url("ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/counterx/jgl009.mtx.gz"))))
}
data(KNex)
\dontrun{
writeHB(KNex$mm, "mmHB.rua")
}
writeMM(KNex$mm, "mmMM.mtx")
}
\keyword{IO}
\keyword{array}
\keyword{algebra}

