\name{ScotsSec}
\alias{ScotsSec}
\docType{data}
\title{Scottish secondary school scores}
\description{
  Scores attained by 3435 Scottish secondary school students on a
  standardized test taken at age 16.  Both the primary school and the
  secondary school that the student attended have been recorded.
}
\usage{data(ScotsSec)}
\format{
  A data frame with 3435 observations on the following 6 variables.
  \describe{
    \item{verbal}{The verbal reasoning score on a test taken by the
      students on entry to secondary school.}
    \item{attain}{The score attained on the standardized test taken at
      age 16.}
    \item{primary}{A factor indicating the primary school that the
      student attended.}
    \item{sex}{A factor with levels \code{M} and \code{F}}
    \item{social}{The student's social class on a numeric scale from low
      to high social class.}
    \item{second}{A factor indicating the secondary school that the
      student attended.}
  }
}
\details{
  These data are an example of cross-classified grouping factors.
}
\source{
  \url{http://multilevel.ioe.ac.uk/softrev/2lev-xc.html}
}
\references{
  Paterson, L. (1991). Socio economic status and educational attainment:
  a multidimensional and multilevel study. \emph{Evaluation and Research
    in Education} \bold{5}: 97-121.
}
\examples{
data(ScotsSec)
str(ScotsSec)
}
\keyword{datasets}
