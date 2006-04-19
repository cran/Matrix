\name{pedigree}
\alias{pedigree}
\title{Pedigree Constructor}
\description{
  Construct an object of class \code{"pedigree"}, more
  conveniently than by \code{new("pedigree", ....)}.
}
\usage{
pedigree(sire, dam, label)
}
\arguments{
    \item{sire}{numeric vector (with some \code{NA} entries) of integer IDs,
      denoting a \emph{previous} entry in the pedigree corresponding to
      the current entry's \dQuote{father}.}
    \item{dam}{similarly as \code{sire} for the \dQuote{mother} of each entry.}
    \item{label}{a vector coercable to \code{"character"} of the same
      length as \code{sire} and \code{dam} giving a unique ID for the
      corresponding entry.}
}
\value{
  an object of formal class \code{"pedigree"}, see \code{\link{pedigree-class}}.
}
\seealso{\code{\link{pedigree-class}}}
\examples{
example("pedigree-class") ## 'p1' pedigree object `the hard way'

ped <- pedigree(sire = c(NA,NA,1, 1,4,5),
                dam  = c(NA,NA,2,NA,3,2), label= 1:6)
## note that 'label' is coerced to character automatically
ped
stopifnot(identical(ped, p1))
}
\keyword{misc}

