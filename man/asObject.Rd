\name{asObject}
\alias{asObject}
\alias{prependClass}
\title{Assign or Modify the Class of an Object}
\description{
    \code{asObject} assigns a class to an object and returns the object.
    \code{prependClass} prepends a class to the set of classes for an
    object and returns the object.
}
\usage{
asObject(x, cl)
prependClass(x, cl)
}
\arguments{
  \item{x}{an object.}
  \item{cl}{a character vector of classes.}
}
\details{
    \code{asObject} is simply a wrapper for the sequence
    \code{class(x) <- cl; x}.  Its main use is in constructing the
    value that a function will return.   It differs from the assignment
    \code{class(x) <- cl} in that the assignment returns the class
    whereas \code{asObject} returns the object with the class.

    \code{prependClass} is a similar function that extends the classes
    (if any) of the object and returns the object.
}
\value{
    the object \code{x} with the class \code{cl} (or \code{c(cl,
	class(x))} for \code{prependClass}).
}
\examples{
x <- asObject(3, "myClass")
print.default(x)
print.default(prependClass(x, "firstClass"))
}
\keyword{methods}
\keyword{classes}
