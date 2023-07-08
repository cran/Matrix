\name{coerce-methods-graph}
\title{Conversions "graph" <--> (sparse) Matrix}
%
\docType{methods}
\keyword{methods}
\keyword{utilities}
%
\alias{coerce,Matrix,graph-method}
\alias{coerce,Matrix,graphNEL-method}
\alias{coerce,TsparseMatrix,graphNEL-method}
\alias{coerce,graph,CsparseMatrix-method}
\alias{coerce,graph,Matrix-method}
\alias{coerce,graph,RsparseMatrix-method}
\alias{coerce,graph,TsparseMatrix-method}
\alias{coerce,graph,sparseMatrix-method}
\alias{coerce,graphAM,TsparseMatrix-method}
\alias{coerce,graphNEL,TsparseMatrix-method}
%
\alias{T2graph}
\alias{graph2T}
%
\description{
  Since 2005, package \pkg{Matrix} has supported coercions to and
  from class \code{\link[graph:graph-class]{graph}} from package
  \href{https://bioconductor.org/packages/graph/}{\pkg{graph}}.
  Since 2013, this functionality has been exposed via functions
  \code{T2graph} and \code{graph2T}, which, unlike methods for
  \code{\link{as}(from, "<Class>")}, support optional arguments.
}
\usage{
graph2T(from, use.weights = )
T2graph(from, need.uniq = is_not_uniqT(from), edgemode = NULL)
}
\arguments{
  \item{from}{for \code{graph2T()}, an \R object of class
    \code{"graph"}; \cr
    for \code{T2graph()}, a sparse matrix inheriting from
    \code{"\linkS4class{TsparseMatrix}"}.}
  \item{use.weights}{logical indicating if weights should be used, i.e.,
    equivalently the result will be numeric, i.e. of class
    \code{\linkS4class{dgTMatrix}}; otherwise the result will be
    \code{\linkS4class{ngTMatrix}} or \code{\linkS4class{nsTMatrix}},
    the latter if the graph is undirected.  The default looks if there
    are weights in the graph, and if any differ from \code{1}, weights
    are used.}
  \item{need.uniq}{a logical indicating if \code{from} may need to be
    internally \dQuote{uniqified}; do not set this and hence rather use
    the default, unless you know what you are doing!}
  \item{edgemode}{one of \code{NULL}, \code{"directed"}, or
    \code{"undirected"}.  The default \code{NULL} looks if the matrix is
    symmetric and assumes \code{"undirected"} in that case.}
}
\value{
  For \code{graph2T()}, a sparse matrix  inheriting from
  \code{"\linkS4class{TsparseMatrix}"}.

  For \code{T2graph()} an \R object of class \code{"graph"}.
}
\seealso{
  Package \CRANpkg{igraph}, which provides similar coercions
  to and from its class \code{igraph} via functions
  \code{graph_from_adjacency_matrix} and \code{as_adjacency_matrix}.
}
\examples{
if(requireNamespace("graph")) {
  n4 <- LETTERS[1:4]; dns <- list(n4,n4)
  show(a1 <- sparseMatrix(i= c(1:4),   j=c(2:4,1),   x = 2,    dimnames=dns))
  show(g1 <- as(a1, "graph")) # directed
  unlist(graph::edgeWeights(g1)) # all '2'

  show(a2 <- sparseMatrix(i= c(1:4,4), j=c(2:4,1:2), x = TRUE, dimnames=dns))
  show(g2 <- as(a2, "graph")) # directed
  # now if you want it undirected:
  show(g3  <- T2graph(as(a2,"TsparseMatrix"), edgemode="undirected"))
  show(m3 <- as(g3,"Matrix"))
  show( graph2T(g3) ) # a "pattern Matrix" (nsTMatrix)
\dontshow{
  stopifnot(
   identical(as(g3,"Matrix"), as(as(a2 + t(a2), "nMatrix"),"symmetricMatrix"))
  ,
   identical(tg3 <- graph2T(g3), graph2T(g3, use.weights=FALSE))
  ,
   identical(as(m3,"TsparseMatrix"), uniqTsparse(tg3))
  )
}
  a. <- sparseMatrix(i=4:1, j=1:4, dimnames=list(n4, n4), repr="T") # no 'x'
  show(a.) # "ngTMatrix"
  show(g. <- as(a., "graph"))
\dontshow{
  stopifnot(graph::edgemode(g.) == "undirected",
            graph::numEdges(g.) == 2,
            all.equal(as(g., "TsparseMatrix"),
                      as(a., "symmetricMatrix"))
)
}
}
}
