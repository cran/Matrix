\name{graph-sparseMatrix}
\title{Conversions "graph" <--> (sparse) Matrix}
\alias{graph2T}
\alias{T2graph}
% graph stuff
\alias{coerce,graph,Matrix-method}
\alias{coerce,graph,sparseMatrix-method}
\alias{coerce,graph,CsparseMatrix-method}
\alias{coerce,graphAM,sparseMatrix-method}
\alias{coerce,graphNEL,CsparseMatrix-method}
\alias{coerce,graphNEL,TsparseMatrix-method}
\alias{coerce,sparseMatrix,graph-method}
\alias{coerce,sparseMatrix,graphNEL-method}
\alias{coerce,TsparseMatrix,graphNEL-method}
\alias{coerce,dgTMatrix,graphNEL-method}
\description{
  The \pkg{Matrix} package has supported conversion from and to
  \code{"\link[graph:graph-class]{graph}"} objects from (Bioconductor) package
  \pkg{graph} since summer 2005, via the usual \code{as(., "<class>")}
  coercion, \preformatted{  as(from, Class)
}
  Since 2013, this functionality is further exposed as the
  \code{graph2T()} and \code{T2graph()} functions (with further
  arguments than just \code{from}), which convert graphs to and from
  the triplet form of sparse matrices (of class
  \code{"\linkS4class{TsparseMatrix}"}) .
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
  Note that the CRAN package \pkg{igraph} also provides conversions from
  and to sparse matrices (of package \pkg{Matrix}) via its
  \code{\link[igraph]{graph.adjacency}()} and
  \code{\link[igraph]{get.adjacency}()}.
}
\examples{
if(isTRUE(try(require(graph)))) { ## super careful .. for "checking reasons"
  n4 <- LETTERS[1:4]; dns <- list(n4,n4)
  show(a1 <- sparseMatrix(i= c(1:4),   j=c(2:4,1),   x = 2,    dimnames=dns))
  show(g1 <- as(a1, "graph")) # directed
  unlist(edgeWeights(g1)) # all '2'

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
   identical(as(m3,"TsparseMatrix"), Matrix:::uniqTsparse(tg3))
  )
}
  a. <- sparseMatrix(i= 4:1, j=1:4, dimnames=list(n4,n4), giveC=FALSE) # no 'x'
  show(a.) # "ngTMatrix"
  show(g. <- as(a., "graph"))
\dontshow{
  stopifnot(edgemode(g.) == "undirected", numEdges(g.) == 2,
            all.equal(as(g., "TsparseMatrix"),
                      as(a., "symmetricMatrix"))
)
}
}% only if( graph )
}
\keyword{graph}
\keyword{utilities}
