## METHODS ENHANCING PACKAGE: graph
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## README: These need package 'graph', from which we do not import,
##         _on purpose_.  Hence we must use graph:: in case 'graph'
##         is only loaded, not attached ...

## NB: undirected graph <==> symmetric matrix


## ~~~~ UTILITIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## NB: These may no longer be needed in future versions of 'graph'

graph.has.weights <- function(g)
    "weight" %in% names(graph::edgeDataDefaults(g))

graph.non.1.weights <- function(g)
    any(unlist(graph::edgeData(g, attr = "weight")) != 1)

graph.wgtMatrix <- function(g)
{
    ## Purpose: work around "graph" package's  as(g, "matrix") bug
    ## ----------------------------------------------------------------------
    ## Arguments: g: an object inheriting from (S4) class "graph"
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, based on Seth Falcon's code;  Date: 12 May 2006

    ## MM: another buglet for the case of  "no edges":
    if(graph::numEdges(g) == 0) {
        p <- length(nd <- graph::nodes(g))
        return( matrix(0, p,p, dimnames = list(nd, nd)) )
    }

    ## Usual case, when there are edges:
    if(has.w <- graph.has.weights(g)) {
        ## graph.non.1.weights(g) :
        w <- unlist(graph::edgeData(g, attr = "weight"))
        has.w <- any(w != 1)
    }
    ## now 'has.w' is TRUE  iff  there are weights != 1
    m <- as(g, "matrix")
    ## now is a 0/1 - matrix (instead of 0/wgts) with the 'graph' bug
    if(has.w) { ## fix it if needed
        tm <- t(m)
        tm[tm != 0] <- w
        t(tm)
    } else m
}

graph2T <- function(from, use.weights =
                              graph.has.weights(from) &&
                              graph.non.1.weights(from))
{
    nd <- graph::nodes(from); dnms <- list(nd,nd)
    dm <- rep.int(length(nd), 2)
    edge2i <- function(e) {
	## return (0-based) row indices 'i'
	rep.int(0:(dm[1]-1L), lengths(e))
    }

    if(use.weights) {
	eWts <- graph::edgeWeights(from); names(eWts) <- NULL
	i <- edge2i(eWts)
	To <- unlist(lapply(eWts, names))
	j <- as.integer(match(To,nd)) - 1L # columns indices (0-based)
	## symm <- symm && <weights must also be symmetric>: improbable
	## if(symm) new("dsTMatrix", .....) else
	new("dgTMatrix", i = i, j = j, x = unlist(eWts),
            Dim = dm, Dimnames = dnms)
    } else { ## no weights: 0/1 matrix -> pattern
	edges <- lapply(from@edgeL[nd], "[[", "edges")
	symm <- graph::edgemode(from) == "undirected"
	if(symm)# each edge appears twice; keep upper triangle only
	    edges <- lapply(seq_along(edges),
                            function(i) {e <- edges[[i]]; e[e >= i]})
	i <- edge2i(edges)
	j <- as.integer(unlist(edges)) - 1L # column indices (0-based)
	## if(symm) {			# symmetric: ensure upper triangle
	##     tmp <- i
	##     flip <- i > j
	##     i[flip] <- j[flip]
	##     j[flip] <- tmp[flip]
	##     new("nsTMatrix", i = i, j = j, Dim = dm, Dimnames = dnms, uplo = "U")
	## } else {
	##     new("ngTMatrix", i = i, j = j, Dim = dm, Dimnames = dnms)
	## }
	new(if(symm) "nsTMatrix" else "ngTMatrix", i = i, j = j,
	    Dim = dm, Dimnames = dnms)# uplo = "U" is default
    }
}

T2graph <- function(from, need.uniq = is_not_uniqT(from), edgemode = NULL)
{
    d <- dim(from)
    if((n <- d[1L]) != d[2L])
	stop("only square matrices can be used as graph incidence matrices")
    if(n == 0L)
        return(new("graphNEL"))
    if(is.null(rn <- dimnames(from)[[1]]))
	rn <- as.character(1:n)
    if(need.uniq) ## Need to 'uniquify' the triplets!
	from <- uniqTsparse(from)

    if(is.null(edgemode))
        edgemode <-
            if(isSymmetric(from)) { # either "symmetricMatrix" or otherwise
                ##-> undirected graph: every edge only once!
                if(!is(from, "symmetricMatrix")) {
                    ## a general matrix which happens to be symmetric
                    ## ==> remove the double indices
                    from <- tril(from)
                }
                "undirected"
            } else {
                "directed"
            }
    ## every edge is there only once, either upper or lower triangle
    ft1 <- cbind(rn[from@i + 1L], rn[from@j + 1L])
    graph::ftM2graphNEL(ft1,
                        W = if(.hasSlot(from,"x")) as.numeric(from@x),
			V = rn,
                        edgemode = edgemode)
}


## ~~~~ COERCIONS FROM graph TO Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("graphAM", "TsparseMatrix",
      function(from) {
          symm <- graph::edgemode(from) == "undirected" &&
              isSymmetric(from@adjMat)
	  if(graph.has.weights(from))
	      .m2sparse(graph.wgtMatrix(from), if(symm) "dsT" else "dgT")
          else ## no weights: 0/1 matrix -> pattern
	      .m2sparse(   as(from, "matrix"), if(symm) "nsT" else "ngT")
      })
setAs("graphNEL", "TsparseMatrix",
      function(from) graph2T(from))

setAs("graph", "CsparseMatrix",
      function(from) .T2C(as(from, "TsparseMatrix")))
setAs("graph", "RsparseMatrix",
      function(from) .T2R(as(from, "TsparseMatrix")))
setAs("graph", "TsparseMatrix",
      function(from) graph2T(as(from, "graphNEL")))
setAs("graph", "sparseMatrix",
      function(from) as(from, "CsparseMatrix"))
setAs("graph", "Matrix",
      function(from) as(from, "CsparseMatrix"))


## ~~~~ COERCIONS FROM Matrix TO graph ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## NB: Since we have a method for TsparseMatrix (below), the method
##     for Matrix can assume that 'from' is _not_ a TsparseMatrix,
##     and therefore expect that as(from, "TsparseMatrix") is "unique"
setAs("TsparseMatrix", "graphNEL",
      function(from) T2graph(from))
setAs("Matrix", "graphNEL",
      function(from) T2graph(as(from, "TsparseMatrix"), need.uniq = FALSE))
setAs("Matrix", "graph",
      function(from) as(from, "graphNEL"))
