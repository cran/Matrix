## "gTMatrix" : Virtual class of general sparse matrices in triplet-format

## Select rows
setMethod("[", signature(x = "gTMatrix", i = "numeric", j = "missing",
			 drop = "logical"),
	  function (x, i, j, ..., drop) { ## select rows
	      storage.mode(i) <- "integer"
              xi <- x@i + 1:1 # 1-indexing
	      sel <- xi %in% i
              nd <- c(length(i), ncol(x))
	      x@Dim <- nd
              x@i <- match(xi[sel], i) - 1:1
              x@j <- x@j[sel]
              x@x <- x@x[sel]
	      if (drop && any(nd == 1)) drop(as(x,"matrix")) else x
	  })


## Select columns
setMethod("[", signature(x = "gTMatrix", i = "missing", j = "numeric",
			 drop = "logical"),
	  function (x, i, j, ..., drop) { ## select columns
	      storage.mode(j) <- "integer"
              xj <- x@j + 1:1 # 1-indexing
	      sel <- xj %in% j
              nd <- c(nrow(x), length(j))
	      x@Dim <- nd
              x@i <-  x@i[sel]
              x@j <- match(xj[sel], j) - 1:1
              x@x <- x@x[sel]
	      if (drop && any(nd == 1)) drop(as(x,"matrix")) else x
	  })


## [.data.frame has : drop = if (missing(i)) TRUE else length(cols) == 1)

setMethod("[", signature(x = "gTMatrix",
			 i = "numeric", j = "numeric", drop = "logical"),
	  function (x, i, j, ..., drop)
      {
	  ## (i,j, drop) all specified
	  storage.mode(i) <- "integer"
	  storage.mode(j) <- "integer"
          xi <- x@i + 1:1
          xj <- x@j + 1:1
	  sel <- (xi %in% i) & (xj %in% j)
          nd <- c(length(i), length(j))
          x@Dim <- nd
          x@i <- match(xi[sel], i) - 1:1
          x@j <- match(xj[sel], j) - 1:1
          x@x <- x@x[sel]
          if (drop && any(nd == 1)) drop(as(x,"matrix")) else x
      })
