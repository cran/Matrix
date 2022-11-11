## METHODS FOR GENERIC: which
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("which", "ndenseMatrix",
          function(x, arr.ind, useNames) {
              wh <- which(isN0(.dense2g(x)@x)) # NA <=> TRUE
              if(arr.ind)
                  arrayInd(wh, x@Dim, dimnames(x), useNames = useNames)
              else wh
          })
setMethod("which", "ldenseMatrix",
          function(x, arr.ind, useNames) {
              wh <- which(.dense2g(x)@x)
              if(arr.ind)
                  arrayInd(wh, x@Dim, dimnames(x), useNames = useNames)
              else wh
          })

setMethod("which", "nsparseMatrix",
	  function(x, arr.ind, useNames = TRUE) {
              if(arr.ind)
                  which(.CR2T(x), arr.ind = TRUE, useNames = useNames)
	      else as(x, "sparseVector")@i
	  })
setMethod("which", "lsparseMatrix", # to be called only for l.[CR]Matrix
	  function(x, arr.ind, useNames = TRUE) {
              if(arr.ind)
                  which(.CR2T(x), arr.ind = TRUE, useNames = useNames)
	      else which(as(x, "sparseVector"))
	  })

setMethod("which", "ldiMatrix",
	  function(x, arr.ind, useNames) {
              d <- x@Dim
              i <- indDiag(n = d[1L], packed = FALSE)
              if(x@diag == "N")
                  i <- i[which(x@x)]
              if(arr.ind)
                  arrayInd(i, d, x@Dimnames, useNames = useNames)
	      else i
          })

## Constructs 'dimnames' as arrayInd(*, useNames = TRUE):
.arr.DN <- function(ind, .dimnames)
    list(.dimnames[[1L]][ind],
         if(any(nzchar(nd <- names(.dimnames)))) nd else c("row", "col"))

.which.ngT <- function(x, arr.ind, useNames = TRUE)
    if(arr.ind) {
        ij <- cbind(x@i, x@j) + 1L
        if(useNames)
            dimnames(ij) <- .arr.DN(ij[, 1L], x@Dimnames)
        ij
    } else as(x, "sparseVector")@i

.which.lgT <- function(x, arr.ind, useNames = TRUE) {
    if(arr.ind) {
	k <- which(x@x)
	ij <- cbind(x@i[k], x@j[k]) + 1L
        if (useNames)
            dimnames(ij) <- .arr.DN(ij[, 1L], x@Dimnames)
        ij
    } else which(as(x, "sparseVector"))
}

setMethod("which", "ngTMatrix", .which.ngT)
setMethod("which", "lgTMatrix", .which.lgT)

setMethod("which", "ntTMatrix",
          function(x, arr.ind, useNames = TRUE)
              .which.ngT(.Call(R_sparse_diag_U2N, x), arr.ind, useNames))
setMethod("which", "ltTMatrix",
          function(x, arr.ind, useNames = TRUE)
              .which.lgT(.Call(R_sparse_diag_U2N, x), arr.ind, useNames))

setMethod("which", "nsTMatrix",
          function(x, arr.ind, useNames = TRUE)
	      .which.ngT(.sparse2g(x), arr.ind, useNames))
setMethod("which", "lsTMatrix",
          function(x, arr.ind, useNames = TRUE)
              .which.lgT(.sparse2g(x), arr.ind, useNames))

setMethod("which", "nsparseVector",
          function(x, arr.ind, useNames) x@i)
setMethod("which", "lsparseVector",
          function(x, arr.ind, useNames) x@i[isT(x@x)])
