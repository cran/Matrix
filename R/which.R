## METHODS FOR GENERIC: which
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("which", "ndenseMatrix",
          function(x, arr.ind = FALSE, useNames = TRUE) {
              wh <- which(.dense2g(x, "l")@x) # NA <=> TRUE
              if(arr.ind)
                  arrayInd(wh, x@Dim, dimnames(x), useNames = useNames)
              else wh
          })
setMethod("which", "ldenseMatrix",
          function(x, arr.ind = FALSE, useNames = TRUE) {
              wh <- which(.dense2g(x, "l")@x)
              if(arr.ind)
                  arrayInd(wh, x@Dim, dimnames(x), useNames = useNames)
              else wh
          })

setMethod("which", "nsparseMatrix",
          function(x, arr.ind = FALSE, useNames = TRUE) {
              wh <- which(as(x, "sparseVector"))
              if(arr.ind)
                  arrayInd(wh, x@Dim, dimnames(x), useNames = useNames)
              else wh
          })
setMethod("which", "lsparseMatrix",
          function(x, arr.ind = FALSE, useNames = TRUE) {
              wh <- which(as(x, "sparseVector"))
              if(arr.ind)
                  arrayInd(wh, x@Dim, dimnames(x), useNames = useNames)
              else wh
           })
setMethod("which", "ldiMatrix",
          function(x, arr.ind = FALSE, useNames = TRUE) {
              d <- x@Dim
              wh <- indDiag(d[1L])
              if(x@diag == "N")
                  wh <- wh[which(x@x)]
              if(arr.ind)
                  arrayInd(wh, d, x@Dimnames, useNames = useNames)
              else wh
          })

setMethod("which", "nsparseVector",
          function(x, arr.ind = FALSE, useNames = TRUE) x@i)
setMethod("which", "lsparseVector",
          function(x, arr.ind = FALSE, useNames = TRUE) x@i[which(x@x)])
