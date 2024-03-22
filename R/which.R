## METHODS FOR GENERIC: which
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("which", c(x = "ndenseMatrix"),
          function(x, arr.ind = FALSE, useNames = TRUE) {
              wh <- which(.M2v(x))
              if(arr.ind)
                  arrayInd(wh, x@Dim, dimnames(x), useNames = useNames)
              else wh
          })
setMethod("which", c(x = "ldenseMatrix"),
          function(x, arr.ind = FALSE, useNames = TRUE) {
              wh <- which(.M2v(x))
              if(arr.ind)
                  arrayInd(wh, x@Dim, dimnames(x), useNames = useNames)
              else wh
          })

setMethod("which", c(x = "nsparseMatrix"),
          function(x, arr.ind = FALSE, useNames = TRUE) {
              wh <- .M2V(x)@i
              if(arr.ind)
                  arrayInd(wh, x@Dim, dimnames(x), useNames = useNames)
              else wh
          })
setMethod("which", c(x = "lsparseMatrix"),
          function(x, arr.ind = FALSE, useNames = TRUE) {
              wh <- { x. <- .M2V(x); x.@i[which(x.@x)] }
              if(arr.ind)
                  arrayInd(wh, x@Dim, dimnames(x), useNames = useNames)
              else wh
           })

setMethod("which", c(x = "ndiMatrix"),
          function(x, arr.ind = FALSE, useNames = TRUE) {
              wh <- .M2V(x)@i
              if(arr.ind)
                  arrayInd(wh, x@Dim, x@Dimnames, useNames = useNames)
              else wh
          })
setMethod("which", c(x = "ldiMatrix"),
          function(x, arr.ind = FALSE, useNames = TRUE) {
              wh <- { x. <- .M2V(x); x.@i[which(x.@x)] }
              if(arr.ind)
                  arrayInd(wh, x@Dim, x@Dimnames, useNames = useNames)
              else wh
          })

setMethod("which", c(x = "nsparseVector"),
          function(x, arr.ind = FALSE, useNames = TRUE) x@i)
setMethod("which", c(x = "lsparseVector"),
          function(x, arr.ind = FALSE, useNames = TRUE) x@i[which(x@x)])

setMethod("which", c(x = "indMatrix"),
          function(x, arr.ind = FALSE, useNames = TRUE) {
              wh <- .M2V(x)@i
              if(arr.ind)
                  arrayInd(wh, x@Dim, x@Dimnames, useNames = useNames)
              else wh
          })
