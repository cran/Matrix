## METHODS FOR CLASS: denseMatrix (virtual)
## dense matrices with unpacked _or_ packed storage
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dim<-", signature(x = "denseMatrix"),
          function(x, value) {
              if(!is.numeric(value) || length(value) != 2L)
                  stop("dimensions must be numeric of length 2")
              if(anyNA(value))
                  stop("dimensions cannot contain NA")
              if(any(value < 0))
                  stop("dimensions cannot contain negative values")
              if(!is.integer(value)) {
                  if(any(value > .Machine$integer.max))
                      stop("dimensions cannot exceed 2^31-1")
                  value <- as.integer(value)
              }
              if(all(value == (d <- x@Dim)))
                  return(x)
              if((pv <- prod(value)) != (pd <- prod(d)))
                  stop(gettextf("assigned dimensions [product %.0f] do not match Matrix length [%.0f]",
                                pv, pd, domain = NA))
              r <- .M2gen(x)
              r@Dim <- value
              r@factors <- list()
              r
          })

setMethod("mean", signature(x = "denseMatrix"),
          function(x, trim = 0, na.rm = FALSE, ...) {
              if(is.numeric(trim) && length(trim) == 1L && !is.na(trim) &&
                 trim == 0) {
                  ## Be fast in this special case :
                  if(isTRUE(na.rm))
                      x <- x[!is.na(x)]
                  sum(x) / length(x)
              } else mean.default(.M2v(x), trim = trim, na.rm = na.rm, ...)
          })

setMethod("rep", "denseMatrix",
          function(x, ...) rep(.M2v(x), ...))

setMethod("show", "denseMatrix",
          function(object) prMatrix(object))

.dense.band <- function(x, k1, k2, ...) .Call(R_dense_band, x, k1, k2)
.dense.triu <- function(x, k = 0L, ...) .Call(R_dense_band, x, k, NULL)
.dense.tril <- function(x, k = 0L, ...) .Call(R_dense_band, x, NULL, k)
for (.cl in c("denseMatrix", "matrix")) {
    setMethod("band", signature(x = .cl), .dense.band)
    setMethod("triu", signature(x = .cl), .dense.triu)
    setMethod("tril", signature(x = .cl), .dense.tril)
}
rm(.dense.band, .dense.triu, .dense.tril, .cl)

## x[] <- value :
setReplaceMethod("[", signature(x = "denseMatrix", i = "missing", j = "missing",
                                value = "ANY"),## double/logical/...
                 function (x, value) {
                     x <- .M2gen(x)
                     x@x[] <- value
                     validObject(x)# check if type and lengths above match
                     x
                 })

## FIXME: 1) These are far from efficient
## -----
setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "missing",
                                value = "replValue"),
                 function (x, i, j, ..., value) {
                     r <- as(x, "matrix")
                     ## message("`[<-` with nargs()= ",nargs())
                     if((na <- nargs()) == 3)
                         r[i] <- value
                     else if(na == 4)
                         r[i, ] <- value
                     else stop(gettextf("invalid nargs()= %d", na), domain=NA)
                     .m2dense(r, paste0(.M.kind(x), "ge"))
                 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "missing", j = "index",
                                value = "replValue"),
                 function (x, i, j, ..., value) {
                     r <- as(x, "matrix")
                     r[, j] <- value
                     .m2dense(r, paste0(.M.kind(x), "ge"))
                 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "index",
                                value = "replValue"),
                 function (x, i, j, ..., value) {
                     r <- as(x, "matrix")
                     r[i, j] <- value
                     as_denseClass(r, class(x)) ## was as(r, class(x))
                 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "matrix",  # 2-col.matrix
                                j = "missing", value = "replValue"),
                 function(x, i, j, ..., value) {
                     r <- as(x, "matrix")
                     r[ i ] <- value
                     .m2dense(r, paste0(.M.kind(x), "ge"))
                 })
