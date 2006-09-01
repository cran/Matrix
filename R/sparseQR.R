## Methods for the sparse QR decomposition

## The signature should change to y = "ddenseMatrix" later
setMethod("qr.qy", signature(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y) .Call(sparseQR_qty, qr, y, FALSE),
          valueClass = "dgeMatrix")

setMethod("qr.qy", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y) .Call(sparseQR_qty, qr, y, FALSE),
          valueClass = "dgeMatrix")

setMethod("qr.qy", signature(qr = "sparseQR", y = "numeric"),
          function(qr, y) .Call(sparseQR_qty, qr, y, FALSE),
          valueClass = "dgeMatrix")

## The signature should change to y = "ddenseMatrix" later
setMethod("qr.qty", signature(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y) .Call(sparseQR_qty, qr, y, TRUE),
          valueClass = "dgeMatrix")

setMethod("qr.qty", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y) .Call(sparseQR_qty, qr, y, TRUE),
          valueClass = "dgeMatrix")

setMethod("qr.qty", signature(qr = "sparseQR", y = "numeric"),
          function(qr, y) .Call(sparseQR_qty, qr, y, TRUE),
          valueClass = "dgeMatrix")

.coef.trunc <- function(qr, res) res[1:ncol(qr@R),,drop=FALSE]

## The signature should change to y = "ddenseMatrix" later
setMethod("qr.coef", signature(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y)
          .coef.trunc(qr, .Call(sparseQR_coef, qr, y)),
          valueClass = "dgeMatrix")

setMethod("qr.coef", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y)
          .coef.trunc(qr, .Call(sparseQR_coef, qr, y)),
          valueClass = "dgeMatrix")

setMethod("qr.coef", signature(qr = "sparseQR", y = "numeric"),
          function(qr, y)
          .coef.trunc(qr, .Call(sparseQR_coef, qr, y)),
          valueClass = "dgeMatrix")

## The signature should change to y = "ddenseMatrix" later
setMethod("qr.resid", signature(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y)
          .Call(sparseQR_resid_fitted, qr, y, TRUE),
          valueClass = "dgeMatrix")

setMethod("qr.resid", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y)
          .Call(sparseQR_resid_fitted, qr, y, TRUE),
          valueClass = "dgeMatrix")

setMethod("qr.resid", signature(qr = "sparseQR", y = "numeric"),
          function(qr, y)
          .Call(sparseQR_resid_fitted, qr, y, TRUE),
          valueClass = "dgeMatrix")

## The signature should change to y = "ddenseMatrix" later
setMethod("qr.fitted", signature(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y, k)
          .Call(sparseQR_resid_fitted, qr, y, FALSE),
          valueClass = "dgeMatrix")

setMethod("qr.fitted", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y, k)
          .Call(sparseQR_resid_fitted, qr, y, FALSE),
          valueClass = "dgeMatrix")

setMethod("qr.fitted", signature(qr = "sparseQR", y = "numeric"),
          function(qr, y, k)
          .Call(sparseQR_resid_fitted, qr, y, FALSE),
          valueClass = "dgeMatrix")
