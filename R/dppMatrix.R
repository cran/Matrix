## METHODS FOR CLASS: dppMatrix
## dense (packed) symmetric positive definite matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.dsp2dpp <- function(from) {
    if(is.null(tryCatch(.Call(dppMatrix_trf, from, 2L),
                        error = function(e) NULL)))
        stop("not a positive definite matrix")
    ## FIXME: check=FALSE
    copyClass(from, "dppMatrix",
              sNames = c("Dim", "Dimnames", "uplo", "x", "factors"))
}

setAs("dspMatrix", "dppMatrix", .dsp2dpp)

setAs("dsyMatrix", "dppMatrix",
      function(from) pack(.dsy2dpo(from)))

setAs("matrix", "dppMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .dsp2dpp(pack(from, symmetric = TRUE))
      })

setAs("Matrix", "dppMatrix",
      function(from) {
          ## still needs as(<dspMatrix>, "dppMatrix") to work
          as(as(as(as(from,"dMatrix"),"symmetricMatrix"),"packedMatrix"),
             "dppMatrix")
      })


## ~~~~ COERCIONS FROM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("dppMatrix", "dpoMatrix", function(from) unpack(from))

## MJ: no longer needed
if(FALSE) {
setAs("dppMatrix", "dpoMatrix",
      function(from) {
          ## FIXME: check=FALSE
          copyClass(.Call(dspMatrix_as_dsyMatrix, from),
                    "dpoMatrix",
                    sNames = c("x", "Dim", "Dimnames", "uplo", "factors"))
      })
} ## MJ

## MJ: redundant, as coercions are inherited from superclass dspMatrix
if(FALSE) {
dpp2sC <- function(from) as(.Call(dspMatrix_as_dsyMatrix, from), "dsCMatrix")
## setAs("dppMatrix", "dsCMatrix", dpp2sC)
setAs("dppMatrix", "CsparseMatrix", dpp2sC)
setAs("dppMatrix", "sparseMatrix", dpp2sC)
} ## MJ

## MJ: no longer needed ... replacement in ./denseMatrix.R
## (was infelicitous anyway because result did not have packed storage)
if(FALSE) {
setAs("dppMatrix", "lMatrix",
      function(from) as(as(from, "dsyMatrix"), "lMatrix"))
setAs("dppMatrix", "nMatrix",
      function(from) as(as(from, "dsyMatrix"), "nMatrix"))
} ## MJ


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: no longer needed ... replacement in ./packedMatrix.R
if(FALSE) {
setMethod("t", signature(x = "dppMatrix"),
          function(x) as(t(as(x, "dspMatrix")), "dppMatrix"),
          valueClass = "dppMatrix")
} ## MJ
