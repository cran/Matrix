## METHODS FOR CLASS: dpoMatrix
## dense (unpacked) symmetric positive definite matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.dsy2dpo <- function(from) {
    if(is.null(tryCatch(.Call(dpoMatrix_trf, from, 2L),
                        error = function(e) NULL)))
        stop("not a positive definite matrix")
    ## FIXME: check=FALSE
    copyClass(from, "dpoMatrix",
              sNames = c("Dim", "Dimnames", "uplo", "x", "factors"))
}

setAs("dsyMatrix", "dpoMatrix", .dsy2dpo)

setAs("dspMatrix", "dpoMatrix",
      function(from) unpack(.dsp2dpp(from)))

setAs("matrix", "dpoMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .dsy2dpo(.M2sym(from))
      })

setAs("Matrix", "dpoMatrix",
      function(from) {
          ## still needs as(<dsyMatrix>, "dpoMatrix") to work
          as(as(as(as(from,"dMatrix"),"symmetricMatrix"),"unpackedMatrix"),
             "dpoMatrix")
      })


## ~~~~ COERCIONS FROM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("dpoMatrix", "dppMatrix", function(from) pack(from))

## MJ: no longer needed ... prefer above
if(FALSE) {
setAs("dpoMatrix", "dppMatrix",
      function(from) {
          ## FIXME: check=FALSE
          copyClass(.Call(dsyMatrix_as_dspMatrix, from), "dppMatrix",
                    sNames = c("x", "Dim", "Dimnames", "uplo", "factors"))
      })
} ## MJ

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
setAs("dpoMatrix", "lMatrix",
      function(from) as(as(from, "dsyMatrix"), "lMatrix"))
setAs("dpoMatrix", "nMatrix",
      function(from) as(as(from, "dsyMatrix"), "nMatrix"))
} ## MJ


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

