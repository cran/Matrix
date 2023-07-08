## METHODS FOR CLASS: dppMatrix
## dense (packed) symmetric positive semidefinite matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.dsp2dpp <- function(from) {
    if(is.null(tryCatch(Cholesky(from, perm = FALSE),
                        error = function(e) NULL)))
        stop("not a positive definite matrix (and positive semidefiniteness is not checked)")
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
