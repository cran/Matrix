## METHODS FOR CLASS: dpoMatrix
## dense (unpacked) symmetric positive semidefinite matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.dsy2dpo <- function(from) {
    if(is.null(tryCatch(Cholesky(from, perm = FALSE),
                        error = function(e) NULL)))
        stop("not a positive definite matrix (and positive semidefiniteness is not checked)")
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
