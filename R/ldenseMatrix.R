setAs("lgeMatrix", "dgeMatrix", l2d_Matrix)
setAs("ltrMatrix", "dtrMatrix", l2d_Matrix)
setAs("ltpMatrix", "dtpMatrix", l2d_Matrix)
setAs("lsyMatrix", "dsyMatrix", l2d_Matrix)
setAs("lspMatrix", "dspMatrix", l2d_Matrix)

setAs("lspMatrix", "lsyMatrix",
      function(from)
      .Call("lspMatrix_as_lsyMatrix", from, PACKAGE = "Matrix"))

setAs("lsyMatrix", "lspMatrix",
      function(from)
      .Call("lsyMatrix_as_lspMatrix", from, PACKAGE = "Matrix"))

setAs("ltpMatrix", "ltrMatrix",
      function(from)
      .Call("ltpMatrix_as_ltrMatrix", from, PACKAGE = "Matrix"))

setAs("ltrMatrix", "ltpMatrix",
      function(from)
      .Call("ltrMatrix_as_ltpMatrix", from, PACKAGE = "Matrix"))

setAs("ldenseMatrix", "matrix",
      function(from) as(as(from, sub("^l", "d", class(from))), "matrix"))

setAs("matrix", "ldenseMatrix",
      function(from) callGeneric(as(from, "lgeMatrix")))

setMethod("t", signature(x = "lgeMatrix"), t_geMatrix)
setMethod("t", signature(x = "ltrMatrix"), t_trMatrix)
setMethod("t", signature(x = "lsyMatrix"), t_trMatrix)
setMethod("t", signature(x = "ltpMatrix"),
          function(x) as(callGeneric(as(x, "ltrMatrix")), "ltpMatrix"))
setMethod("t", signature(x = "lspMatrix"),
          function(x) as(callGeneric(as(x, "lsyMatrix")), "lspMatrix"))

setMethod("!", "ltrMatrix",
          function(e1) {
              e1@x <- !e1@x
              ## And now we must fill in the '!FALSE' results :

              ## FIXME: the following should be .Call using
              ##        a variation of make_array_triangular:
              r <- as(e1, "lgeMatrix")
              n <- e1@Dim[1]
              coli <- rep(1:n, each=n)
              rowi <- rep(1:n, n)
              Udiag <- e1@diag == "U"
              log.i <-
                  if(e1@uplo == "U") {
                      if(Udiag) rowi >= coli else rowi > coli
                  } else {
                      if(Udiag) rowi <= coli else rowi < coli
                  }
              r[log.i] <- TRUE
              r
          })

setMethod("!", "ltpMatrix", function(e1) !as(x, "ltrMatrix"))

## for the other ldense* ones:
setMethod("!", "ldenseMatrix",
          function(e1) { e1@x <- !e1@x ; e1 })
