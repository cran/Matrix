setMethod("t", signature(x = "dtCMatrix"),
          function(x) .Call("tsc_transpose", x, PACKAGE = "Matrix"),
          valueClass = "dtCMatrix")

setAs("dtCMatrix", "ltCMatrix",
      function(from) new("ltCMatrix", i = from@i, p = from@p,
                         uplo = from@uplo, diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("matrix", "dtCMatrix",
      function(from) as(as(from, "dtTMatrix"), "dtCMatrix"))

setAs("dtCMatrix", "dgCMatrix",
      function(from) {
          if(from@diag == "U") { ## add diagonal of 1's
              ##FIXME: do this smartly - directly {in C or R}
              as(as(from, "dgTMatrix"), "dgCMatrix")
          }
          else
              new("dgCMatrix",
                  i = from@i, p = from@p, x = from@x,
                  Dim = from@Dim, Dimnames = from@Dimnames)
      })

setAs("dtCMatrix", "dgTMatrix",
      function(from)
      .Call("tsc_to_dgTMatrix", from, PACKAGE = "Matrix"))

setAs("dtCMatrix", "dgeMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgeMatrix"))

## These are all needed because cholmod doesn't support triangular:
## (see end of ./Csparse.R )
setAs("dtCMatrix", "dtTMatrix",
      function(from) {# and this is not elegant:
          x <- as(from, "dgTMatrix")
          ## FIXME: if(from@diag == "U") should drop diagonal entries:
          new("dtTMatrix", x = x@x, i = x@i, j = x@j,
              Dim = x@Dim, Dimnames = x@Dimnames,
              uplo = from@uplo, diag = "N")
      })

setAs("dtCMatrix", "TsparseMatrix", function(from) as(from, "dtTMatrix"))
setAs("dtCMatrix", "dtrMatrix",
      function(from) as(as(from, "dtTMatrix"), "dtrMatrix"))
