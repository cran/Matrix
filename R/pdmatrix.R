setAs("matrix", "pdmatrix",
      function(from) new("pdmatrix", from))

setAs("pdfactor", "pdmatrix",
      function(from) new("pdmatrix", crossprod(from)))

setAs("pdmatrix", "pdfactor",
      function(from) {
          val <- new("pdfactor",
                     .Call("nlme_Chol", as(from, "pdmatrix")))
          val@logDet <- sum(log(diag(val)))
          val
      },
      function(from, value) {
          as(from, "pdmatrix") <- crossprod(value)
          from
      })

setAs("pdmatrix", "corrmatrix",
      function(from) {
          ss = sqrt(diag(from))
          new("corrmatrix", t(from/ss)/ss, stdDev = ss)
      })
