setAs("matrix", "pdmatrix",
      function(from) new("pdmatrix", from))

setAs("pdmatrix", "corrmatrix",
      function(from) {
          ss = sqrt(diag(from))
          new("corrmatrix", t(from/ss)/ss, stdDev = ss)
      })
