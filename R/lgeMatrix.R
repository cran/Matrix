setAs("matrix", "lgeMatrix",
      function(from) {
	  new("lgeMatrix",
	      x = as.logical(from),
	      Dim = as.integer(dim(from)),
	      Dimnames =
	      if(!is.null(dn <- dimnames(from))) dn else list(NULL,NULL)
	      )
      })

setAs("lgeMatrix", "matrix",
      function(from) {
	  array(from@x, dim = from@Dim, dimnames = from@Dimnames)
      })


## dense |-> compressed :
setAs("lgeMatrix", "lgTMatrix",
      function(from) {
          ##  cheap but not so efficient:
          ij <- which(as(from,"matrix"), arr.ind = TRUE) - 1:1
          new("lgTMatrix", i = ij[,1], j = ij[,2],
              Dim = from@Dim, Dimnames = from@Dimnames,
              factors = from@factors)
      })

setMethod("as.vector", signature(x = "lgeMatrix", mode = "missing"),
          function(x) x@x)
