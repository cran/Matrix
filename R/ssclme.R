setAs("ssclme", "dsCMatrix",
      function(from)
      new("dsCMatrix", i = from@i, p = from@p, Dim = from@Dim, x = from@x,
	  uplo = "U"))

setAs("ssclme", "dtCMatrix",
      function(from)
      new("dtCMatrix", i = from@Li, p = from@Lp, Dim = from@Dim, x = from@Lx,
	  uplo = "L", diag = "U"))

setReplaceMethod("coef", signature(object = "ssclme", value = "numeric"),
		 function(object, unconst = FALSE, ..., value)
		 .Call("ssclme_coefGets", object, as.double(value), unconst))

setMethod("coef", signature(object = "ssclme"),
	  function(object, unconst = FALSE, ...)
	  .Call("ssclme_coef", object, unconst) )

