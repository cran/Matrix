setMethod("expand", signature(x = "denseLU"),
	  function(x, ...) .Call(LU_expand, x))
