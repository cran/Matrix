setMethod("expand", signature(x = "LU"),
          function(x, ...) .Call("LU_expand", x))
