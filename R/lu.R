setMethod("expand", signature(x = "denseLU"),
	  function(x, ...) .Call(LU_expand, x))

setMethod("solve", signature(a = "denseLU", b = "missing"),
	  function(a, b, ...) {
	      ll <- expand(a) #-> list(L, U, P); orig  x = P %*% L %*% U
	      ## too expensive: with(lapply(ll, solve), U %*% L %*% P)
	      solve(ll$U, solve(ll$L, ll$P))
	  })

