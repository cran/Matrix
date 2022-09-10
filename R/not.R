## METHODS FOR GENERIC: ! (not)
## logical negation
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Divert everything other than [ln]Matrix to lMatrix:
setMethod("!", "Matrix", function(x) !as(x, "lMatrix"))

## -- diagonalMatrix --

setMethod("!", "ldiMatrix",
          function(x) {
              r <- new("lspMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if((n <- d[1L]) > 0L) {
                  r@x <- rep.int(TRUE, 0.5 * n * (n + 1))
                  r@x[indDiag(n = n, upper = TRUE, packed = TRUE)] <-
                      if(x@diag == "N") !x@x else FALSE
              }
              r
          })

## -- [ln]sparseMatrix -- [[ FALSE->TRUE  =>  dense result ]]

setMethod("!", "lsparseMatrix", function(x) !.sparse2dense(x))
setMethod("!", "nsparseMatrix", function(x) !.sparse2dense(x))

## -- [ln]denseMatrix --

for(.cl in paste0("l", c("ge", "sy", "sp"), "Matrix"))
    setMethod("!", .cl,
              function(x) {
                  x@x <- !x@x
                  x
              })

for(.cl in paste0("n", c("ge", "sy", "sp"), "Matrix"))
    setMethod("!", .cl,
              function(x) {
                  x@x <- !(is.na(x@x) | x@x) # NA <=> TRUE
                  x
              })

for(.cl in paste0("l", c("tr", "tp"), "Matrix"))
    setMethod("!", .cl,
              function(x) {
                  r <- .dense2g(x)
                  r@x <- !r@x
                  r
              })

for(.cl in paste0("n", c("tr", "tp"), "Matrix"))
    setMethod("!", .cl,
              function(x) {
                  r <- .dense2g(x)
                  r@x <- !(is.na(r@x) | r@x) # NA <=> TRUE
                  r
              })

### -- sparseVector --

setMethod("!", "sparseVector",
	  function(x) {
	      n <- x@length
	      if(2 * length(x@i) <= n)
		  !sp2vec(x)
	      else { ## sparse result
		  ii <- seq_len(n)[-x@i]
		  if((has.x <- !is(x, "nsparseVector"))) {
		      xx <- rep.int(TRUE, length(ii))
		      if((.na <- any(x.na <- is.na(x@x))) |
			 (.fa <- any(x.f <- !x.na & !x@x))) {
			  ## deal with 'FALSE' and 'NA' in  x slot
			  if(.na) {
			      ii <- c(ii, x@i[x.na])
			      xx <- c(xx, x@x[x.na])
			  }
			  if(.fa) { ## any(x.f)
			      x.f <- x.f & !x.na
			      ii <- c(ii, x@i[x.f])
			      xx <- c(xx, rep.int(TRUE, sum(x.f)))
			  }
			  ## sort increasing in index:
			  i.s <- sort.list(ii)
			  ii <- ii[i.s]
			  xx <- xx[i.s]
		      }
		  }
		  if(has.x)
		      newSpV("lsparseVector", x = xx, i = ii, length = n)
		  else new("nsparseVector", i = ii, length = n)
	      }
	  })
