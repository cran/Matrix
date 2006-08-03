#### "pedigree" class methods

## Simple constructor; main point are the 'as.*' -> prettier calls
pedigree <- function(sire, dam, label) {
    n <- length(sire)
    stopifnot(n == length(dam), n == length(label))
    sire <- as.integer(sire); dam <- as.integer(dam)
    sire[sire < 1 | sire > n] <- NA
    dam[dam < 1 | dam > n] <- NA
    new("pedigree", sire = sire, dam = dam,
	label = as.character(label))
}

setAs("pedigree", "dtCMatrix", # representation as T^{-1}
      function(from) {
	  sire <- from@sire
	  n <- length(sire)
	  animal <- seq(along = sire)
	  j <- c(sire, from@dam)
	  ind <- !is.na(j)
	  as(new("dtTMatrix", i = rep.int(animal, 2)[ind] - 1:1,
		 j = j[ind] - 1:1, x = rep.int(-0.5, sum(ind)),
		 Dim = c(n,n), Dimnames = list(from@label, NULL),
		 uplo = "L", diag = "U"), "dtCMatrix")
      })

## these data frames are now storage efficient but print less nicely
setAs("pedigree", "data.frame",
      function(from)
      data.frame(sire = from@sire, dam = from@dam,
		 row.names = from@label))

ped2DF <- function(x) {
    lab <- x@label
    lev <- seq(along = lab)
    data.frame(sire = factor(x@sire, levels = lev, labels = lab),
	       dam  = factor(x@dam,  levels = lev, labels = lab),
	       row.names = lab)
}

setMethod("show", signature(object = "pedigree"),
	  function(object) print(ped2DF(object)))

setMethod("head", "pedigree", function(x, ...)
	  do.call("head", list(x = ped2DF(x), ...)))

setMethod("tail", "pedigree", function(x, ...)
	  do.call("tail", list(x = ped2DF(x), ...)))

setMethod("chol", "pedigree",
          function(x, pivot, LINPACK) {
              ttrans <- solve(t(as(x, "dtCMatrix")))
              .Call(pedigree_chol, x,
                    as(diagU2N(t(ttrans)), "dtCMatrix"))
          })

