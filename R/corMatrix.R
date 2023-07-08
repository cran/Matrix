## METHODS FOR CLASS: p?corMatrix
## dense correlation matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.dpo2cor <- function(from) {
    if(!is.null(to <- from@factors$correlation))
        return(to)
    sd <- sqrt(diag(from, names = FALSE))

    to <- new("corMatrix")
    to@Dim <- d <- from@Dim
    to@Dimnames <- from@Dimnames
    to@uplo <- from@uplo
    to@sd <- sd

    n <- d[1L]
    x <- from@x / sd / rep(sd, each = n)
    x[indDiag(n)] <- 1
    to@x <- x

    .set.factor(from, "correlation", to)
}

.dpp2pcor <- function(from) {
    if(!is.null(to <- from@factors$correlation))
        return(to)
    sd <- sqrt(diag(from, names = FALSE))

    to <- new("pcorMatrix")
    to@Dim <- d <- from@Dim
    to@Dimnames <- from@Dimnames
    to@uplo <- uplo <- from@uplo
    to@sd <- sd

    n <- d[1L]
    u <- uplo == "U"
    if(u) {
        r <- seq_len(n)
        s <- 1L
    } else {
        r <- seq.int(to = 1L, by = -1L, length.out = n)
        s <- seq_len(n)
    }
    x <-  from@x / rep.int(sd, r) / sd[sequence.default(r, s)]
    x[indDiag(n, upper = u, packed = TRUE)] <- 1
    to@x <- x

    .set.factor(from, "correlation", to)
}

.M2cor <- function(from) .dpo2cor(as(from, "dpoMatrix"))

.M2pcor <- function(from) .dpp2pcor(as(from, "dppMatrix"))

setAs("dpoMatrix",  "corMatrix", .dpo2cor)
setAs(   "Matrix",  "corMatrix", .M2cor)
setAs(   "matrix",  "corMatrix", .M2cor)

setAs("dppMatrix", "pcorMatrix", .dpp2pcor)
setAs(   "Matrix", "pcorMatrix", .M2pcor)
setAs(   "matrix", "pcorMatrix", .M2pcor)

if(TRUE) {
## Needed to bypass S4 quirk/bug ... without it we see the behaviour below (??)
setAs("dsyMatrix",  "corMatrix", .M2cor)
setAs("dspMatrix", "pcorMatrix", .M2pcor)
} else {
library(Matrix)
body(selectMethod("coerce", c("dsyMatrix", "corMatrix")))
## .dpo2cor(as(from, "dpoMatrix"))
as(new("dsyMatrix"), "corMatrix")
## 0 x 0 Matrix of class "corMatrix"
## <0 x 0 matrix>
body(selectMethod("coerce", c("dsyMatrix", "corMatrix")))
## {
##     obj <- new("corMatrix")
##     as(obj, "dsyMatrix") <- from
##     obj
## }
}

rm(.M2cor, .M2pcor)
