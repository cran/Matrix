#### "corMatrix" (was "correlation" in 2005) ---
#### ----------- correlation matrices, inheriting from  "dpoMatrix"

.dpo2cor <- function(from) {
    if(!is.null(r <- from@factors$correlation))
        return(r)
    sd <- sqrt(diag(from, names = FALSE))
    Is <- Diagonal(x = 1 / sd)
    r <- new("corMatrix", Dim = from@Dim, Dimnames = from@Dimnames,
             uplo = from@uplo, x = (Is %*% from %*% Is)@x, sd = sd)
    .set.factors(from, "correlation", r)
}
.M2cor <- function(from) .dpo2cor(as(from, "dpoMatrix"))

setAs("dpoMatrix", "corMatrix", .dpo2cor)
setAs("dppMatrix", "corMatrix", .M2cor)
setAs(   "matrix", "corMatrix", .M2cor)
setAs(   "Matrix", "corMatrix", .M2cor)

## The 'setAs' call below is necessary to override the _implicitly defined_
## dsy->cor coercion (see ?setAs). Without it, we get:
##
## > selectMethod("coerce", c("dsyMatrix", "corMatrix"))
## function (from, to)
## {
##     obj <- new("corMatrix")
##     as(obj, "dsyMatrix") <- from
##     obj
## }
##
## which is incorrect!
setAs("dsyMatrix", "corMatrix", .M2cor)
rm(.M2cor)

## MJ: no longer needed ... prefer variant above which is much faster
if(FALSE) {
.dpo2cor <- function(from) {
    if(!is.null(cm <- from@factors$correlation))
        return(cm)
    sd <- sqrt(diag(from))
    if(is.null(names(sd)) && !is.null(nms <- from@Dimnames[[1L]]))
        names(sd) <- nms
    Is <- Diagonal(x = 1 / sd)
    .set.factors(from, "correlation",
                 new("corMatrix",
                     as(forceSymmetric(Is %*% from %*% Is), "dpoMatrix"),
                     sd = unname(sd)))
}
} ## MJ
