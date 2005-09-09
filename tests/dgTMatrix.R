library(Matrix)

### Use ``non-unique'' versions of dgTMatrix objects

N <- 200
set.seed(1)
i <- as.integer(round(runif (N, 0, 100)))
j <- as.integer(3* rpois (N, lam=15))
x <- round(rnorm(N), 2)
which(duplicated(cbind(i,j))) # 8 index pairs are duplicated

m1 <- new("dgTMatrix", Dim = c(max(i)+1:1, max(j)+1:1), i = i, j = j, x = x)
mc <- as(m1, "dgCMatrix")
m2 <- as(mc, "dgTMatrix")## the same as 'm1' but without duplicates

stopifnot(!isTRUE(all.equal(m1, m2)),
          all.equal(as(m1,"matrix"), as(m2,"matrix"), tol=1e-15),
          all.equal(crossprod(m1), crossprod(m2), tol=1e-15),
          identical(mc, as(m2, "dgCMatrix")))


if(FALSE)
uniqify_gT <- function(x)
{
    ## Purpose: produce a *unique* triplet representation:
    ##		by having (i,j) sorted and unique
    ## ----------------------------------------------------------------------
    ## Arguments: a "gT" Matrix
    stopifnot(is(x, "gTMatrix"))
    ii <- order(x@i, x@j)
    x@i <- x@i[ii]
    x@j <- x@j[ii]
    x@x <- x@x[ii]
    if(any(duplicated(cbind(x@i, x@j))))
        sum.the.x.etc() ## UNFINISHED - FIXME
    ### We should use an exported utility for this which uses  .Call(.)
}

uniq2 <- function(x) {
    stopifnot(is(x, "gTMatrix"))
    if(is(x,"dgTMatrix")) as(as(x, "dgCMatrix"), "dgTMatrix")
    else if(is(x,"lgTMatrix")) as(as(x, "lgCMatrix"), "lgTMatrix")
    else stop("not implemented for class", class(x))
}

(t2 <- system.time(um2 <- uniq2(m1)))

if(FALSE) {
 t1 <- system.time(um1 <- uniqify_gT(m1))
 stopifnot(identical(um1, m2),
           identical(um2, m2))
}
