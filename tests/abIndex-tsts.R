#### Testing consistency  of  "abIndex" == "abstract-indexing vectors" class :
library(Matrix)

source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

validObject(ab <- new("abIndex"))
str(ab)

set.seed(1)
ex. <- list(2:1000, 0:10, sample(100), c(-3:40, 20:70),
            c(1:100,77L, 50:40, 10L), c(17L, 3L*(12:3)))
## we know which kinds will come out: "compressed" for all but random:
rD <- "rleDiff"; kinds <- c(rD,rD,"int32", rD, rD, rD)
isCmpr <- kinds == rD
ab. <- lapply(ex., as, Class = "abIndex")
nu. <- lapply(ab., as, Class = "numeric")
in. <- lapply(ab., as, Class = "integer")
rles <- lapply(ab.[isCmpr], function(u) u@rleD@rle)
r.x <-  lapply(ex.[isCmpr], function(.) rle(diff(.)))

stopifnot(sapply(ab., validObject),
          identical(ex., nu.),
          identical(ex., in.),
          ## Check that the relevant cases really *are* "compressed":
          sapply(ab., slot, "kind") == kinds,
          ## Using rle(diff(.)) is equivalent to using our C code:
          identical(rles, r.x),
          ## Checking Group Methods - "Summary" :
          sapply(ab., range) == sapply(ex., range),
          sapply(ab., any) == sapply(ex., any),
          TRUE)

for(n in 1:120) {
    cat(".")
    k <- 1 + 4*rpois(1, 5) # >= 1
    ## "random" rle -- NB: consecutive values *must* differ (for uniqueness)
    v <- as.integer(1+ 10*rnorm(k))
    while(any(dv <- duplicated(v)))
        v[dv] <- v[dv] + 1L
    rl <- structure(list(lengths = as.integer(1 + rpois(k, 10)), values  = v),
                    class = "rle")
    ai <- new("abIndex", kind = "rleDiff",
              rleD = new("rleDiff", first = rpois(1, 20), rle = rl))
    validObject(ai)
    ii <- as(ai, "numeric")
    stopifnot(is.numeric(ii), ii == round(ii),
              identical(ai, as(ii, "abIndex")))
    if(n %% 40 == 0) cat(n,"\n")
}

cat('Time elapsed: ', (.pt <- proc.time()),'\n') # "stats"
