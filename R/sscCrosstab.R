sscCrosstab <- function(flist, upper = TRUE)
{
    .Call("sscCrosstab", lapply(as(flist, "list"),
                                function(x) as(x, "factor")[drop = TRUE]),
          as(upper, "logical")[1])
}


"groupedPerm" <- function(i, j)
{
    fr <- subset(data.frame(i = as.integer(i), j = as.integer(j)), i < j)
    p1 <- character(0)
    p2 <- character(0)
    ti <- table(fr$i)
    ni <- names(ti)
    while (nrow(fr)) {
        jj <- names(which.max(table(subset(fr,
                                          i %in%
                                          names(which(ti == min(ti))))$j)))
        p2 <- c(jj, p2)
        fr <- subset(fr, !(j %in% p2))
        ti <- table(fr$i)
        p1 <- c(p1, setdiff(ni, names(ti)))
        ni <- names(ti)
    }
    iperm <- perm <- as.integer(c(p1,p2))
    iperm[perm+1] <- seq(along=perm) - 1
    as.integer(iperm)
}

setMethod("image", signature(x = "sscCrosstab"),
          function(x, ...) {
              x <- as(x, "tripletMatrix")
              callGeneric()
          })
