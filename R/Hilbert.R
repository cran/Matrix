Hilbert <- function(n)
{   ## generate the Hilbert matrix of dimension n
    i <- 1:n
    new("poMatrix", x = c(1/outer(i - 1, i, "+")), Dim = as.integer(c(n,n)))
}
