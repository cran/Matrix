hilbert <- function(n)
{   ## generate the Hilbert matrix of dimension n
    i <- 1:n
    Matrix(1 / outer(i - 1, i, "+"))
}
