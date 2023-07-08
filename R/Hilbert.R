## Generate the Hilbert matrix of dimension 'n' :
Hilbert <- function(n) {
    n <- as.integer(n)
    i <- seq_len(n)
    new("dpoMatrix", Dim = c(n, n), x = c(1/outer(i - 1L, i, `+`)))
}
