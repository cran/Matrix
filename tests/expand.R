### Testing expansions of factorizations

library(Matrix)

(m1 <- round(Matrix(rnorm(25), 5), 2))
(lul <- expand(lu(m1)))
stopifnot(all.equal(as(m1, "matrix"),
                    as(lul$P %*% (lul$L %*% lul$U), "matrix")))
