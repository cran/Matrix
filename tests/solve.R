library(Matrix)

## Solve() methods, first for "dgeMatrix" and all kinds of RHS :

(m6 <- 1 + as(diag(0:5), "dgeMatrix"))
rcond(m6)
I6 <- as(diag(6), "dgeMatrix")
stopifnot(all.equal(I6, m6 %*% solve(m6)),
          all.equal(I6, solve(m6) %*% m6) )

(i6 <- solve(m6, Matrix(1:6)))
stopifnot(identical(i6, as(cbind(c(-4, rep(1,5))), "dgeMatrix")),
          identical(i6, solve(m6, 1:6),
          identical(i6, solve(m6, matrix(1:6))),
          identical(i6, solve(m6, matrix(c(1,2,3,4,5,6)))),
          )
