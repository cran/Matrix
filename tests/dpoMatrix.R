### Testing positive definite matrices

library(Matrix)


h9 <- Hilbert(9)
stopifnot(c(0,0) == dim(Hilbert(0)),
          c(9,9) == dim(h9))
str(h9)
all.equal(c(determinant(h9)$modulus), -96.7369456, tol= 2e-8)
stopifnot(0 == length(h9@factors))# nothing yet
round(ch9 <- chol(h9), 3) ## round() preserves 'triangular' !
str(f9 <- as(chol(h9), "dtrMatrix"))
## h9 now has factorization
stopifnot(names(h9@factors) == "Cholesky",
          all.equal(rcond(h9), 9.0938e-13),
          all.equal(rcond(f9), 9.1272e-7, tol = 1e-6))# more precision fails
str(h9)# has 'factors'
options(digits=4)
(cf9 <- crossprod(f9))# looks the same as  h9 :
stopifnot(all.equal(as.matrix(h9),
                    as.matrix(cf9), tol= 1e-15))

str(hp9 <- as(h9, "dppMatrix"))

s9 <- solve(hp9, seq(nrow(hp9)))
signif(t(s9)/10000, 4)# only rounded numbers are platform-independent
(I9 <- hp9 %*% s9)
m9 <- matrix(1:9, dimnames = list(NULL,NULL))
stopifnot(all.equal(m9, as.matrix(I9), tol = 2e-9))

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

