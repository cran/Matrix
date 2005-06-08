### Testing positive definite matrices

library(Matrix)


h9 <- Hilbert(9)
stopifnot(c(0,0) == dim(Hilbert(0)),
          c(9,9) == h9)
str(h9)
all.equal(determinant(h9)$modulus, -96.7369450737858, tol= 1e-15)
##-> TRUE only on some platforms; seen relative difference of 10^-8
stopifnot(0 == length(h9@factors))# nothing yet
str(f9 <- as(chol(h9), "dtrMatrix"))
## h9 now has factorization
stopifnot(names(h9@factors) == "Cholesky",
          all.equal(rcond(h9), 9.0938e-13),
          all.equal(rcond(f9), 9.1272e-7, tol = 1e-6))# more precision fails
str(h9)# has 'rcond' and 'factors'
options(digits=4)
(cf9 <- crossprod(f9))# looks the same as  h9 :
stopifnot(all.equal(as.matrix(h9),
                    as.matrix(cf9), tol= 1e-15))

str(hp9 <- as(h9, "dppMatrix"))

s9 <- solve(hp9, seq(nrow(hp9)))
signif(t(s9)/10000, 4)# only rounded numbers are platform-independent
(I9 <- hp9 %*% s9)
stopifnot(all.equal(cbind(1:9), as.matrix(I9), tol = 2e-9))

