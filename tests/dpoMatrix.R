### Testing positive definite matrices

library(Matrix)

stopifnot(c(0,0) == dim(Hilbert(0)))

h9 <- Hilbert(9)
str(h9)
all.equal(determinant(h9)$modulus, -96.7369450737858, tol= 1e-15)
stopifnot(0 == length(h9@factors))# nothing yet
str(f9 <- as(chol(h9), "dtrMatrix"))
## h9 now has factorization
stopifnot(names(h9@factors) == "Cholesky")
stopifnot(all.equal(rcond(h9), 9.0938e-13))
stopifnot(all.equal(rcond(f9), 9.1272e-7, tol = 1e-6))# more precision fails
str(h9)# has 'rcond' and 'factors'
options(digits=4)
(cf9 <- crossprod(f9))# looks the same as  h9 :
stopifnot(all.equal(as(h9, "matrix"),
                    as(cf9,"matrix"), tol= 1e-15))

