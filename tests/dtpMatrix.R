### triangular packed
library(Matrix)

cp6 <- chol(Hilbert(6))
tp6 <- as(cp6,"dtpMatrix")
round(tp6, 3)## round() is "Math2" group method
1/tp6        ## "Arith" group : gives 'dgeMatrix'
str(tp6)
## arithmetic with a mix of dimnames / no dimnames
tp <- tp6; dimnames(tp) <- list(LETTERS[1:6], letters[11:16])
stopifnot(all(tp - tp6 == tp6 - tp),
          0 == as.matrix(tp - tp6)) # FIXME: fails w/o 'as.matrix'

stopifnot(validObject(tp6),
          all.equal(tp6 %*% diag(6), as(tp6, "dgeMatrix")),
          validObject(tp6. <- diag(6) %*% tp6),
          class((tt6 <- t(tp6))) == "dtpMatrix",
          identical(t(tt6), tp6),
          tp6@uplo == "U" && tt6@uplo == "L")

all.equal(as(tp6.,"matrix"),
          as(tp6, "matrix"), tol= 1e-15)
(tr6 <- as(tp6, "dtrMatrix")) ## prints using wrong class name
D. <- determinant(tp6)
rc <- rcond(tp6)
stopifnot(all.equal(c(D.$modulus), -6.579251212),
          all.equal(rc, 1.791511257e-4),
          all.equal(norm(tp6, "I") , 2.45),
          all.equal(norm(tp6, "1") , 1),
          all.equal(norm(tp6, "F") , 1.37047826623)
          )
object.size(tp6)
object.size(as(tp6, "dtrMatrix"))
object.size(as(tp6, "matrix"))
D6 <- as(diag(6), "dgeMatrix")
ge6 <- as(tp6, "dgeMatrix")
stopifnot(all.equal(D6 %*% tp6, ge6),
          all.equal(tp6 %*% D6, ge6))

## larger case
set.seed(123)
rl <- new("dtpMatrix", uplo="L", diag="N", Dim = rep.int(1000:1000,2),
          x = rnorm(500*1001))
validObject(rl)
str(rl)
sapply(c("I", "1", "F"), function(type) norm(rl, type=type))
rcond(rl)# 0 !
stopifnot(all.equal(as(rl %*% diag(1000),"matrix"),
                    as(rl, "matrix")))
object.size(rl) ## 4 MB
object.size(as(rl, "dtrMatrix"))# 8 MB
object.size(as(rl, "matrix"))# ditto
determinant(rl)


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
