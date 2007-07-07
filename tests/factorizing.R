#### Matrix Factorizations  --- of all kinds

library(Matrix)

source(system.file("test-tools.R", package = "Matrix"))# identical3() etc


### "sparseQR" : Check consistency of methods
##   --------
data(KNex); mm <- KNex$mm; y <- KNex$y
stopifnot(is((Y <- Matrix(y)), "dgeMatrix"))
md <- as(mm, "matrix")                  # dense

system.time(mmq <- qr(mm))
system.time(mdq <- qr(md))# much (~ 150 x) slower

## qr.qy and qr.qty should be inverses
stopifnot(all.equal(qr.qy (mmq, qr.qty(mmq, y))@x, y),
          all.equal(qr.qty(mmq, qr.qy (mmq, y))@x, y),
          all.equal(qr.qty(mmq, y), qr.qty(mmq, Y)) )

## consistency of results dense and sparse
stopifnot(is.all.equal3(qr.coef  (mdq, y), qr.coef  (mmq,y)@x, qr.coef  (mmq,Y)@x) ,
          is.all.equal3(qr.resid (mdq, y), qr.resid (mmq,y)@x, qr.resid (mmq,Y)@x) ,
          is.all.equal3(qr.fitted(mdq, y), qr.fitted(mmq,y)@x, qr.fitted(mmq,Y)@x) )


### "denseLU"

## Testing expansions of factorizations {was ./expand.R, then in simple.R }

set.seed(1)
(m1 <- round(Matrix(rnorm(25), 5), 2))
str(lu1 <- lu(m1))
(luX <- expand(lu1))
stopifnot(all.equal(as(m1, "matrix"),
                    as(luX$P %*% (luX$L %*% luX$U), "matrix")))

### "sparseLU"
por1 <- readMM(system.file("external/pores_1.mtx", package = "Matrix"))
lu1 <- lu(por1)
pm <- as(por1, "CsparseMatrix")
(pmLU <- lu(pm)) # -> show(<MatrixFactorization>)
## identical only as long as we don't keep the original class info:
stopifnot(identical(lu1, pmLU))

## permute rows and columns of original matrix
ppm <- pm[pmLU@p + 1:1, pmLU@q + 1:1]
Ppm <- pmLU@L %*% pmLU@U
## these two should be the same, and `are' in some ways:
assert.EQ.mat(ppm, as(Ppm, "matrix"), tol = 1e-14)
## *however*
length(ppm@x)# 180
length(Ppm@x)# 317 !
table(Ppm@x == 0)# (194, 123) - has 123 "zero" and 14 ``almost zero" entries

## FIXME:  expand(pmLU)


