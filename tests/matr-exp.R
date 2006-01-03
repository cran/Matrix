library(Matrix)

## Matrix Exponential

source(system.file("test-tools.R", package = "Matrix"))

## e ^ 0 = 1  - for matrices:
assert.EQ.mat(expm(Matrix(0, 3,3)), diag(3), tol = 0)# exactly
## e ^ diag(.) = diag(e ^ .):
assert.EQ.mat(expm(as(diag(-1:4), "dgeMatrix")), diag(exp(-1:4)))
set.seed(1)
rE <- replicate(100,
            { x <- rlnorm(12)
              relErr(as(expm(as(diag(x), "dgeMatrix")),
                        "matrix"),
                     diag(exp(x))) })
stopifnot(mean(rE) < 1e-15,
          max(rE)  < 1e-14)
summary(rE)

## Some small matrices

m1 <- Matrix(c(1,0,1,1), nc = 2)
e1 <- expm(m1)
assert.EQ.mat(e1, cbind(c(exp(1),0), exp(1)))

m2 <- Matrix(c(-49, -64, 24, 31), nc = 2)
e2 <- expm(m2)
## The true matrix exponential is 'te2':
e_1 <-  exp(-1)
e_17 <- exp(-17)
te2 <- rbind(c(3*e_17 - 2*e_1, -3/2*e_17 + 3/2*e_1),
             c(4*e_17 - 4*e_1, -2  *e_17 + 3  *e_1))
assert.EQ.mat(e2, te2, tol = 1e-13)
## See the (average relative) difference:
all.equal(as(e2,"matrix"), te2, tol = 0) # 1.48e-14 on "lynne"

## The ``surprising identity''      det(exp(A)) == exp( tr(A) )
## or                           log det(exp(A)) == tr(A) :
stopifnot(all.equal(c(determinant(e2)$modulus), sum(diag(m2))))

m3 <- Matrix(cbind(0,rbind(6*diag(3),0)), nc = 4)#  sparse
e3 <- expm(m3)
E3 <- expm(Matrix(m3, sparse=FALSE))
stopifnot(identical(e3, E3))
e3. <- rbind(c(1,6,18,36),
             c(0,1, 6,18),
             c(0,0, 1, 6),
             c(0,0, 0, 1))
assert.EQ.mat(e3, e3.)

proc.time() # for ``statistical reasons''
