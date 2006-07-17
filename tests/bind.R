#### Testing  cbind() & rbind()

library(Matrix)

source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

### --- Dense Matrices ---

m1 <- m2 <- m <- Matrix(1:12, 3,4)
dimnames(m2) <- list(LETTERS[1:3],
                     letters[1:4])
dimnames(m1) <- list(NULL,letters[1:4])

stopifnot(identical(cbind ( m, 10*m) -> R,
                    cbind2( m, 10*m))); R
stopifnot(identical(cbind (m1,100+m1) -> R,
                    cbind2(m1,100+m1))); R
stopifnot(identical(cbind (m1, 10*m2) -> R,
                    cbind2(m1, 10*m2))); R

## TODO: m1+m2 "warning" - improve dimnames() automatism
stopifnot(identical(cbind (m2, m1+m2) -> R,
                    cbind2(m2, m1+m2))); R

cbind(m2, 10*m2[nrow(m2):1 ,])# keeps the rownames from the first

(im <- cbind(I = 100, m))
str(im)
(mi <- cbind(m2, I = 1000))
str(mi)
(m1m <- cbind(m,I=100,m2))

### --- Sparse Matrices ---

m <- Matrix(c(0, 0, 2:0), 3, 5)
(mC <- as(m, "dgCMatrix"))
(mT <- as(m, "dgTMatrix"))
stopifnot(identical(mT, as(mC, "dgTMatrix")))
cbind(0, mC)
cbind(0, mT)
cbind(diag(3), mT)
(cc <- cbind(mC, 0,7,0, diag(3), 0))
stopifnot(identical3(cc, cbind(mT, 0,7,0, diag(3), 0),
                     as( cbind(m, 0,7,0, diag(3), 0), "dgCMatrix")))

cbind(mC, 1, 100*mC, 0, 0:2)
cbind(mT, 1, 0, mT+10*mT, 0, 0:2)

## print() / show() of  non-structural zeros:
(m <- Matrix(c(0, 0, 2:0), 3, 5))
(m2 <- cbind(m,m))
(m4 <- rbind(m2,m2))
diag(m4)

for(i in 1:6) {
    m4[i, i ] <- i
    m4[i,i+1] <- 0
}
m4 ## now show some non-structural zeros:


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
