#### Testing  cbind() & rbind()
if(FALSE)
library(Matrix, lib="/u/maechler/R/Pkgs/Matrix.Rcheck")
library(Matrix)

if(paste(R.version$major, R.version$minor, sep=".") < "2.2")
 q('no')

## else : R 2.2.0 and later --- and when Matrix was built with R >= 2.2.0:

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

