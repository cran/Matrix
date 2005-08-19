#### Thanks to  our as.matrix() maniplation in base namespace,  ../R/zzz.R ,
#### all the functions (in 'base' or namespaces that import) and
#### start with something like
####	" x <- as.matrix(x) "
#### will work for 'Matrix'-matrices

library(Matrix)

str(m1 <- mm[1:500, 1:200])
m11 <- m1[1:100, 1:20]
## These now work thanks to using our as.matrix():
str(D1 <- dist(m11))
str(rs <- apply(m1, 1, sum))

kappa(Matrix(2:5, 2))
## used to seg.fault, PR#7984,
## because qr() was calling the wrong as.matrix()

## also matplot() or pairs().
