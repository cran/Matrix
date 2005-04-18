#### For both 'Extract' ("[") and 'Replace' ("[<-") Method testing

library(Matrix)

identical3 <- function(x,y,z)	identical(x,y) && identical (y,z)
identical4 <- function(a,b,c,d) identical(a,b) && identical3(b,c,d)

m <- Matrix(1:28, nrow = 7)

## TODO: not yet for dense matrices


### Sparse Matrices

m <- 1:800
set.seed(101) ; m[sample(800, 600)] <- 0
m <- Matrix(m, nrow = 40)
mm <- as(m, "matrix")
dimnames(mm) <- NULL ## << workaround: as(<sparse>, "matrix") has NULL dimnames
str(mC <- as(m, "dgCMatrix"))
str(mT <- as(m, "dgTMatrix"))
stopifnot(identical(mT, as(mC, "dgTMatrix")),
	  identical(mC, as(mT, "dgCMatrix")))

mC[,1]
mC[1:2,]
mC[7, drop = FALSE]

mT[,c(2,4)]
mT[1,]
mT[4, drop = FALSE]
stopifnot(identical3(mm[,1], mC[,1], mT[,1]),
	  identical3(mm[3,], mC[3,], mT[3,]),
	  identical3(mT[2,3], mC[2,3], 0),
	  identical(mT[], mT),
	  ## TODO: identical4() with  m[c(3,7), 2:4]
	  identical3(as(mC[c(3,7), 2:4],"matrix"), mm[c(3,7), 2:4],
		     as(mT[c(3,7), 2:4],"matrix")))


