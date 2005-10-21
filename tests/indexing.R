#### For both 'Extract' ("[") and 'Replace' ("[<-") Method testing

library(Matrix)

source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

### Dense Matrices

m <- Matrix(1:28, nrow = 7)
validObject(m) ; m@x <- as.double(m@x) ; validObject(m)
stopifnot(identical(m, m[]),
          identical(m[2, 3],  16), # simple number
          identical(m[2, 3:4], c(16,23))) # simple numeric of length 2

m[2, 3:4, drop=FALSE] # sub matrix of class 'dgeMatrix'
m[-(4:7), 3:4]        # dito; the upper right corner of 'm'

## rows or columns only:
m[1,]     # first row, as simple numeric vector
m[,2]     # 2nd column
m[,1:2]   # sub matrix of first two columns
m[-(1:6),, drop=FALSE] # not the first 6 rows, i.e. only the 7th

## logical indexing
stopifnot(identical(m[2,3], m[(1:nrow(m)) == 2, (1:ncol(m)) == 3]),
          identical(m[2,], m[(1:nrow(m)) == 2, ]),
          identical(m[,3:4], m[, (1:4) >= 3]))

## dimnames indexing:
mn <- m
dimnames(mn) <- list(paste("r",letters[1:nrow(mn)],sep=""),
                     LETTERS[1:ncol(mn)])
mn["rd", "D"]
stopifnot(identical(mn["rc", "D"], mn[3,4]),
          identical(mn[, "A"], mn[,1]),
          identical(mn[c("re", "rb"), "B"], mn[c(5,2), 2])
          )

mo <- m
m[2,3] <- 100
m[1:2, 4] <- 200
m[, 1] <- -1
m[1:3,]

## TODO: more --- particularly once we have "m > 10" working!


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
mC[7,  drop = FALSE]
stopifnot(identical(mC[7,  drop = FALSE],
                    mC[7,, drop = FALSE]))
mT[,c(2,4)]
mT[1,]
mT[4, drop = FALSE]
stopifnot(identical3(mm[,1], mC[,1], mT[,1]),
	  identical3(mm[3,], mC[3,], mT[3,]),
	  identical3(mT[2,3], mC[2,3], 0),
	  identical(mT[], mT),
	  ## TODO: identical4() with  m[c(3,7), 2:4] - fail because of 'dimnames'
	  identical3(as(mC[c(3,7), 2:4],"matrix"), mm[c(3,7), 2:4],
		     as(mT[c(3,7), 2:4],"matrix")))

## --- negative indices ----------
mc <- mC[1:5, 1:7]
mt <- mT[1:5, 1:7]
## sub matrix
stopifnot(identical(mc[-(3:5), 0:2], mC[1:2, 0:2]),
          identical(mt[-(3:5), 0:2], mT[1:2, 0:2]))
## sub vector
stopifnot(identical4(mc[-(1:4), ], mC[5, 1:7],
                     mt[-(1:4), ], mT[5, 1:7]))
stopifnot(identical4(mc[-(1:4), -(2:4)], mC[5, c(1,5:7)],
                     mt[-(1:4), -(2:4)], mT[5, c(1,5:7)]))

## mixing of negative and positive must give error
assertError(mT[-1:1,])

## At least these now give a nicely understandable error:
try(mT[1, 4] <- -99)
try(mT[2:3, ] <- 0)
