## For both 'Extract' ("[") and 'Replace' ("[<-") Method testing

library(Matrix)

source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

### Dense Matrices

m <- Matrix(1:28 +0, nrow = 7)
validObject(m)
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
m[integer(0),] #-> 0 x 4 Matrix
m[2:4, numeric(0)] #-> 3 x 0 Matrix

## logical indexing
stopifnot(identical(m[2,3], m[(1:nrow(m)) == 2, (1:ncol(m)) == 3]),
          identical(m[2,], m[(1:nrow(m)) == 2, ]),
          identical(m[,3:4], m[, (1:4) >= 3]))

## dimnames indexing:
mn <- m
dimnames(mn) <- list(paste("r",letters[1:nrow(mn)],sep=""),
                     LETTERS[1:ncol(mn)])
mn["rd", "D"]
stopifnot(identical(mn["rc", "D"], mn[3,4]), mn[3,4] == 24,
          identical(mn[, "A"], mn[,1]), mn[,1] == 1:7,
          identical(mn[c("re", "rb"), "B"], mn[c(5,2), 2])
          )

mo <- m
m[2,3] <- 100
m[1:2, 4] <- 200
m[, 1] <- -1
m[1:3,]

m. <- as.matrix(m)

## m[ cbind(i,j) ] indexing:
ij <- cbind(1:6, 2:3)
stopifnot(identical(m[ij], m.[ij]))

## testing operations on logical Matrices rather more than indexing:
g10 <- m [ m > 10 ]
stopifnot(18 == length(g10))
stopifnot(10 == length(m[ m <= 10 ]))
sel <- (20 <  m) & (m <  150)
sel.<- (20 <  m.)& (m.<  150)
nsel <-(20 >= m) | (m >= 150)
(ssel <- as(sel, "sparseMatrix"))
stopifnot(is(sel, "lMatrix"), is(ssel, "lsparseMatrix"),
	  identical3(as.mat(sel.), as.mat(sel), as.mat(ssel)),
	  identical3(!sel, !ssel, nsel), # !<sparse> is typically dense
	  identical3(m[ sel],  m[ ssel], as.matrix(m)[as.matrix( ssel)]),
	  identical3(m[!sel],  m[!ssel], as.matrix(m)[as.matrix(!ssel)])
	  )

## more sparse Matrices --------------------------------------

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
assert.EQ.mat(mC[1:2,], mm[1:2,])
stopifnot(all.equal(mC[,3], mm[,3]),
	  identical(mC[ij], mm[ij]))
assert.EQ.mat(mC[7, , drop=FALSE], mm[7, , drop=FALSE])

stopifnot(dim(mC[numeric(0), ]) == c(0,20), # used to give warnings
          dim(mC[, integer(0)]) == c(40,0),
          identical(mC[, integer(0)], mC[, FALSE]),
          identical(mC[7,  drop = FALSE],
                    mC[7,, drop = FALSE]))
validObject(print(mT[,c(2,4)]))
stopifnot(all.equal(mT[2,], mm[2,]),
          ## row or column indexing in combination with t() :
          identical(mT[2,], t(mT)[,2]),
          identical(mT[-2,], t(t(mT)[,-2])),
          identical(mT[c(2,5),], t(t(mT)[,c(2,5)]))
          )
assert.EQ.mat(mT[4,, drop = FALSE], mm[4,, drop = FALSE])
stopifnot(identical3(mm[,1], mC[,1], mT[,1]),
	  identical3(mm[3,], mC[3,], mT[3,]),
	  identical3(mT[2,3], mC[2,3], 0),
	  identical(mT[], mT),
          identical4(       mm[c(3,7), 2:4],  as.mat( m[c(3,7), 2:4]),
                     as.mat(mT[c(3,7), 2:4]), as.mat(mC[c(3,7), 2:4]))
          )

x.x <- crossprod(mC)
stopifnot(class(x.x) == "dsCMatrix",
          class(x.x. <- round(x.x / 10000)) == "dsCMatrix",
          identical(x.x[cbind(2:6, 2:6)],
                    diag(x.x [2:6, 2:6])))
head(x.x.) # Note the *non*-structural 0's printed as "0"
tail(x.x., -3) # all but the first three lines

lx.x <- as(x.x, "lsCMatrix") # FALSE only for "structural" 0
(l10 <- lx.x[1:10, 1:10])# "lsC"
(l3 <-  lx.x[1:3, ])
m.x <- as(x.x, "matrix")
stopifnot(class(l10) == "lsCMatrix", # symmetric indexing -> symmetric !
          identical(as.mat(lx.x), m.x != 0),
          identical(as.logical(lx.x), as.logical(m.x)),
          identical(as.mat(l10), m.x[1:10, 1:10] != 0),
          identical(as.mat(l3 ), m.x[1:3, ] != 0)
          )

## used to fail
n <- 5 ## or much larger
sm <- new("dsTMatrix", i=as.integer(1),j=as.integer(1),
          Dim=as.integer(c(n,n)), x = 1)
(cm <- as(sm, "CsparseMatrix"))
sm[2,]
stopifnot(sm[2,] == c(0:1, rep.int(0,ncol(sm)-2)),
	  sm[2,] == cm[2,],
	  sm[,3] == sm[3,],
	  all(sm[,-(1:3)] == t(sm[-(1:3),])), # all(<lge.>)
	  all(sm[,-(1:3)] == 0)
	  )

### Diagonal -- Sparse:
m0 <- Diagonal(5)
(m1 <- as(m0, "sparseMatrix"))  # dtTMatrix
(m2 <- as(m0, "CsparseMatrix")) # dtCMatrix (with an irrelevant warning)

M <- m0; M[1,] <- 0
stopifnot(identical(M, Diagonal(x=c(0, rep(1,4)))))
M <- m0; M[,3] <- 3 ; M ; stopifnot(is(M, "sparseMatrix"), M[,3] == 3)
validObject(M)
M <- m0; M[1:3, 3] <- 0 ;M
T <- m0; T[1:3, 3] <- 10
stopifnot(identical(M, Diagonal(x=c(1,1, 0, 1,1))),
          is(T, "triangularMatrix"), identical(T[,3], c(10,10,10,0,0)))

M <- m1; M[1,] <- 0 ; M ; assert.EQ.mat(M, diag(c(0,rep(1,4))), tol=0)
M <- m1; M[,3] <- 3 ; stopifnot(is(M,"sparseMatrix"), M[,3] == 3)
validObject(M)
M <- m1; M[1:3, 3] <- 0 ;M
assert.EQ.mat(M, diag(c(1,1, 0, 1,1)), tol=0)
T <- m1; T[1:3, 3] <- 10; validObject(T)
stopifnot(is(T, "dtTMatrix"), identical(T[,3], c(10,10,10,0,0)))

M <- m2; M[1,] <- 0 ; M ; assert.EQ.mat(M, diag(c(0,rep(1,4))), tol=0)
M <- m2; M[,3] <- 3 ; stopifnot(is(M,"sparseMatrix"), M[,3] == 3)
validObject(M)
M <- m2; M[1:3, 3] <- 0 ;M
assert.EQ.mat(M, diag(c(1,1, 0, 1,1)), tol=0)
T <- m2; T[1:3, 3] <- 10; validObject(T)
stopifnot(is(T, "dtCMatrix"), identical(T[,3], c(10,10,10,0,0)))



## --- negative indices ----------
mc <- mC[1:5, 1:7]
mt <- mT[1:5, 1:7]
## sub matrix
assert.EQ.mat(mC[1:2, 0:3], mm[1:2, 0:3]) # test 0-index
stopifnot(identical(mc[-(3:5), 0:2], mC[1:2, 0:2]),
          identical(mt[-(3:5), 0:2], mT[1:2, 0:2]),
          identical(mC[2:3, 4],      mm[2:3, 4]))
assert.EQ.mat(mC[1:2,], mm[1:2,])
## sub vector
stopifnot(identical4(mc[-(1:4), ], mC[5, 1:7],
                     mt[-(1:4), ], mT[5, 1:7]))
stopifnot(identical4(mc[-(1:4), -(2:4)], mC[5, c(1,5:7)],
                     mt[-(1:4), -(2:4)], mT[5, c(1,5:7)]))

## mixing of negative and positive must give error
assertError(mT[-1:1,])

## Sub *Assignment* ---- now works (partially):
mt0 <- mt
mt[1, 4] <- -99
mt[2:3, 1:6] <- 0
mt
m2 <- mt+mt
m2[1,4] <- -200
m2[c(1,3), c(5:6,2)] <- 1:6
stopifnot(m2[1,4] == -200,
          as.vector(m2[c(1,3), c(5:6,2)]) == 1:6)
mt[,3] <- 30
mt[2:3,] <- 250
mt[1:5 %% 2 == 1, 3] <- 0
mt[3:1, 1:7 > 5] <- 0
mt

tt <- as(mt,"matrix")
ii <- c(0,2,5)
jj <- c(2:3,5)
tt[ii, jj] <- 1:6 # 0 is just "dropped"
mt[ii, jj] <- 1:6
assert.EQ.mat(mt, tt)

mt[1:5, 2:6]
as((mt0 - mt)[1:5,], "dsparseMatrix")# [1,5] and lines 2:3

mt[c(2,4), ] <- 0; stopifnot(as(mt[c(2,4), ],"matrix") == 0)
mt[2:3, 4:7] <- 33
validObject(mt)
mt

mc[1,4] <- -99 ; stopifnot(mc[1,4] == -99)
mc[1,4] <-  00 ; stopifnot(mc[1,4] ==  00)
mc[1,4] <- -99 ; stopifnot(mc[1,4] == -99)
mc[1:2,4:3] <- 4:1; stopifnot(as.matrix(mc[1:2,4:3]) == 4:1)

mc[-1, 3] <- -2:1 # 0 should not be entered; 'value' recycled
mt[-1, 3] <- -2:1
stopifnot(mc@x != 0, mt@x != 0,
	  mc[-1,3] == -2:1, mt[-1,3] == -2:1) ## failed earlier

mc0 <- mc
set.seed(1)
for(i in 1:20) {
    mc <- mc0
    ev <- 1:5 %% 2 == round(runif(1))# 0 or 1
    j <- sample(ncol(mc), 1 + round(runif(1)))
    nv <- rpois(sum(ev) * length(j), lambda = 1)
    mc[ev, j] <- nv
    if(i < 5) print(mc[ev,j, drop = FALSE])
    stopifnot(as.vector(mc[ev, j]) == nv) ## failed earlier...
    validObject(mc)
}

mc # no longer has non-structural zeros
mc[ii, jj] <- 1:6
mc[c(2,5), c(3,5)] <- 3.2
validObject(mc)
m. <- mc
mc[4,] <- 0
mc

H <- Hilbert(9)
Hc <- as(round(H, 3), "dsCMatrix")# a sparse matrix with no 0 ...
(trH <- tril(Hc[1:5, 1:5]))
stopifnot(is(trH, "triangularMatrix"), trH@uplo == "L")

i <- c(1:2, 4, 6:7); j <- c(2:4,6)
H[i,j] <- 0
(H. <- round(as(H, "sparseMatrix"), 3)[ , 2:7])
Hc. <- Hc
Hc.[i,j] <- 0 ## now "works", but setting "non-structural" 0s
stopifnot(as.matrix(Hc.[i,j]) == 0)
Hc.[, 1:6]

## an example that failed long
sy3 <- new("dsyMatrix", Dim = as.integer(c(2, 2)), x = c(14, -1, 2, -7))
validObject(dm <- kronecker(Diagonal(2), sy3))
(s2 <- as(dm, "sparseMatrix"))
validObject(st <- as(s2, "TsparseMatrix"))
validObject(s.32  <- st[1:3,1:2]) ## 3 x 2 - and *not* dsTMatrix
validObject(s2.32 <- s2[1:3,1:2])
I <- c(1,4:3)
stopifnot(is(s2.32, "generalMatrix"),
          is(s.32,  "generalMatrix"),
          identical(as.mat(s.32), as.mat(s2.32)),
          identical3(dm[1:3,-1], asD(s2[1:3,-1]), asD(st[1:3,-1])),
          identical4(2, dm[4,3], s2[4,3], st[4,3]),
          identical3(diag(dm), diag(s2), diag(st)),
          is((cI <- s2[I,I]), "dsCMatrix"),
          is((tI <- st[I,I]), "dsTMatrix"),
          identical4(as.mat(dm)[I,I], as.mat(dm[I,I]), as.mat(tI), as.mat(cI))
          )

## now sub-assign  and check for consistency
## symmetric subassign should keep symmetry
st[I,I] <- 0; validObject(st); stopifnot(is(st,"symmetricMatrix"))
s2[I,I] <- 0; validObject(s2); stopifnot(is(s2,"symmetricMatrix"))

m <- as.mat(st)
 m[2:1,2:1] <- 4:1
st[2:1,2:1] <- 4:1
s2[2:1,2:1] <- 4:1
stopifnot(identical(m, as.mat(st)),
	  1:4 == as.vector(s2[1:2,1:2]),
	  identical(m, as.mat(s2)))

## m[cbind(i,j)] <- value:
m.[ cbind(3:5, 1:3) ] <- 1:3
stopifnot(m.[3,1] == 1, m.[4,2] == 2)
x.x[ cbind(2:6, 2:6)] <- 12:16
validObject(x.x)
stopifnot(class(x.x) == "dsCMatrix",
	  12:16 == as.mat(x.x)[cbind(2:6, 2:6)])

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
