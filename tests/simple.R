#### Currently a collection of simple tests
##	(since 'Matrix' takes long to load, rather have fewer source files!)

library(Matrix)

source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

options(verbose = TRUE)# to show message()s

### Matrix() ''smartness''
(d4 <- Matrix(diag(4)))
(z4 <- Matrix(0*diag(4)))
(o4 <- Matrix(1+diag(4)))
(tr <- Matrix(cbind(1,0:1)))
(m4 <- Matrix(cbind(0,rbind(6*diag(3),0))))
dm4 <- Matrix(m4, sparse = FALSE)
class(mN <-  Matrix(NA, 3,4)) # NA *is* logical
stopifnot(validObject(d4), validObject(z4), validObject(o4),
          validObject(m4), validObject(dm4), validObject(mN))
assert.EQ.mat(dm4, as(m4, "matrix"))
assert.EQ.mat(mN, matrix(NA, 3,4))
sL <- Matrix(, 3,4, sparse=TRUE)# -> "lgC"
trS <- Matrix(tr, sparse=TRUE)# failed in 0.9975-11
stopifnot(is(d4, "diagonalMatrix"),   is(z4,  "diagonalMatrix"),
          is(tr, "triangularMatrix"), is(trS, "triangularMatrix"),
          all(is.na(sL@x)), ## not yet:  all(is.na(sL)),
          !any(sL, na.rm=TRUE), all(!sL, na.rm=TRUE),
          validObject(Matrix(c(NA,0), 4, 3, byrow = TRUE)),
          validObject(Matrix(c(NA,0), 4, 4)),
          is(Matrix(c(NA,0,0,0), 4, 4), "sparseMatrix"))

## large sparse ones: these now directly "go sparse":
str(m0 <- Matrix(0,     nrow=100, ncol = 1000))
str(l0 <- Matrix(FALSE, nrow=100, ncol = 200))
stopifnot(all(!l0),
          identical(FALSE, any(l0)))

## with dimnames:
m. <- matrix(c(0, 0, 2:0), 3, 5)
dimnames(m.) <- list(LETTERS[1:3], letters[1:5])
(m0 <- m <- Matrix(m.))
m@Dimnames[[2]] <- m@Dimnames[[1]]
## not valid anymore:
(val <- validObject(m, test=TRUE))
stopifnot(is.character(val))
dm <- as(m0, "denseMatrix")
stopifnot(identical(rcond(dm), rcond(as.matrix(dm))),
	  all.equal(rcond(dm), 0.4899474520656))
rm(m)

###--  Sparse Triangular :

g5 <- new("dgCMatrix", Dim = c(5L, 5L),
          x = c(10, 1, 3, 10, 1, 10, 1, 10, 10),
          i = c(0L,2L,4L, 1L, 3L,2L,4L, 3L, 4L),
          p = c(0L, 3L, 5L, 7:9))
t5 <- as(g5, "triangularMatrix") # works fine (but slowly) FIXME
stopifnot(class(t5) == "dtCMatrix",
          identical(t5, tril(g5)))

(t1 <- new("dtTMatrix", x= c(3,7), i= 0:1, j=3:2,
           Dim= as.integer(c(4,4))))
## Diagonal  o  Sparse
I4 <- Diagonal(4)
validObject(t2  <-   t1  + I4)
validObject(tt2 <- t(t1) + I4)
validObject(t1c <- as(t1, "CsparseMatrix"))
validObject(t2c <- as(t2, "CsparseMatrix"))
stopifnot(validObject(t1),
          identical(t1, t(t(t1))),
          identical(t1c, t(t(t1c))),
##           is(t2,"triangularMatrix"), is(t2c,"triangularMatrix"),
##           is(t1c,"triangularMatrix"),
          is(t1c + I4,"triangularMatrix"), is(t2c + I4,"triangularMatrix"),
          c(class(t2), class(t1c), class(t2c), class(tt2)) == "dtCMatrix",
          identical(t(tt2), t2))
assert.EQ.mat(t1, as(t1c, "matrix"))

## as(<diag>, <anything>) :
str(cls <- names(getClass("Matrix")@subclasses))# all Matrix classes
D4 <- Diagonal(4, x=1:4)

for(cl in cls)
    if(canCoerce(I4, cl)) {
	cat(cl,":")
	M  <- as(I4, cl)
	M. <- as(D4, cl)
        stopifnot(diag(4) == as(M,"matrix"),
                  if(is(cl,"dMatrix")) diag(x=1:4) == as(M.,"matrix") else TRUE)
	cat(" [Ok]\n")
    }

## from  0-diagonal to unit-diagonal {low-level step}:
tu <- t1 ; tu@diag <- "U"
tu
cu <- as(tu, "dtCMatrix")
validObject(cnu <- Matrix:::diagU2N(cu))# <- testing diagU2N
stopifnot(validObject(cu), validObject(tu. <- as(cu, "dtTMatrix")),
          validObject(tt <- as(cu, "TsparseMatrix")),
	  ## NOT: identical(tu, tu.), # since T* is not unique!
	  identical(cu, as(tu., "dtCMatrix")),
          length(cnu@i) == length(cu@i) + nrow(cu),
          identical(cu, Matrix:::diagN2U(cnu)),# <- testing diagN2U
	  all(cu >= 0, na.rm = TRUE), all(cu >= 0),
	  any(cu >= 7),
	  validObject(tcu <- t(cu)),
	  validObject(ttu <- t(tu)))

assert.EQ.mat(cu, as(tu,"matrix"), tol=0)
assert.EQ.mat(cnu, as(tu,"matrix"), tol=0)

## <sparse> o <numeric> (of length > 1):
stopifnot(is(tm <- tu * 1:8, "sparseMatrix"),
          identical4(tm, cu * 1:8, 1:8 * cu, 1:8 * tu))

cu[1,2] <- tu[1,2] <- NA
mu <- as(tu,"matrix")
stopifnot(is(cu, "CsparseMatrix"), is(cu, "triangularMatrix"),
          is(tu, "TsparseMatrix"), is(tu, "triangularMatrix"),
          identical(cu * 1:8, tu * 1:8), # but are no longer triangular
          all(cu >= 0, na.rm=TRUE), !all(cu >= 1), is.na(all(tu >= 0)))
assert.EQ.mat(cu * 1:8, mu * 1:8)

## tu. is diag "U", but tu2 not:
tu2 <- as(as(tu., "generalMatrix"), "triangularMatrix")
assert.EQ.mat(cu, mu, tol=0)
stopifnot(identical3(cu[cu > 1],  tu [tu > 1], mu [mu > 1]),
          identical3(cu <= 1, tu <= 1, as(mu <= 1, "lMatrix")),# all lgeMatrix
	  identical3(cu[cu <= 1], tu[tu <= 1], mu[mu <= 1]),
	  identical3(cu , triu(cu ), t(t(cu))),
	  identical3(tu , triu(tu ), t(t(tu))),
	  identical3(tu., triu(tu.), t(t(tu.))),
	  identical(tu2, triu(tu2)),
	  identical(tcu , tril(tcu)),
	  identical(ttu , tril(ttu)),
	  identical(t(tu), tril(t(tu)))
          )
assert.EQ.mat(triu(cu),   as.matrix(triu(as.matrix(cu))))
for(k in -1:1)
    assert.EQ.mat(tril(cu,k), as.matrix(tril(as.matrix(cu),k)))

(dtr <- Matrix(local({m <- diag(2); m[1,2] <- 3;m})))
identical(dtr, triu(dtr))
assert.EQ.mat(tril(dtr), diag(2))


(t4 <- new("dgTMatrix", i = 3:0, j = 0:3, x = rep(1,4), Dim = as.integer(c(4,4))))
c4 <- as(t4, "CsparseMatrix")
## the same but "dsT" (symmetric)
M <- Matrix(c(0, rep(c(0,0:1),4)), 4,4)#ok warning
tt <- as(M, "TsparseMatrix")
stopifnot(all.equal(triu(t4) + tril(t4), c4),
          all.equal(triu(tt) + tril(tt), c4))


###-- Numeric Dense: Crossprod & Solve

set.seed(123)
mm. <- mm <- Matrix(rnorm(500 * 150), nc = 150)
stopifnot(validObject(mm))
xpx <- crossprod(mm)
stopifnot(identical(mm, mm.),# once upon a time, mm was altered by crossprod()
          validObject(xpx))
str(mm) # 'dge*"
str(xpx)# 'dpo*"
xpy <- crossprod(mm, rnorm(500))
res <- solve(xpx, xpy)
str(xpx)# now with Cholesky factor
stopifnot(validObject(xpx),
          validObject(xpy),
          validObject(res))
stopifnot(all.equal(xpx %*% res, xpy, tol= 1e-12))
lp <- xpx >= 1
slp <- as(lp, "sparseMatrix")

ltlp  <-  lp[ lower.tri(lp) ]
sltlp <- slp[ lower.tri(slp) ]
dim(ij <- which(lower.tri(lp), arr.ind = TRUE))
ss <- slp[ij] # now fast (!)
stopifnot(identical4(lp[ij], ltlp, sltlp, as(lp, "matrix")[ij]),
          identical(ss, sltlp),
          is(lp, "lsyMatrix"), lp@uplo == "U")

###-- more solve() methods  {was ./solve.R }

## first for "dgeMatrix" and all kinds of RHS :
(m6 <- 1 + as(diag(0:5), "dgeMatrix"))
rcond(m6)
I6 <- as(diag(6), "dgeMatrix")
stopifnot(all.equal(I6, m6 %*% solve(m6)),
          all.equal(I6, solve(m6) %*% m6) )

(i6 <- solve(m6, Matrix(1:6)))
stopifnot(identical(i6, as(cbind(c(-4, rep(1,5))), "dgeMatrix")),
          identical(i6, solve(m6, 1:6)),
          identical(i6, solve(m6, matrix(1:6))),
          identical(i6, solve(m6, matrix(c(1,2,3,4,5,6))))
          )

## solve(<sparse>)
(m <- t1+ t(t1) + Diagonal(4))
i.m <- solve(as.mat(m))
I1 <- m %*% i.m
o4 <- diag(I1)
im <- solve(m)
(I2 <- m %*% im)
(ms <- as(m, "symmetricMatrix"))
## solve(<sparse>, <sparse>):
s.mm <-  solve(m,m)
s.mms <- solve(m, ms)
## these now work "fully-sparse"
s.ms2 <- solve(ms, ms)
s.msm <- solve(ms, m)
I4c <- as(Matrix(diag(4),sparse=TRUE), "generalMatrix")
stopifnot(is(im, "Matrix"), is(I2, "Matrix"), class(I4c) == "dgCMatrix",
          all.equal(I1, I2, tol = 1e-14),
          all.equal(diag(4), as.mat(I2), tol = 1e-12),
          all.equal(s.mm,  I2, tol = 1e-14),
          all.equal(s.mms, I2, tol = 1e-14),
          all.equal(s.ms2, s.msm, tol = 4e-15),
          all.equal(s.ms2, I4c  , tol = 4e-15),
          abs(o4 - 1) < 1e-14)

###-- row- and column operations  {was ./rowcolOps.R }

set.seed(321)
(m1 <- round(Matrix(rnorm(25), 5), 2))
m1k <- Matrix(round(rnorm(1000), 2), 50, 20)
m.m <- as(m1k, "matrix")
stopifnot(all.equal(colMeans(m1k), colMeans(m.m)),
          all.equal(colSums (m1k), colSums (m.m)),
          all.equal(rowMeans(m1k), rowMeans(m.m)),
          all.equal(rowSums (m1k), rowSums (m.m)))

###-- kronecker for nonsparse uses Matrix(.):
stopifnot(is(kr <- kronecker(m1, m6), "Matrix"))
assert.EQ.mat(kr,
              kronecker(as(m1, "matrix"),
                        as(m6, "matrix")), tol = 0)

## sparse:
(kt1 <- kronecker(t1, tu))
kt2 <- kronecker(t1c, cu)
stopifnot(identical(Matrix:::uniq(kt1), Matrix:::uniq(kt2)))
## but kt1 and kt2, both "dgT" are different since entries are not ordered!
ktf <- kronecker(as.matrix(t1), as.matrix(tu))
if(FALSE) # FIXME? our kronecker treats "0 * NA" as "0" for structural-0
assert.EQ.mat(kt2, ktf, tol= 0)
(cs1 <- colSums(kt1))
NA.or.True <- function(x) is.na(x) | x
eq <- (cs1 == colSums(as(kt1, "matrix")))
stopifnot(NA.or.True(eq), identical(is.na(eq), is.na(cs1)))
nt1 <- as(kt1, "nMatrix") # no NA's anymore
(ng1 <- as(as(nt1, "generalMatrix"),"CsparseMatrix")) # ngC
dg1 <- as(ng1, "dMatrix")# dgC
lt1 <- kt1 > 5
nt1 <- as(lt1, "nMatrix")
(colSums(nt1, sparseResult = TRUE))
(colSums(kt1, sparseResult = TRUE)) # dsparse, with NA
(colSums(lt1, sparseResult = TRUE)) # isparse, with NA
(colSums(lt1, sparseResult = TRUE, na.rm = TRUE))
(colSums(nt1, sparseResult = TRUE)) # isparse, no NA
## check correct sparseness of both:
for(M in list(kt1, nt1, ng1, dg1, lt1, nt1)) {
    m <- as(M, "matrix")
    for(na.rm in c(FALSE,TRUE)) {
	cs  <- colSums(M, na.rm = na.rm)
	cs. <- colSums(M, na.rm = na.rm, sparseResult = TRUE)
	rs  <- rowSums(M, na.rm = na.rm)
	rs. <- rowSums(M, na.rm = na.rm, sparseResult = TRUE)
	stopifnot(identical(cs, as(cs., "vector")),
		  identical(rs, as(rs., "vector")),
		  {eq <- cs == colSums(m, na.rm = na.rm) ; ineq <- is.na(eq)
		   all(ineq | eq) && identical(ineq, is.na(cs)) },
		  {eq <- rs == rowSums(m, na.rm = na.rm) ; ineq <- is.na(eq)
		   all(ineq | eq) && identical(ineq, is.na(rs)) },
		  is(cs., "sparseVector"),
		  is(rs., "sparseVector"))
    }
}

## coercion from "dpo" or "dsy"
xx <- as(xpx, "dsyMatrix")
stopifnot(isSymmetric(xxS  <- as(xx,  "sparseMatrix")),
          isSymmetric(xpxS <- as(xpx, "sparseMatrix")))

tm <- matrix(0, 8,8)
tm[cbind(c(1,1,2,7,8),
         c(3,6,4,8,8))] <- c(2,-30,15,20,80)
(tM <- Matrix(tm))                ## dtC
(mM <- Matrix(m <- (tm + t(tm)))) ## dsC
mT <- as(mM, "dsTMatrix")
gC <- as(as(mT, "dgTMatrix"), "dgCMatrix")
## Check that 'mT' and gC print properly :
pr.mT <- capture.output(mT)
nn <- unlist(strsplit(gsub(" +\\.", "", sub("^....", "", pr.mT[-(1:2)])), " "))
stopifnot(as.numeric(nn[nn != ""]) == m[m != 0],
          capture.output(gC)[-1] == pr.mT[-1])
assert.EQ.mat(tM, tm, tol=0)
assert.EQ.mat(gC, m,  tol=0)
assert.EQ.mat(mT, m,  tol=0)
stopifnot(is(mM, "dsCMatrix"), is(tM, "dtCMatrix")
	  , identical(mT, as(mM, "TsparseMatrix"))
	  , identical(gC, as(mM, "generalMatrix"))
	  ## coercions	general <-> symmetric
	  , identical(as(as(mM, "generalMatrix"), "symmetricMatrix"), mM)
	  , identical(as(as(mM, "dgTMatrix"),     "symmetricMatrix"), mT)
	  , identical(as(as(tM, "generalMatrix"),"triangularMatrix"), tM)
          , identical(tM + Diagonal(8), tMD <- Diagonal(8) + tM)
          , is(tMD, "dtCMatrix")
	  )
eM <- eigen(mM) # works thanks to base::as.matrix hack in ../R/zzz.R
stopifnot(all.equal(eM$values,
                { v <- c(162.462112512353, 30.0665927567458)
                  c(v, 15, 0, 0, 160-v[1], -15, -v[2])}, tol=1e-14))

##--- symmetric -> pos.def. needs valid test:
m5 <- Matrix(diag(5) - 1)
if(FALSE) { # FIXME: this as(.,.) happily "works"
 assertError(mpo <- as(m5, "dpoMatrix"))
 validObject(mpo) #-> TRUE  FIXME?  it is *not* really pos.definite!
}

###-- dense nonzero pattern:
class(m <- Matrix(TRUE,2,2)) # lsy
(n <- as(m, "nMatrix")) # nsy
validObject(n)

## 1)
as(n,"CsparseMatrix") # used to give CHOLMOD error: invalid xtype...
ls2 <- as(m, "CsparseMatrix") # works fine
## and really  'm' and 'n' are interally slot identical (!!!)

as(n,"sparseMatrix")
as(m,"sparseMatrix")

### -- now when starting with nsparse :
nT <- new("ngTMatrix",
          i = as.integer(c(0, 1, 0)),
          j = as.integer(c(0, 0, 1)), Dim = as.integer(c(2,2)),
          Dimnames = list(NULL, NULL))
(nC <- as(nT, "ngCMatrix"))
str(nC)# of course, no 'x' slot

stopifnot(identical(tt <- as(nT,"denseMatrix"), # lge
		    as(as(nT, "lMatrix"),"denseMatrix")))
tt
as(nC,"denseMatrix")


###-- sparse nonzero pattern : ----------

(nkt <- as(as(as(kt1, "generalMatrix"), "CsparseMatrix"), "ngCMatrix"))# ok
dkt <- as(nkt, "denseMatrix")
(clt <- crossprod(nkt))
stopifnot(is(nkt, "ngCMatrix"), is(clt, "nsCMatrix"))
crossprod(clt) ## a warning: crossprod() of symmetric

## a Csparse with *repeated* entry is not valid!
assertError(new("ngCMatrix", p = c(0L,2L), i = c(0L,0L), Dim = 2:1))


### "d" <-> "l"  for (symmetric) sparse :
data(KNex)
mm <- KNex$mm
xpx <- crossprod(mm)
## extract nonzero pattern
nxpx <- as(xpx, "nsCMatrix")
show(nxpx) ## now ok, since subsetting works
r <- nxpx[1:2,]

lmm <- as(mm, "lgCMatrix")
nmm <- as(lmm, "nMatrix")
xlx <- crossprod(lmm)
x.x <- crossprod(nmm)

## now A = lxpx and B = xlx should be close, but not quite the same
## since <x,y> = 0 is well possible when x!=0 and y!=0 .
## However,  A[i,j] != 0 ==> B[i,j] != 0:
A <- as(as(nxpx, "lMatrix"), "TsparseMatrix")
B <- as(as(xlx,  "lMatrix"), "TsparseMatrix")
ij <- function(a) a@i + ncol(a) * a@j
stopifnot(all(ij(A) %in% ij(B)))

l3 <- upper.tri(matrix(,3,3))
(c3 <- as(l3, "CsparseMatrix"))# lgC
stopifnot(validObject(c3), is(c3, "CsparseMatrix"), is(c3, "lMatrix"))
(M <- Matrix(l3))  # -> "ltCMatrix"
M2 <- M %x% M
stopifnot(is(M, "ltCMatrix"), validObject(M2), dim(M2) == c(9,9),
          identical(M2, kronecker(M,M)),
          is(M2, "triangularMatrix")) # is "dtT" (why not "dtC" ?)
M3 <- M %x% M2 #ok
(cM3 <- colSums(M3, sparse=TRUE))
identical(as.vector(cM3),
          as(rev(rowSums(M3, sparse=TRUE)), "vector"))
M. <- M2 %x% M # gave infinite recursion

## diagonal, sparse & interactions
stopifnot(is(as(Diagonal(3), "TsparseMatrix"), "TsparseMatrix"),
          is(X <- Diagonal(7) + 1.5 * tM[1:7,1:7], "sparseMatrix"),
          is(X, "triangularMatrix"),
          is(XX <- X - chol(crossprod(X)), "triangularMatrix"))
X
XX
XX <- as(drop0(XX), "dsCMatrix")
stopifnot(identical(XX, Matrix(0, nrow(X), ncol(X))))

M <- Matrix(m., sparse = FALSE)
(sM <- Matrix(m.))
class(dlM <- M >= 1)
stopifnot(identical(dlM, !(M < 1)),
	  is(sM, "sparseMatrix"),
	  is(dlM, "denseMatrix"))
(lM  <- as(dlM, "sparseMatrix"))
lM2 <- as(dlM, "CsparseMatrix") #-> now ok
lM0 <- Matrix:::as_Csparse(dlM)
stopifnot(identical3(lM, lM2, lM0))

selectMethod("coerce",	c("lgeMatrix", "CsparseMatrix"),
	     useInherited = c(from = TRUE, to = FALSE))

ms0 <- Matrix(c(0,1,1,0), 2,2)
(ms <- as(ms0, "TsparseMatrix"))
cs <- as(ms, "CsparseMatrix")
ll <- as(ms, "lMatrix")
lt <- as(ll, "lgTMatrix")
nn <- as(cs, "nsparseMatrix")
l2 <- as(cs, "lsparseMatrix")
nt <- triu(nn)
n3 <- as(nt, "lsparseMatrix")
da <- nt + t(nt)
dm <- nt * t(nt) + da
stopifnot(as(ms0,"matrix") == as(ll, "matrix"), # coercing num |-> log
	  as(lt, "matrix") == as(ll, "matrix"),
	  identical(ms, as(ll, "dMatrix")), is(ms, "dsTMatrix"),
	  identical4(as(ll, "CsparseMatrix"), as(cs, "lMatrix"),# lsC*
		     as(nn, "lsparseMatrix"), l2),
	  identical3(da, dm, as(cs, "generalMatrix")),		# dgC*
	  identical(as(da, "lMatrix"), as(lt, "CsparseMatrix")) # lgC*
	  )


cat('Time elapsed: ', proc.time(),'\n') # "stats"

if(!interactive()) warnings()

