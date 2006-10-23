library(Matrix)

### Matrix Products including  cross products

source(system.file("test-tools.R", package = "Matrix"))

m5 <- 1 + as(diag(-1:4)[-5,], "dgeMatrix")
## named dimnames:
dimnames(m5) <- list(Rows= LETTERS[1:5], paste("C", 1:6, sep=""))
m. <- as(m5, "matrix")
stopifnot(dim(m5) == 5:6,
          class(cm5 <- crossprod(m5)) == "dpoMatrix")
assert.EQ.mat((c.m5 <- t(m5) %*% m5), as(cm5, "matrix"))
## crossprod() with numeric vector RHS and LHS
## not sensical for tcrossprod() because of 'vec' --> cbind(vec) promotion:
assert.EQ.mat( crossprod(rep(1,5), m5),  rbind( colSums(m5)))
assert.EQ.mat( crossprod(rep(1,5), m.),  rbind( colSums(m5)))
assert.EQ.mat( crossprod(m5, rep(1,5)),  cbind( colSums(m5)))
assert.EQ.mat( crossprod(m., rep(1,5)),  cbind( colSums(m5)))

## classes differ
tc.m5 <- m5 %*% t(m5)    # "dge*", no dimnames (FIXME)
(tcm5 <- tcrossprod(m5)) # "dpo*"  w/ dimnames
assert.EQ.mat(tc.m5, mm5 <- as(tcm5, "matrix"))
## tcrossprod(x,y) :
assert.EQ.mat(tcrossprod(m5, m5), mm5)
assert.EQ.mat(tcrossprod(m5, m.), mm5)
assert.EQ.mat(tcrossprod(m., m5), mm5)

## simple cases with 'scalars' treated as 1x1 matrices:
d <- Matrix(1:5)
d %*% 2
10 %*% t(d)
assertError(3 %*% d)             # must give an error , similar to
assertError(5 %*% as.matrix(d))  # -> error

## right and left "numeric" and "matrix" multiplication:
(p1 <- m5 %*% c(10, 2:6))
(p2 <- c(10, 2:5) %*% m5)
(pd1 <- m5 %*% diag(1:6))
(pd. <- m5 %*% Diagonal(x = 1:6))
(pd2 <- diag (10:6)        %*% m5)
(pd..<- Diagonal(x = 10:6) %*% m5)
stopifnot(dim(crossprod(t(m5))) == c(5,5),
          c(class(p1),class(p2),class(pd1),class(pd2),
            class(pd.),class(pd..)) == "dgeMatrix")
assert.EQ.mat(p1, cbind(c(20,30,33,38,54)))
assert.EQ.mat(pd1, m. %*% diag(1:6))
assert.EQ.mat(pd2, diag(10:6) %*% m.)
assert.EQ.mat(pd., as(pd1,"matrix"))
assert.EQ.mat(pd..,as(pd2,"matrix"))

## check that 'solve' and '%*%' are inverses
set.seed(1)
A <- Matrix(rnorm(25), nc = 5)
y <- rnorm(5)
all.equal((A %*% solve(A, y))@x, y)
Atr <- new("dtrMatrix", Dim = A@Dim, x = A@x, uplo = "U")
all.equal((Atr %*% solve(Atr, y))@x, y)

## sparse matrix products

data(KNex); mm <- KNex$mm
M <- mm[1:500, 1:200]
MT <- as(M, "TsparseMatrix")
cpr <- t(mm) %*% mm
showMethods("%*%", class=class(M))

v1 <- rep(1, ncol(M))
str(r <-  M %*% Matrix(v1))
str(rT <- MT %*% Matrix(v1))
stopifnot(identical(r, rT))
str(r. <- M %*% cbind(v1))
stopifnot(identical4(r, r., rT, M %*% as(v1, "matrix")))

v2 <- rep(1,nrow(M))
r2 <- t(Matrix(v2)) %*% M
r2T <- v2 %*% MT
str(r2. <- v2 %*% M)
stopifnot(identical4(r2, r2., rbind(v2) %*% M, t(as(v2, "matrix")) %*% M))


###--- "logical" Matrices : ---------------------

## Robert's Example, a bit more readable
fromTo <- rbind(c(2,10),
                c(3, 9))
N <- 10
nrFT <- nrow(fromTo)
rowi <- rep.int(1:nrFT, fromTo[,2]-fromTo[,1] + 1) - 1:1
coli <- unlist(lapply(1:nrFT, function(x) fromTo[x,1]:fromTo[x,2])) - 1:1

## "n" --- nonzero pattern Matrices
sM  <- new("ngTMatrix", i = rowi, j=coli, Dim=as.integer(c(N,N)))
sM # nice

sm <- as(sM, "matrix")
sM %*% sM
assert.EQ.mat(sM %*% sM,        sm %*% sm)
assert.EQ.mat(t(sM) %*% sM,
              (t(sm) %*% sm) > 0, tol=0)
crossprod(sM)
tcrossprod(sM)
stopifnot(identical(as( crossprod(sM), "ngCMatrix"), t(sM) %*%   sM),
          identical(as(tcrossprod(sM), "ngCMatrix"),  sM  %*% t(sM)))

assert.EQ.mat( crossprod(sM),  crossprod(sm) > 0)
assert.EQ.mat(tcrossprod(sM), as(tcrossprod(sm),"matrix") > 0)

## "l" --- logical Matrices -- use usual 0/1 arithmetic
nsM <- sM
sM  <- as(sM, "lMatrix")
sm <- as(sM, "matrix")
stopifnot(identical(sm, as.matrix(nsM)))
sM %*% sM
assert.EQ.mat(sM %*% sM,        sm %*% sm)
assert.EQ.mat(t(sM) %*% sM,
              t(sm) %*% sm, tol=0)
crossprod(sM)
tcrossprod(sM)
stopifnot(identical( crossprod(sM), as(t(sM) %*% sM, "dsCMatrix")),
          identical(tcrossprod(sM), as(sM %*% t(sM), "dsCMatrix")))
assert.EQ.mat( crossprod(sM),  crossprod(sm))
assert.EQ.mat(tcrossprod(sM), as(tcrossprod(sm),"matrix"))


## A sparse example - with *integer* matrix:
M <- Matrix(cbind(c(1,0,-2,0,0,0,0,0,2.2,0),
                  c(2,0,0,1,0), 0, 0, c(0,0,8,0,0),0))
t(M)
(-4:5) %*% M
stopifnot(as.vector(print(t(M %*% 1:6))) ==
          c(as(M,"matrix") %*% 1:6))
(M.M <- crossprod(M))
MM. <- tcrossprod(M)
stopifnot(class(MM.) == "dsCMatrix",
          class(M.M) == "dsCMatrix")


## even simpler
m <- matrix(0, 4,7); m[c(1, 3, 6, 9, 11, 22, 27)] <- 1
(mm <- Matrix(m))
(cm <- Matrix(crossprod(m)))
stopifnot(identical(crossprod(mm), cm))
(tm1 <- Matrix(tcrossprod(m))) #-> had bug in 'Matrix()' !
(tm2 <- tcrossprod(mm))
stopifnot(class(tm1) == class(tm2),
	  class(tm1) == "dsCMatrix")# but they differ by "uplo"

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

