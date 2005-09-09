library(Matrix)

### Matrix Products including  cross products

source(system.file("test-tools.R", package = "Matrix"))

m5 <- 1 + as(diag(-1:4)[-5,], "dgeMatrix")
## named dimnames:
dimnames(m5) <- list(Rows= LETTERS[1:5], paste("C", 1:6, sep=""))
stopifnot(dim(m5) == 5:6,
          class(cm5 <- crossprod(m5)) == "dpoMatrix")
assert.EQ.mat((c.m5 <- t(m5) %*% m5), as(cm5, "matrix"))
## but the 'dimnames' are not the same (and are *both*) wrong -- FIXME

## right and left "numeric" and "matrix" multiplication:
(p1 <- m5 %*% c(10, 2:6))
(p2 <- c(10, 2:5) %*% m5)
(pd1 <- m5 %*% diag(1:6))
(pd2 <- diag(10:6) %*% m5)
stopifnot(dim(crossprod(t(m5))) == c(5,5),
          c(class(p1),class(p2),class(pd1),class(pd2)) == "dgeMatrix"
          )
assert.EQ.mat(p1, cbind(c(20,30,33,38,54)))
assert.EQ.mat(pd1, as(m5,"matrix") %*% diag(1:6))
assert.EQ.mat(pd2, diag(10:6) %*% as(m5,"matrix"))

## <FIXME>  -- maybe several things:
M <- mm[1:500, 1:200] # dgT*

showMethods("%*%", class=class(M)) ## really is empty ``at first''

v1 <- rep(1, ncol(M))
if(FALSE) ## infinite recursion !
r <- M %*% Matrix(v1)
if(FALSE) ## infinite recursion !
r2 <- Matrix(rep(1,nrow(M))) %*% M
try(r3 <- M %*% cbind(v1))# dispatches to old ("S3") default method!
## same result as {but cbind() slightly differs: has dimnames!}:
try(r3. <- M %*% as(v1, "matrix"))

## </FIXME>

###--- "logical" Matrices : ---------------------

## Robert's Example, a bit more readable
fromTo <- rbind(c(2,10),
                c(3, 9))
N <- 10
nrFT <- nrow(fromTo)
rowi <- rep.int(1:nrFT, fromTo[,2]-fromTo[,1] + 1) - 1:1
coli <- unlist(lapply(1:nrFT, function(x) fromTo[x,1]:fromTo[x,2])) - 1:1
sM  <- new("lgTMatrix", i = rowi, j=coli, Dim=as.integer(c(N,N)))
sM # nice

sm <- as(sM, "matrix")
assert.EQ.mat(sM %*% sM,        sm %*% sm)
all.equal(as(t(sM) %*% sM, "matrix"),
          (t(sm) %*% sm) > 0, tol=0)

## is it intended that these are so different?
crossprod(sM)
t(sM) %*% sM

## These two are *not* TRUE !
## assert.EQ.mat(crossprod(sM),    crossprod(sm))
## assert.EQ.mat(tcrossprod(sM),   tcrossprod(sm))


proc.time() # for ``statistical reasons''
