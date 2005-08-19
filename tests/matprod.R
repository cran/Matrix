library(Matrix)

### Matrix Products including  cross products

## checking;  'show' is for convenience of the developer
assert.EQ.mat <- function(M, m, tol = if(show) 0 else 1e-15, show=FALSE) {
    ## temporary fix for R-2.0.1
    MM <- as(M, "matrix")
    attr(MM, "dimnames") <- attr(m, "dimnames") <- NULL
    if(show) all.equal(MM, m, tol = tol)
    else stopifnot(all.equal(MM, m, tol = tol))
}

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


proc.time() # for ``statistical reasons''
