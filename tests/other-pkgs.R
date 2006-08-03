####--------- Test interfaces to other non-standard Packages ---------------

###-- 1)  'graph' (from Bioconductor) ---------------------------
###-- ==  =======                     ---------------------------
if(require(graph)) {

if(packageDescription("graph")$Version <= "1.10.2") {
    ## graph 1.10.x for x <= 2  had too much problems in  as(<graph>, "matrix")
    cat("Version of 'graph' is too old --- no tests done here!\n")
    q('no')
}

## else do things

library(Matrix)

example("graphNEL-class", echo = FALSE)
nodes(gR)

## 1) undirected
sm.g <- as(gR, "sparseMatrix")
str(sm.g)## 'dsT' - fine; even has Dimnames!
validObject(sm.g)
sm.g # should show the Dimnames - at least row ones

## 2) directed
gU <- gR; edgemode(gU) <- "directed"
sgU <- as(gU, "sparseMatrix")
str(sgU)## 'dgT' with dimnames
validObject(sgU)
sgU # should now show the Dimnames!

### Reverse :  sparseMatrix -> graph
if(FALSE) { ## not yet (FIXME!)
    gmg  <- as(sm.g, "graph")
    gmg2 <- as(sm.g, "graphNEL")
    gmgU <- as(sgU, "graph")
    gmgU <- as(sgU, "graphNEL")
}

}## end{graph}

###-- 2)  'SparseM' ---------------------------------------------
###-- ==  ========  ---------------------------------------------

if(require(SparseM)) {

set.seed(1)
a <- round(rnorm(5*4), 2)
a[abs(a) < 0.7] <- 0
A <- matrix(a,5,4)
print(M <- Matrix(A))
stopifnot(
          validObject(A.csr <- as.matrix.csr(A)),
          validObject(At.csr <- as.matrix.csr(t(A))),
          identical(At.csr, t(A.csr)),
          identical(A, as.matrix(A.csr)),
          identical(M, as(A.csr, "dgCMatrix")),
          identical(t(M), as(At.csr, "dgCMatrix"))
          )

## TODO: More tests; in particular for triplets !

}## end{SparseM}
