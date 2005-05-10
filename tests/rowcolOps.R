library(Matrix)

set.seed(321)
mm <- Matrix(round(rnorm(1000), 2), 50, 20)
m1 <- as(mm, "matrix")
stopifnot(all.equal(colMeans(mm), colMeans(m1)),
          all.equal(colSums(mm), colSums(m1)),
          all.equal(rowMeans(mm), rowMeans(m1)),
          all.equal(rowSums(mm), rowSums(m1)))
