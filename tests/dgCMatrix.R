library(Matrix)

data(mm)
is(mm)
stopifnot(dim(mm) == c(1850, 712))
dimnames(mm) # empty  {but currently "NULL" instead of list(NULL,NULL)}
str(mm)
tmm <- t(mm)
str(tmm)
validObject(tmm)
stopifnot(all.equal(as(tmm, "matrix"), t(as(mm, "matrix"))))
