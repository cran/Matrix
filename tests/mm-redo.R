library(Matrix)

### Rebuild the 'mm' example matrix
### Use this if classes are changed
data(KNex)
mmT <- as(KNex$mm, "dgTMatrix")
str(mmT)
mm3 <- cbind(i = mmT@i, j = mmT@j, x = mmT@x)
write.table(mm3, file = "mm-Matrix.tab", row.names=FALSE)# -> ASCII version

str(mmr <- read.table("mm-Matrix.tab", header = TRUE))
mmr$i <- as.integer(mmr$i)
mmr$j <- as.integer(mmr$j)

mmN <- with(mmr, new("dgTMatrix", Dim = c(max(i)+1:1,max(j)+1:1),
                     i = i, j = j, x = x))

stopifnot(identical(mmT, mmN)) # !!
## weaker (and hence TRUE too):
stopifnot(all.equal(as(mmN, "matrix"),
                    as(mmT, "matrix"), tol=0))

mm <- as(mmN, "dgCMatrix")
stopifnot(all.equal(mm, KNex$mm))
## save(mm, file = "....../Matrix/data/mm.rda", compress = TRUE)
