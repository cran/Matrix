#### "corMatrix" (was "correlation" in 2005) ---
#### ----------- correlation matrices, inheriting from  "dpoMatrix"

## dpo* -> cor* is in ./dpoMatrix.R
.M.2cor <- function(from) as(as(from, "dpoMatrix"), "corMatrix")

setAs("Matrix", "corMatrix", .M.2cor)
setAs("matrix", "corMatrix", .M.2cor)
