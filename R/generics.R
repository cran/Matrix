det <- function(x, ...) UseMethod("det")
det.Matrix <-
expand <- function(x, ...) UseMethod("expand")
expand.default <- function(x, ...) x

facmul<- function(x, factor, y, transpose = F, left = T, ...) UseMethod("facmul")

lu <- function(x, ...) UseMethod("lu")

norm <- function(x, ...) UseMethod("norm")
norm.Matrix <- function(x, type = "M")
    .Call("R_LapackPP_norm", x, as.character(type))
norm.default <- function(x, type = "M")
    .Call("R_LapackPP_norm", as.matrix(x), as.character(type))

rcond <- function(x, ...) UseMethod("rcond")

schur <- function(x, ...) UseMethod("schur")

unpack <- function(x, ...) UseMethod("unpack")
unpack.default <- function(x, ...) x
