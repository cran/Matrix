#determinant <- function(x, ...) UseMethod("determinant")
eigen <- function(x, ...) UseMethod("eigen")
expand <- function(x, ...) UseMethod("expand")
expand.default <- function(x, ...) x

facmul<- function(x, factor, y, transpose = FALSE, left = TRUE, ...)
    UseMethod("facmul")

lu <- function(x, ...) UseMethod("lu")

norm <- function(x, ...) UseMethod("norm")
norm.Matrix <- function(x, type = "M", ...)
    .Call("R_LapackPP_norm", x, as.character(type), PACKAGE="Matrix")
norm.default <- function(x, type = "M", ...)
    .Call("R_LapackPP_norm", as.matrix(x),
          as.character(type), PACKAGE="Matrix")

rcond <- function(x, ...) UseMethod("rcond")
rcond.Matrix <- function(x, type = "O", ...)
    .Call("R_LapackPP_rcond", x, as.character(type), PACKAGE="Matrix")
rcond.default <- function(x, type = "O", ...)
    .Call("R_LapackPP_rcond", as.matrix(x),
          as.character(type), PACKAGE="Matrix")

schur <- function(x, ...) UseMethod("schur")

unpack <- function(x, ...) UseMethod("unpack")
unpack.default <- function(x, ...) x

asObject <- function(x, cl) {class(x) <- as.character(cl); x}

prependClass <- function(x, cl) {class(x) <- c(as.character(cl), class(x)); x}
