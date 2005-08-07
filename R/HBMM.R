## Utilities for the Harwell-Boeing and MatrixMarket formats

readHB <- function(filename)
    .Call("Matrix_readHarwellBoeing", as.character(filename),
          PACKAGE = "Matrix")

readMM <- function(filename)
    .Call("Matrix_readMatrixMarket", as.character(filename),
          PACKAGE = "Matrix")

