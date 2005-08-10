## Utilities for the Harwell-Boeing and MatrixMarket formats

readI <- function(conn, nlines, fmt)
{
    if (!grep("\\\([[:digit:]]+I[[:digit:]]+\\\)", fmt))
        stop("Not a valid integer format")
    Iind <- regexpr('I', fmt)
    nper <- as.integer(substr(fmt, regexpr('\\\(', fmt) + 1, Iind - 1))
    iwd <- as.integer(substr(fmt, Iind + 1, regexpr('\\\)', fmt) - 1))
    expnd <- gsub(paste('(', paste(rep('.', iwd), collapse = ''), ')',
                      sep = ''), "\\1 ",
                paste(substr(readLines(conn, nlines, ok = FALSE), 1, nper *
                        iwd), collapse = ''))
    scan(textConnection(expnd), what = integer(0), quiet = TRUE)
}

readF <- function(conn, nlines, fmt)
{
    if (!grep("\\\([[:digit:]]+[DEFG][[:digit:]]+\\\.[[:digit:]]+\\\)", fmt))
        stop("Not a valid floating point format")
    Iind <- regexpr('[DEFG]', fmt)
    nper <- as.integer(substr(fmt, regexpr('\\\(', fmt) + 1, Iind - 1))
    iwd <- as.integer(substr(fmt, Iind + 1, regexpr('\\\.', fmt) - 1))
    expnd <- gsub(paste('(', paste(rep('.', iwd), collapse = ''), ')',
                      sep = ''), "\\1 ",
                paste(substr(readLines(conn, nlines, ok = FALSE), 1, nper *
                        iwd), collapse = ''))
    scan(textConnection(expnd), what = double(0), quiet = TRUE)
}

readHB <- function(filename)
{
    ## _FIXME: Modify this to accept a connection
    if (!file.exists(filename))
        stop(paste("file:", filename, "does not exist"))
    ff <- file(filename, "r")
    hdr <- readLines(ff, 4, ok = FALSE)
    Title <- sub('[[:space:]]+$', '', substr(hdr[1], 1, 72))
    Key <- sub('[[:space:]]+$', '', substr(hdr[1], 73, 80))
    totln <- as.integer(substr(hdr[2], 1, 14))
    ptrln <- as.integer(substr(hdr[2], 15, 28))
    indln <- as.integer(substr(hdr[2], 29, 42))
    valln <- as.integer(substr(hdr[2], 43, 56))
    rhsln <- as.integer(substr(hdr[2], 57, 70))
    if (!(t1 <- substr(hdr[3], 1, 1)) %in% c('C', 'R', 'P'))
        stop(paste("Invalid storage type:", t1))
    if (t1 != 'R') stop("Only numeric sparse matrices allowed")
    ## _FIXME: Patterns should also be allowed
    if (!(t2 <- substr(hdr[3], 2, 2)) %in% c('H', 'R', 'S', 'U', 'Z'))
        stop(paste("Invalid storage format:", t2))
    if (!(t3 <- substr(hdr[3], 3, 3)) %in% c('A', 'E'))
        stop(paste("Invalid assembled indicator:", t3))
    nr <- as.integer(substr(hdr[3], 15, 28))
    nc <- as.integer(substr(hdr[3], 29, 42))
    nz <- as.integer(substr(hdr[3], 43, 56))
    nel <- as.integer(substr(hdr[3], 57, 70))
    ptrfmt <- sub('[[:space:]]+$', '', substr(hdr[4], 1, 16))
    indfmt <- sub('[[:space:]]+$', '', substr(hdr[4], 17, 32))
    valfmt <- sub('[[:space:]]+$', '', substr(hdr[4], 33, 52))
    rhsfmt <- sub('[[:space:]]+$', '', substr(hdr[4], 53, 72))
    if (rhsln > 0) {
        h5 <- readLines(ff, 1, ok = FALSE)
    }        
    ptr <- readI(ff, ptrln, ptrfmt)
    ind <- readI(ff, indln, indfmt)
    vals <- readF(ff, valln, valfmt)
    close(ff)
    if (t2 == 'S')
        new("dsCMatrix", uplo = "L", p = ptr - 1:1,
            i = ind - 1:1, x = vals, Dim = c(nr, nc))
    else 
        new("dgCMatrix", p = ptr - 1:1,
            i = ind - 1:1, x = vals, Dim = c(nr, nc))

}

readMM <- function(filename)
    .Call("Matrix_readMatrixMarket", as.character(filename),
          PACKAGE = "Matrix")

