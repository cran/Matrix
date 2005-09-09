## Utilities for the Harwell-Boeing and MatrixMarket formats

readone <- function(ln, iwd, nper, conv)
{
    ln <- gsub("D", "E", ln)
    inds <- seq(0, by = iwd, length = nper + 1)
    (conv)(substring(ln, 1 + inds[-length(inds)], inds[-1]))
}

readmany <- function(conn, nlines, nvals, fmt, conv)
{
    if (!grep("[[:digit:]]+[DEFGI][[:digit:]]+", fmt))
        stop("Not a valid format")
    Iind <- regexpr('[DEFGI]', fmt)
    nper <- as.integer(substr(fmt, regexpr('[[:digit:]]+[DEFGI]', fmt), Iind - 1))
    iwd <- as.integer(substr(fmt, Iind + 1, regexpr('[\\\.\\\)]', fmt) - 1)) 
    rem <- nvals %% nper
    full <- nvals %/% nper
    ans <- vector("list", nvals %/% nper)
    for (i in seq(len = full))
        ans[[i]] <- readone(readLines(conn, 1, ok = FALSE),
                            iwd, nper, conv)
    if (!rem) return(unlist(ans))
    c(unlist(ans),
      readone(readLines(conn, 1, ok = FALSE), iwd, rem, conv))
}

readHB <- function(file)
{
    if (is.character(file)) 
        if (file == "") 
            file <- stdin()
        else 
            file <- file(file)
    if (!inherits(file, "connection")) 
        stop("'file' must be a character string or connection")
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
    hdr <- readLines(file, 4, ok = FALSE)
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
    ptrfmt <- toupper(sub('[[:space:]]+$', '', substr(hdr[4], 1, 16)))
    indfmt <- toupper(sub('[[:space:]]+$', '', substr(hdr[4], 17, 32)))
    valfmt <- toupper(sub('[[:space:]]+$', '', substr(hdr[4], 33, 52)))
    rhsfmt <- toupper(sub('[[:space:]]+$', '', substr(hdr[4], 53, 72)))
    if (!is.na(rhsln) && rhsln > 0) {
        h5 <- readLines(file, 1, ok = FALSE)
    }        
    ptr <- readmany(file, ptrln, nc + 1, ptrfmt, as.integer)
    ind <- readmany(file, indln, nz, indfmt, as.integer)
    vals <- readmany(file, valln, nz, valfmt, as.numeric)
    if (t2 == 'S')
        new("dsCMatrix", uplo = "L", p = ptr - 1:1,
            i = ind - 1:1, x = vals, Dim = c(nr, nc))
    else 
        new("dgCMatrix", p = ptr - 1:1,
            i = ind - 1:1, x = vals, Dim = c(nr, nc))

}

readMM <- function(file)
{
    if (is.character(file)) 
        if (file == "") 
            file <- stdin()
        else 
            file <- file(file)
    if (!inherits(file, "connection")) 
        stop("'file' must be a character string or connection")
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
    if ((hdr <- scan(file, nmax = 1, what = character(0), quiet = TRUE))
        != "%%MatrixMarket")
        stop("file is not a MatrixMarket file")
    typ <- tolower(scan(file, nmax = 1, what = character(0), quiet = TRUE))
    if (!typ %in% "matrix")
        stop(paste("type '", typ, "' not recognized", sep = ""))
    repr <- tolower(scan(file, nmax = 1, what = character(0), quiet = TRUE))
    if (!repr %in% c("coordinate", "array"))
        stop(paste("representation '", repr, "' not recognized", sep = ""))
    elt <- tolower(scan(file, nmax = 1, what = character(0), quiet = TRUE))
    if (!elt %in% c("real", "complex", "integer", "pattern"))
        stop(paste("element type '", elt, "' not recognized", sep = ""))
    sym <- tolower(scan(file, nmax = 1, what = character(0), quiet = TRUE))
    if (!sym %in% c("general", "symmetric", "skew-symmetric", "hermitian"))
        stop(paste("symmetry form '", sym, "' not recognized", sep = ""))
    nr <- scan(file, nmax = 1, what = integer(0),
               comment.char = "%", quiet = TRUE)
    nc <- scan(file, nmax = 1, what = integer(0), quiet = TRUE)
    nz <- scan(file, nmax = 1, what = integer(0), quiet = TRUE)
    if (repr == "coordinate" && elt == "real") {
        els <- scan(file, nmax = nz,
                    what = list(i = integer(0), j = integer(0),
                    x = numeric(0)), quiet = TRUE)
        if (sym == "general")
            return(new("dgTMatrix", Dim = c(nr, nc), i = els$i - 1:1,
                       j = els$j - 1:1, x = els$x))
        if (sym == "symmetric")
            return(new("dsTMatrix", uplo = "L", Dim = c(nr, nc),
                       i = els$i - 1:1, j = els$j - 1:1, x = els$x))
    }
}
