## METHODS FOR CLASS: sparseVector (virtual)
## sparse vectors
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("diff", signature(x = "sparseVector"),
          ## Mostly cut and paste of base::diff.default :
          function(x, lag = 1L, differences = 1L, ...) {
              if(length(lag) != 1L || length(differences) != 1L ||
                  lag < 1L || differences < 1L)
                  stop(gettextf("'%s' and '%s' must be positive integers",
                                "lag", "differences"),
                       domain = NA)
              if(lag * differences >= length(x))
                  return(x[0L])
              i1 <- -seq_len(lag)
              for(i in seq_len(differences))
                  x <- x[i1] - x[-length(x):-(length(x) - lag + 1L)]
              x
          })

setMethod("mean", signature(x = "sparseVector"),
          function(x, trim = 0, na.rm = FALSE, ...) {
              kind <- .M.kind(x)
              if(kind == "z" && trim > 0)
                  stop("trimmed means are not defined for complex data")
              n <- length(x)
              if(kind != "n" && n > 0L && anyNA(x@x)) {
                  if(!na.rm)
                      return(NA_real_)
                  n <- n - sum(is.na(x@x))
              }
              if(n == 0L)
                  return(if(kind == "z") NaN * 0i else NaN)
              if(kind == "n") {
                  nnz <- length(x@i)
                  if(trim <= 0)
                      return(nnz / n)
                  ntrim <- trunc(n * min(trim, 0.5))
                  if(nnz < ntrim)
                      0
                  else if(nnz == ntrim) {
                      if(n - 2 * ntrim > 0)
                           0
                      else 0.5
                  } else {
                      if(n - 2 * ntrim > 0)
                          (nnz - ntrim - max(ntrim - (n - nnz), 0)) /
                              (n - 2 * ntrim)
                      else 1
                  }
              } else {
                  if(trim <= 0)
                      return(sum(x@x, na.rm = na.rm) / n)
                  ntrim <- trunc(n * min(trim, 0.5))
                  x <- .V.sort(x, na.last = NA)[(ntrim + 1):(n - ntrim)]
                  sum(x@x) / length(x)
              }
          })

.V.rep.each <- function(x, each) {
    each <- as.double(each)
    if(length(each) != 1L) {
        warning(gettextf("first element used of '%s' argument", "each"),
                domain = NA)
        each <- each[1L]
    }
    if(!is.finite(each) || each <= -1)
        stop(gettextf("invalid '%s' argument", "each"), domain = NA)
    if(each < 1)
        return(x[0L])
    if(each < 2)
        return(x)
    n <- length(x)
    each <- trunc(each)
    if(n * each > 0x1p+53)
        stop(gettextf("%s length cannot exceed %s", "sparseVector", "2^53"),
             domain = NA)
    else if(n * each > .Machine$integer.max) {
        a <- as.double
        one <- 1
    } else {
        each <- as.integer(each)
        a <- as.integer
        one <- 1L
    }
    x@length <- n * each
    x@i <- rep(each * (a(x@i) - one), each = each) + seq_len(each)
    if(.M.kind(x) != "n")
        x@x <- rep(x@x, each = each)
    x
}

.V.rep.int  <- function(x, times) {
    times <- as.double(times)
    if(length(times) != 1L) {
        ## FIXME: support length(times) == length(x)
        warning(gettextf("first element used of '%s' argument", "times"),
                domain = NA)
        times <- times[1L]
    }
    if(!is.finite(times) || times <= -1)
        stop(gettextf("invalid '%s' argument", "times"), domain = NA)
    if(times < 1)
        return(x[0L])
    if(times < 2)
        return(x)
    n <- length(x)
    times <- trunc(times)
    if(n * times > 0x1p+53)
        stop(gettextf("%s length cannot exceed %s", "sparseVector", "2^53"),
             domain = NA)
    else if(n * times > .Machine$integer.max) {
        a <- as.double
        zero <- 0
    } else {
        times <- as.integer(times)
        a <- as.integer
        zero <- 0L
    }
    x@length <- n * times
    x@i <- rep(a(seq.int(from = zero, by = n, length.out = times)),
               each = length(x@i)) + x@i
    if(.M.kind(x) != "n")
        x@x <- rep.int(x@x, times)
    x
}

.V.rep.len  <- function(x, length.out) {
    length.out <- as.double(length.out)
    if(length(length.out) != 1L) {
        warning(gettextf("first element used of '%s' argument", "length.out"),
                domain = NA)
        length.out <- length.out[1L]
    }
    if(!is.finite(length.out) || length.out <= -1)
        stop(gettextf("invalid '%s' argument", "length.out"), domain = NA)
    if(length.out > 0x1p+53)
        stop(gettextf("%s length cannot exceed %s", "sparseVector", "2^53"),
             domain = NA)
    n <- length(x)
    length.out <-
        if(length.out - 1 < .Machine$integer.max)
            as.integer(length.out)
        else trunc(length.out)
    if(length.out > n && n > 0L) {
        x <- .V.rep.int(x, ceiling(length.out / n))
        n <- length(x)
    }
    x@length <- length.out
    if(length.out < n) {
        head <- x@i <= length.out
        x@i <- x@i[head]
        if(.M.kind(x) != "n")
            x@x <- x@x[head]
    } else if(length.out > n && n == 0L) {
        x@i <- seq_len(length.out)
        if(.M.kind(x) != "n")
            x@x <- rep.int(x@x[NA_integer_], length.out)
    }
    x
}

setMethod("rep", signature(x = "sparseVector"),
          function(x, times, length.out, each, ...) {
              if(!missing(each))
                  x <- .V.rep.each(x, each)
              if(!missing(length.out))
                  x <- .V.rep.len (x, length.out)
              else if(!missing(times))
                  x <- .V.rep.int (x, times)
              x
          })

.V.sort <- function(x, decreasing = FALSE, na.last = NA, ...) {
    nnz <- length(x@i)
    if(nnz == 0L)
        return(x)
    n <- length(x)
    kind <- .M.kind(x)
    if(kind == "n") {
        x@i <- if(decreasing)
                   seq_len(nnz)
               else seq.int(to = n, length.out = nnz)
        return(x)
    }
    x@x <- y <- sort.int(x@x, na.last = na.last,
                         decreasing = decreasing, ...)
    if(!is.na(na.last))
        nna <- if(anyNA(y)) sum(is.na(y)) else 0L
    else {
        x@length <- n <- n - (nnz - length(y))
        nna <- 0L
        nnz <- length(y)
    }
    nnn <- switch(kind,
                  "l" = nnz - nna,
                  "i" = sum(y >= 0L, na.rm = TRUE),
                  "d" = sum(y >= 0 , na.rm = TRUE),
                  "z" =
                      {
                          arg <- Arg(y)
                          hpi <- 0.5 * pi
                          sum(arg > -hpi & arg <= hpi, na.rm = TRUE)
                      },
                  stop("should never happen ..."))
    if(nna > 0L && decreasing != na.last)
        nnn <- nnn + nna
    x@i <-
        if(nnn < nnz) {
            if(decreasing)
                c(seq_len(nnn), seq.int(to = n, length.out = nnz - nnn))
            else
                c(seq_len(nnz - nnn), seq.int(to = n, length.out = nnn))
        } else {
            if(decreasing)
                seq_len(nnn)
            else
                seq.int(to = n, length.out = nnn)
        }
    x
}

if(FALSE) {
## MJ: once 'sort' becomes implicit generic in package 'methods' :
setMethod("sort", signature(x = "sparseVector"), .V.sort)
## TODO: parallel method for internal generic 'xtfrm'
}

setMethod("t", signature(x = "sparseVector"),
          function(x) .tCRT(.V2C(x)))

setMethod("toeplitz", signature(x = "sparseVector"),
          function(x, symmetric = TRUE, repr = c("C", "R", "T"),
                   giveCsparse, ...) {
              n <- length(x)
              if(n > .Machine$integer.max)
                  stop(gettextf("dimensions cannot exceed %s", "2^31-1"),
                       domain = NA)
              nn <- c(n, n)
              r <- spV2M(x[as.integer(abs(.col(nn) - .row(nn))) + 1L],
                         nrow = n, ncol = n, symmetric = symmetric,
                         check = FALSE)
              repr <- # keep in sync with sparseMatrix
                  if(missing(giveCsparse))
                      match.arg(repr)
                  else if(!missing(repr)) {
                      warning(gettextf("'%s' is deprecated; using '%s' instead",
                                       "giveCsparse", "repr"),
                              domain = NA)
                      match.arg(repr)
                  } else if(giveCsparse) {
                      "C"
                  } else {
                      warning(gettextf("'%s' is deprecated; setting %s=\"%s\"",
                                       "giveCsparse", "repr", "T"),
                              domain = NA)
                      "T"
                  }
              switch(repr, "C" = .M2C(r), "R" = .M2R(r), "T" = r)
          })
