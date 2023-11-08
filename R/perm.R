invertPerm <- function(p, off = 1L, ioff = 1L)
    .Call(R_invertPerm, as.integer(p), as.integer(off), as.integer(ioff))

signPerm <- function(p, off = 1L)
    .Call(R_signPerm, as.integer(p), as.integer(off))

isPerm <- function(p, off = 1L)
    .Call(R_isPerm, as.integer(p), as.integer(off))

asPerm <- function(pivot, off = 1L, ioff = 1L, n = length(pivot))
    .Call(R_asPerm, as.integer(pivot), as.integer(off), as.integer(ioff),
          as.integer(n))

invPerm <- function(p, zero.p = FALSE, zero.res = FALSE)
    invertPerm(p, if(zero.p) 0L else 1L, if(zero.res) 0L else 1L)
