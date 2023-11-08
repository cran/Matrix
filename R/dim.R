validDim <- function(dim)
    .Call(R_Dim_validate, dim)

validDimGetsValue <- function(value, mn) {
    if(mode(value) != "numeric")
        gettextf("assigned dimensions are not of type \"%s\" or \"%s\"",
                 "integer", "double")
    else if(length(value) != 2L)
        gettextf("assigned dimensions do not have length %d", 2L)
    else if(anyNA(value))
        gettext("assigned dimensions are NA")
    else if(any(value < 0L))
        gettext("assigned dimensions are negative")
    else if(is.double(value) && any(trunc(value) > .Machine$integer.max))
        gettextf("assigned dimensions exceed %s", "2^31-1")
    else if((p <- prod(value)) != mn)
        gettextf("assigned dimensions [product %.0f] do not match object length [%.0f]",
                 p, as.double(mn))
    else TRUE
}

validDN <- function(dn, dim)
    .Call(R_DimNames_validate, dn, dim)

fixupDN <- function(dn)
    .Call(R_DimNames_fixup, dn)

fixupDN.if.valid <- function(dn, dim) {
    if(is.character(s <- validDim(dim)) || is.character(s <- validDN(dn, dim)))
        stop(s, domain = NA)
    fixupDN(dn)
}

symDN <- function(dn)
    .Call(R_symDN, dn)

symmetrizeDN <- function(x) {
    if(isS4(x)) # assuming is(x, "Matrix")
        `dimnames<-`(x, symDN(x@Dimnames))
    else if(!is.null(dn <- dimnames(x))) # assuming list of length 2
        `dimnames<-`(x, symDN(dn))
    else x
}

isSymmetricDN <- function(dn)
    .Call(R_DimNames_is_symmetric, dn)

is.null.DN <- function(dn) {
    if(is.null(dn))
        return(TRUE)
    if(!is.null(names(dn)))
        names(dn) <- NULL
    ch0 <- character(0L)
    identical(dn, list(NULL, NULL)) ||
    identical(dn, list( ch0, NULL)) ||
    identical(dn, list(NULL,  ch0)) ||
    identical(dn, list( ch0,  ch0))
}


## METHODS FOR GENERIC: dim
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dim", signature(x = "Matrix"),
          function(x) x@Dim)

setMethod("dim", signature(x = "MatrixFactorization"),
          function(x) x@Dim)


## METHODS FOR GENERIC: dim<-
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dim<-", signature(x = "denseMatrix"),
          function(x, value) {
              if(is.character(s <- validDimGetsValue(value, prod(d <- x@Dim))))
                 stop(s, domain = NA)
              value <- as.integer(value)
              if(all(value == d))
                  return(x)
              r <- .M2gen(x)
              r@Dim <- value
              if(length(r@factors))
                  r@factors <- list()
              r
          })

setMethod("dim<-", signature(x = "sparseMatrix"),
          function(x, value) {
              if(is.character(s <- validDimGetsValue(value, prod(d <- x@Dim))))
                 stop(s, domain = NA)
              value <- as.integer(value)
              if(all(value == d))
                  return(x)
              r <- spV2M(.M2V(x), nrow = value[1L], ncol = value[2L])
              switch(.M.repr(x), "C" = .M2C(r), "R" = .M2R(r), r)
          })

setMethod("dim<-", signature(x = "sparseVector"),
          function(x, value) {
              if(is.character(s <- validDimGetsValue(value, length(x))))
                 stop(s, domain = NA)
              value <- as.integer(value)
              spV2M(x, nrow = value[1L], ncol = value[2L])
          })


## METHODS FOR GENERIC: length
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("length", "Matrix",
          function(x)
              if((r <- prod(x@Dim)) > .Machine$integer.max)
                  r
              else as.integer(r))

setMethod("length", "MatrixFactorization",
          function(x)
              if((r <- prod(x@Dim)) > .Machine$integer.max)
                  r
              else as.integer(r))

setMethod("length", "sparseVector",
          function(x)
              if(is.integer(r <- x@length) || r > .Machine$integer.max)
                  r
              else as.integer(r))


## METHODS FOR GENERIC: dimnames
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dimnames", signature(x = "Matrix"),
          function(x) x@Dimnames)

setMethod("dimnames", signature(x = "symmetricMatrix"),
          function(x) symDN(x@Dimnames))

setMethod("dimnames", signature(x = "MatrixFactorization"),
          function(x) x@Dimnames)


## METHODS FOR GENERIC: dimnames<-
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dimnames<-", signature(x = "Matrix", value = "NULL"),
          function(x, value) {
              x@Dimnames <- list(NULL, NULL)
              x
          })

setMethod("dimnames<-", signature(x = "compMatrix", value = "NULL"),
          function(x, value) {
              if(length(x@factors))
                  x@factors <- list()
              x@Dimnames <- list(NULL, NULL)
              x
          })

setMethod("dimnames<-", signature(x = "MatrixFactorization", value = "NULL"),
          function(x, value) {
              x@Dimnames <- list(NULL, NULL)
              x
          })

setMethod("dimnames<-", signature(x = "Matrix", value = "list"),
          function(x, value) {
              x@Dimnames <- fixupDN.if.valid(value, x@Dim)
              x
          })

setMethod("dimnames<-", signature(x = "compMatrix", value = "list"),
          function(x, value) {
              if(length(x@factors))
                  x@factors <- list()
              x@Dimnames <- fixupDN.if.valid(value, x@Dim)
              x
          })

setMethod("dimnames<-", signature(x = "MatrixFactorization", value = "list"),
          function(x, value) {
              x@Dimnames <- fixupDN.if.valid(value, x@Dim)
              x
          })


## METHODS FOR GENERIC: unname
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("unname", signature(obj = "Matrix"),
          function(obj, force = FALSE) {
              obj@Dimnames <- list(NULL, NULL)
              obj
          })

setMethod("unname", signature(obj = "MatrixFactorization"),
          function(obj, force = FALSE) {
              obj@Dimnames <- list(NULL, NULL)
              obj
          })


## METHODS FOR GENERIC: drop
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("drop", signature(x = "Matrix"),
          function(x) if(any(x@Dim == 1L)) drop(.M2m(x)) else x)
