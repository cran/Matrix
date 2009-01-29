## SPQR currently experimental, but possibly THE method with potential
## ----
## From spqr_user_guide.pdf :
## ordering : string describing what ordering method to use. Let [m2 n2]=size(S)
## where S is obtained by removing singletons from A. The singleton permutation places
## A*P in the form [A11 A12 ; 0 S] where A11 is upper triangular with diagonal entries
## all greater than tol.
## The default is to use COLAMD if m2<=2*n2; otherwise try AMD. Let f be the flops for
## chol((S*P)’*(S*P)) with the ordering P found by AMD. Then if f/nnz(R) >= 500
## and nnz(R)/nnz(S) >= 5 then try METIS, and take the best ordering found (AMD
## or METIS); otherwise use AMD without trying METIS. If METIS is not installed then
## the default ordering is to use COLAMD if m2<=2*n2 and to use AMD otherwise.
## The available orderings are:

## ’default’: the default ordering.
## ’amd’: use amd(S’*S).
## ’colamd’: use colamd(S).
## ’metis’: use metis(S’*S), only if METIS is installed.
## ’best’: try all three (AMD, COLAMD, METIS) and take the best.
## ’bestamd’: try AMD and COLAMD and take the best.
## ’fixed’: use P=I; this is the only option if P is not present in the output.
## ’natural’: singleton removal only.

#define SPQR_ORDERING_FIXED 0
#define SPQR_ORDERING_NATURAL 1
#define SPQR_ORDERING_COLAMD 2
#define SPQR_ORDERING_GIVEN 3       /* only used for C/C++ interface */
#define SPQR_ORDERING_CHOLMOD 4     /* CHOLMOD best-effort (COLAMD, METIS,...)*/
#define SPQR_ORDERING_AMD 5         /* AMD(A'*A) */
#define SPQR_ORDERING_METIS 6       /* metis(A'*A) */
#define SPQR_ORDERING_DEFAULT 7     /* SuiteSparseQR default ordering */
#define SPQR_ORDERING_BEST 8        /* try COLAMD, AMD, and METIS; pick best */
#define SPQR_ORDERING_BESTAMD 9     /* try COLAMD and AMD; pick best */

## here (instead of ./AllClass.R for now):
setClass("SPQR", contains = "MatrixFactorization",
	 representation(Q = "dgCMatrix", R = "dgCMatrix",
			p = "integer", rank = "integer"),
	 validity = function(object) {
	     dQ <- dim(object@Q)
	     dR <- dim(object@R)
	     if(length(rnk <- object@rank) != 1 || rnk < 0 || rnk > dR[2])
		 "invalid rank"
	     else if(dQ[2] != dR[1])
		 "ncol(Q) != nrow(R)"
	     else if(length(object@p) != dR[2])
		 "length(p) != ncol(R)"
	     else TRUE
	 })


## Using spqr() is as Tim Davis has named his Matlab interface;
## but the function arguments will differ anyway
setGeneric("spqr", function(x, ...) standardGeneric("spqr"))

setMethod("spqr", signature(x = "dgCMatrix"),
	  function(x, econ = 0,
                   ordering = c("default", "fixed","natural","amd",
                                "colamd", "metis", "best", "bestamd"),
                   tol = -2)
      {
          ordering <- as.integer(
                                 c("fixed" = 0,
                                   "natural" = 1,
                                   "colamd" = 2,
                                   "given" = 3,
                                   "cholmod" = 4,
                                   "amd" = 5    ,
                                   "metis" = 6  ,
                                   "default" = 7,
                                   "best" = 8   ,
                                   "bestamd" = 9)[match.arg(ordering)])
          .Call(dgCMatrix_SPQR, x, ordering, econ, tol)
      })

setMethod("spqr", signature(x = "TsparseMatrix"),
	  function(x, ...) spqr(as(x, "CsparseMatrix"), ...))
setMethod("spqr", signature(x = "CsparseMatrix"),
	  function(x, ...) spqr(as(as(x, "generalMatrix"), "dgCMatrix"), ...))


## an alternative would be (a version of) the following interface:
## (--> ./dgCMatrix.R )
if(FALSE)
setMethod("qr", signature(x = "dgCMatrix"),
	  function(x, tol = 1e-07, LAPACK = FALSE, SPQR = FALSE)
          ## SPQR currently experimental, but possibly THE method with potential
          if(SPQR).Call(dgCMatrix_SPQR, x, ncol(x)) else
	  .Call(dgCMatrix_QR, x, TRUE))
