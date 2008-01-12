library(Matrix)
source(system.file("test-tools.R", package = "Matrix"))# identical3() etc
if(interactive()) options(error = recover)

#### Automatically display the class inheritance structure
#### possibly augmented with methods

allCl <- getClasses("package:Matrix")
cat("actual and virtual classes:\n")
tt <- table( isVirt <- sapply(allCl, isVirtualClass) )
names(tt) <- c('"actual"', "virtual")
tt
## The "actual" Matrix classes:
aCl <- allCl[!isVirt]
(aMcl <- aCl[grep("Matrix$", aCl)]) # length 48
aMc2 <-  aCl[sapply(aCl, extends, class2 = "Matrix")]
stopifnot(all( aMcl %in% aMc2 ))
aMc2[!(aMc2 %in% aMcl)] ## only 4 : p?Cholesky & p?BunchKaufman

## Really nice would be to construct an inheritance graph and display
## it.  Following things are computational variations on the theme..

## We use a version of  canCoerce()  that works with two *classes* instead of
## canCoerce <- function (object, Class)
classCanCoerce <- function (class1, class2)
{
    extends(class1, class2) ||
    !is.null(selectMethod("coerce", optional = TRUE,
			  signature    = c(from = class1, to = class2),
			  useInherited = c(from = TRUE,	  to = FALSE)))
}
.dq <- function(ch) paste0('"', ch, '"')
for(n in allCl) {
    if(isVirtualClass(n))
        cat("Virtual class", .dq(n),"\n")
    else {
        cat("\"Actual\" class", .dq(n),":\n")
        x <- new(n)
        for(m in allCl)
            if(classCanCoerce(n,m)) {
                ext <- extends(n, m)
                if(ext) {
                    cat(sprintf("   extends  %20s %20s \n", "", .dq(m)))
                } else {
                    cat(sprintf("   can coerce: %20s -> %20s: ", .dq(n), .dq(m)))
                    tt <- try(as(x, m), silent = TRUE)
                    if(inherits(tt, "try-error")) {
                        cat("\t *ERROR* !!\n")
                    } else {
                        cat("as() ok; validObject: ")
                        vo <- validObject(tt, test = TRUE)
                        cat(if(isTRUE(vo)) "ok" else paste("OOOOOOPS:", vo), "\n")
                    }
                }
            }
        cat("---\n")
    }
}

cat('Time elapsed: ', proc.time(),'\n') # for the above "part I"


if(!interactive()) { # don't want to see on source()

cat("All classes in the 'Matrix' package:\n")
for(cln in allCl) {
    cat("\n-----\n\nClass", dQuote(cln),":\n      ",
        paste(rep("~",nchar(cln)),collapse=''),"\n")
    ## A smarter version would use  getClass() instead of showClass(),
    ## build the "graph" and only then display.
    ##
    showClass(cln)
}

cat("\n\n")

## One could extend the `display' by using (something smarter than)
## are the "coerce" methods showing more than the 'Extends' output above?
cat("All (S4) methods in the 'Matrix' package:\n")
showMethods(where="package:Matrix")

} # end{non-interactive}

## 1-indexing instead of 0-indexing for direct "dgT" should give error:
ii <- as.integer(c(1,2,2))
jj <- as.integer(c(1,1,3))
assertError(new("dgTMatrix",  i=ii, j=jj,        x= 10*(1:3), Dim=2:3))
assertError(new("dgTMatrix",  i=ii, j=jj - 1:1,  x= 10*(1:3), Dim=2:3))
assertError(new("dgTMatrix",  i=ii - 1:1, j=jj,  x= 10*(1:3), Dim=2:3))
(mm <- new("dgTMatrix",  i=ii - 1:1, j=jj - 1:1, x= 10*(1:3), Dim=2:3))
validObject(mm)

### Sparse Logical:
m <- Matrix(c(0,0,2:0), 3,5)
mT <- as(mC <- as(m, "CsparseMatrix"), "TsparseMatrix")
stopifnot(identical(as(mT,"CsparseMatrix"), mC))
(mC. <- as(mT[1:2, 2:3], "CsparseMatrix"))
(mlC <- as(mC. , "lMatrix"))
as(mlC,"ltCMatrix")



### Test all classes:  validObject(new( * )) should be fulfilled -----------

## need stoplist for now:
Rcl.struc <- c("gR", "sR", "tR")
(dR.classes <- paste0(paste0("d", Rcl.struc[Rcl.struc != "gR"]),   "Matrix"))
(.R.classes <- paste0(sort(outer(c("l", "n"), Rcl.struc, paste0)), "Matrix"))
                                        # have only stub implementation

## not.ok..: are left out almost completely
not.ok.classes <- NULL  ## was  .R.classes
## From the rest, those that don't show {have no coerce to "dge":}
no.show.classes <- NULL ## was  dR.classes
##
Mat.MatFact <- c("Cholesky", "pCholesky",
                 "BunchKaufman", "pBunchKaufman")##, "LDL"
no.t.etc <- c(no.show.classes, .R.classes, dR.classes, Mat.MatFact)
no.t.classes <- c(no.t.etc)     # no t() available
no.norm.classes <- no.t.classes
not.Ops      <- no.show.classes # "Ops", e.g. "+" fails
not.coerce0  <- no.show.classes # not coercable to   "matrix" & "dgeMatrix"
not.coerce1  <- no.t.etc        # not coercable from "dgeMatrix"
not.coerce2  <- no.t.etc        # not coercable from "matrix"

tstMatrixClass <-
    function(cl, mM = Matrix(c(2,1,1,2) + 0, 2,2,
                 dimnames=rep( list(c("A","B")), 2)), # dimnames: *symmetric*
             mm = as(mM, "matrix"), recursive = TRUE, offset = 0)
{
    ## Purpose: Test 'Matrix' class {and do this for all of them}
    ## ----------------------------------------------------------------------
    ## Arguments: cl: class object of a class that extends "Matrix"
    ##            mM: a "Matrix"-matrix which will be coerced to class 'cl'
    ##            mm: a S3-matrix       which will be coerced to class 'cl'
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler

    ## This is sfsmisc::bl.string():
    bl.string <- function (no) paste(rep(" ", no), collapse = "")

    ## Compute a few things only once :
    mM <- as(mM, "dgeMatrix")
    trm <- mm; trm[lower.tri(mm)] <- 0
    summList <- lapply(getGroupMembers("Summary"), get,
                       envir = asNamespace("Matrix"))
    if(recursive)
        cList <- character(0)

    extraValid <- function(m, cl = class(m)) {
        sN <- slotNames(cl)
        sN <- sN[sN != "factors"]
	for(nm in sN)
	    if(!is.null(a <- attributes(slot(m, nm))))
		stop(sprintf("slot '%s' with %d attributes, named: ",
			     nm, length(a)), paste(names(a), collapse=", "))
        invisible(TRUE)
    }


    ## This is the recursive function
    dotestMat <- function(cl, offset)
    {
        cat. <- function(...) cat(bl.string(offset), ...)

        clNam <- cl@subClass
        cat("\n")
        cat.(clNam)
        ##---------
        if(isVirtualClass(clNam)) {
            cat(" - is virtual\n")
            if(recursive) {
                cat.("----- begin{class :", clNam, "}----new subclasses----\n")
                for(ccl in getClass(clNam)@subclasses) {
                    cclN <- ccl@subClass
                    if(cclN %in% cList)
                        cat.(cclN,": see above\n")
                    else {
                        cList <<- c(cList, cclN)
                        dotestMat(ccl, offset = offset + 3)
                    }
                }
                cat.("----- end{class :", clNam, "}---------------------\n")
            }
        } else { ## --- actual class ---
            genC <- extends(clNam, "generalMatrix")
            symC <- extends(clNam, "symmetricMatrix")
            triC <- extends(clNam, "triangularMatrix")
            diaC <- extends(clNam, "diagonalMatrix")
            if(!(genC || symC || triC || diaC))
                stop("does not extend one of 'general', 'symmetric', 'triangular', or 'diagonal'")
            sparseC <- extends(clNam, "sparseMatrix")
            denseC  <- extends(clNam, "denseMatrix")
            if(!(sparseC || denseC))
                stop("does not extend either 'sparse' or 'dense'")
	    cat("; new(..): ")
	    m <- new(clNam) ; cat("ok; ")
            is_p <- extends(clNam, "pMatrix")
            is_cor <- (clNam == "corMatrix") # has diagonal divided out
	    if(canCoerce(mm, clNam)) { ## replace 'm' by `non-empty' version
		cat("canCoerce() ")
		m0 <- {
		    if(triC) trm
		    else if(is_p)
			mm == 1 # logical *and* "true" permutation
		    else mm
		}
		if(extends(clNam, "lMatrix") ||
		   extends(clNam, "nMatrix"))
		    storage.mode(m0) <- "logical"
		else if(extends(clNam, "zMatrix"))
		    storage.mode(m0) <- "complex"
		validObject(m) ## validity of trivial 'm' before replacing
		m <- as(m0, clNam)
		if(is_cor)
                    m0 <- cov2cor(m0)
	    } else {
                m0 <- vector(Matrix:::.type.kind[Matrix:::.M.kindC(clNam)])
                dim(m0) <- c(0L,0L)
            }
	    ## m0 is the 'matrix' version of our 'Matrix' m
	    m. <- if(is_p) as.integer(m0) else m0
            EQ <- if(is_cor) all.equal else identical
	    stopifnot(EQ(m0[FALSE], m[FALSE]),
		      EQ(m.[TRUE],  m[TRUE]),
		      if(length(m) >= 2) EQ(m.[2:1], m[2:1]) else TRUE)

	    if(any(clNam == not.ok.classes)) {
		cat("in 'stop list' - no validity\n")
	    } else {
		cat("valid: ", validObject(m), extraValid(m, clNam))

		## This can only work as long as 'm' has no NAs :
                ## not yet -- have version in not.Ops below
## once we have is.na():
## 		stopifnot(all(m == m | is.na(m))) ## check all() and "==" [Compare]
## 		if(any(m != m && !is.na(m)))
		stopifnot(all(m == m)) ## check all() and "==" [Compare]
		if(any(m != m))
		    stop(" any (m != m) should not be true")

                if(clNam %in% no.t.classes) {
                    cat(" in t()-'stop list'\n")
                } else {
                    cat("; t(t(m)) ==?== m :")
                    stopifnot(Qidentical(m, t(t(m))))
                    cat(" ok\n")
                }
                if(all(clNam != no.show.classes))
                    show(m)
                  ## improve: cat.(  captureOutput(show(m) ) )
                else cat.("	-- no show() yet \n")

                if(clNam %in% no.norm.classes) {
                    cat.(" in norm()-'stop list'\n")
                } else {
                    cat.(sprintf(" norm(m [%d x %d]) :", nrow(m), ncol(m)))
                    for(typ in c("1","I","F","M")) {
                        cat('', typ, '')
                        stopifnot(all.equal(norm(m,typ), norm(m0,typ)))
                    }
                    cat(" ok\n")
                }
		cat.(" Summary: ")
		for(f in summList) {
		    #cat(f@generic, "() ", sep='')
		    stopifnot(if(clNam == "corMatrix") all.equal(f(m), f(m0))
			      else f(m) == f(m0))
		}
		cat(" ok\n")

                if(all(clNam != not.coerce0)) {## coerce to 'matrix'
                    m.m <- as(m, "matrix")
                    ## and test 'dim()' as well:
		    stopifnot(identical(dim(m.m), dim(m)),
			      if(clNam %in% no.t.classes)
			      TRUE else identical(diag(m), diag(t(m))),
			      ## TODO: also === diag(band(m,0,0))
			      diag(m) == diag(m.m),
			      nnzero(m) == sum(m.m != 0))
                }
                else stopifnot(length(dim(m)) == 2)

### FIXME: organize differently :
### 1) produce 'mM'  and 'mm' for the other cases,
### 2) use identical code for all cases

                ## use non-square matrix when "allowed":

		## "!" should work (via as(*, "l...")) :
                m11 <- as(as(!!m,"CsparseMatrix"), "lMatrix")
                m12 <- as(as(  m, "lMatrix"),"CsparseMatrix")
                if(!identical(m11, m12))
                    stopifnot(identical(as(m11, "generalMatrix"),
                                        as(m12, "generalMatrix")))

		if(is(m, "dMatrix") && all(clNam != not.Ops)) {
                    ## makes sense with non-trivial m (!)
		    cat("2*m =?= m+m: ")
		    if(identical(2*m, m+m)) cat("identical\n")
		    else {
			stopifnot(as(2*m,"matrix") == as(m+m, "matrix"))
			cat("ok\n")
		    }
                    ## m == m etc, now for all, see above
		    cat("m >= m for all: ")
		    stopifnot(all(m >= m)); cat("ok\n")
		    cat("m < m for none: ")
		    stopifnot(all(! m < m)); cat("ok\n")
		}

                if(is(m, "dMatrix") && is(m, "compMatrix")) {
                    if(any(clNam == not.coerce1))
                        cat.("not coercable_1\n")
		    else {
			cat.("as(dge*, <(super)class>): ")
			if(canCoerce(mM, clNam))
			    m2 <- as(mM, clNam)
			else { ## find superclass to which to coerce
			    if(extends(clNam, "sparseMatrix")) {
				if(is.na(newcl <- Matrix:::.sp.class(clNam)))
				    stop("internal failure from .sp.class()")
				m2 <- as(mM, newcl)
			    } else { ## ddense & (general or symmetric)
				stop("don't know what to coerce <dge> to - error test-logic")
			    }
			}
			cat("valid:", validObject(m2), "\n")
			if(!is_cor) ## as.vector()
			    stopifnot(as.vector(m2) == as.vector(mM))
			cat.("[cr]bind2():"); mm2 <- cbind2(m2,m2)
			stopifnot(dim(rbind2(m2,m2)) == 2:1 * dim(mM)); cat(" ok")
			if(genC && class(mm2) == clNam)## non-square matrix when "allowed"
			    m2 <- mm2
			dd <- diag(m2)
			cat("; `diag<-` ")
			diag(m2) <- 10*dd
			stopifnot(is_cor || identical(dd, diag(mM)),
				  identical(10*dd, diag(m2))); cat("ok ")
                    }
                    if(all(clNam != not.coerce2)) {
                        cat.("as(matrix, <class>): ")
                        m3 <- as(mm, clNam)
                        cat("valid:", validObject(m3), "\n")
                    }
                }
                else { ## not numeric composite: logical / triangular/diagonal ..
                    if(any(clNam == not.coerce1))
                        cat.("not coercable_1\n")
                    else {
                        ## FIXME: also add tests for these
                        if(is(m, "lMatrix")) { ## should fulfill even with NA:
			    stopifnot(identical(m, m & TRUE),
				      identical(m, FALSE | m),
				      all(m | !m), !any(!m & m))
                        }
                        else if(is(m, "triangularMatrix")) {
                            mm. <- mm
                            i0 <- if(m@uplo == "L")
                                upper.tri(mm.) else lower.tri(mm.)
                            mm.[i0] <- 0
			    cat.("as(triang.matrix, <class>): ")
			    m3 <- as(mm., clNam)
                            cat("valid:", validObject(m3), "\n")
                        }
                        else { ## diagonal (only one)?
                            ## TODO
                        }
                    }
                }

                ## if(is(m, "denseMatrix")) {
                ##     ## .........
                ##     cat.("as dsparse* ")
                ##     msp <- as(m, "dsparseMatrix")
                ##     cat.("; valid coercion: ", validObject(msp), "\n")
                ## } else if(is(m, "sparseMatrix")) {

                ## } else cat.("-- not dense nor sparse -- should not happen(!?)\n")

                if(is(m, "dsparseMatrix")) {
                    if(any(clNam == not.coerce1))
                        cat.("not coercable_1\n")
                    else {
			## make sure we can coerce to dgT* -- needed, e.g. for "image"
			## change: use Tsparse instead of dgT, unless it *is* Tsparse:
			isT <- is(m, "TsparseMatrix")
			prefix <- if(isT) "dgT" else "Tsparse"
			Tcl <- paste(prefix, "Matrix", sep='')
			cat.(sprintf("as %s* ", prefix))
			mgT <- as(m, Tcl)
			cat(sprintf("; valid %s* coercion: %s\n",
				    prefix, validObject(mgT)))
                    }
                }
            }
        }
    } # end{dotestMat}

    for(scl in getClass(cl)@subclasses)
        dotestMat(scl, offset + 1)
}

tstMatrixClass("Matrix")
if(FALSE)## or just a sub class
tstMatrixClass("triangularMatrix")


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

if(!interactive()) warnings()
