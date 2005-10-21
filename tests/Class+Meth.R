library(Matrix)
source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

#### Automatically display the class inheritance structure
#### possibly augmented with methods

allCl <- getClasses("package:Matrix")

## Really nice would be to construct an inheritance graph and display
## it.  The following is just a cheap first step.

cat("All classes in the 'Matrix' package:\n")
for(cln in allCl) {
    cat("\n-----\n\nClass", dQuote(cln),":\n      ",
        paste(rep("~",nchar(cln)),collapse=''),"\n")
    ## A smarter version would use  getClass() instead of showClass(),
    ## build the "graph" and only then display.

    showClass(cln)
}

cat("\n\n")

## One could extend the `display' by using (something smarter than)
## are the "coerce" methods showing more than the 'Extends' output above?
cat("All (S4) methods in the 'Matrix' package:\n")
showMethods(where="package:Matrix")


### Sparse Logical:
m <- Matrix(c(0,0,2:0), 3,5)
mT <- as(mC <- as(m, "dgCMatrix"), "dgTMatrix")
stopifnot(identical(as(mT,"dgCMatrix"), mC))
(mlC <- as(as(mT[1:2, 2:3], "dgCMatrix"), "lgCMatrix"))

if(FALSE) ## ltC no longer extends lgC -- want coercion possibility FIXME
as(mlC,"ltCMatrix")


### Test all classes:  validObject(new( * )) should be fulfilled -----------

## need stoplist for now:
not.ok.classes <- paste(c("lgR", # only stub implementation
			  "lsR", # dito
			  "ltR", # dito

			  "ltT", # ltTMatrix_validate missing; as(*,"matrix")
			  "lsT", # lsTMatrix_validate  "	"

			  ""), "Matrix", sep='')
## From the rest, those that don't show :
no.show.classes <- paste(c("dgR", # only stub implementation
			   "dsR", # dito
			   "dtR", #  "
			   ), "Matrix", sep='')

no.t.classes <- no.show.classes # for the moment

mM <- Matrix(1:4 >= 4, 2,2)
mm <- as(mM, "matrix")
for(cl in getClass("Matrix")@subclasses) {
    clNam <- cl@subClass
    cat(clNam)
    if(isVirtualClass(clNam)) {
	cat(" - is virtual\n")
    } else {
	cat("; new(..):")
	m <- new(clNam)

	if(any(clNam == not.ok.classes)) {
	    cat(" in 'stop list' - no validity\n")
	} else {
	    cat("valid: ", validObject(m))

            if(any(clNam == no.t.classes)) {
                cat(" in t()-'stop list'\n")
            } else {
                cat("; t(t(m)) ==?== m :")
                stopifnot(Qidentical(m, t(t(m))))
                cat(" ok\n")
            }

	    ## The show() method implicitly tests as( <obj> , "matrix"):
	    if(all(clNam != no.show.classes))
		show(m)
	    else cat("	-- no show() yet \n")

	    if(is(m, "dMatrix")) {
                if(FALSE) { ## (FIXME) ?
                    cat("as(dge*, <class>): ")
                    m2 <- as(mM, clNam)
                    cat("valid:", validObject(m2), "\n")
                }
		if(FALSE) { ## FIXME or use another stoplist; fails for 'dsy'
		    cat("as(matrix, <class>): ")
		    m3 <- as(mm, clNam)
		    cat("valid:", validObject(m3), "\n")
		}
	    }

##             if(is(m, "denseMatrix")) {
##                 ## .........
##                 cat("as dsparse* ")
##                 msp <- as(m, "dsparseMatrix")
##                 cat("; valid coercion: ", validObject(msp), "\n")
##             } else if(is(m, "sparseMatrix")) {

##             } else cat("-- not dense nor sparse -- should not happen(!?)\n")

            if(is(m, "dsparseMatrix")) {
                ## make sure that we can coerce to  dgT* -- is needed, e.g. for "image"
                cat("as dgT* ")
                mgT <- as(m, "dgTMatrix")
                cat("; valid dgT* coercion: ", validObject(mgT), "\n")
            }
	}
    }
}
