library(Matrix)

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
as(as(as(mT[1:2, 2:3], "dgCMatrix"), "lgCMatrix"),"ltCMatrix")


### Test all classes:  validObject(new( * )) should be fulfilled -----------

## need stoplist for now:
not.ok.classes <- paste(c("lgR", # only stub implementation
			  "lsR", # dito
			  "ltR", # dito

			  "ltT", # ltTMatrix_validate missing
			  "lsT", # lsTMatrix_validate missing

			  ""), "Matrix", sep='')
## From the rest, those that don't show :
no.show.classes <- paste(c("dgR", # only stub implementation
			   "dsR", # dito
			   "dtR"  # dito
			   ), "Matrix", sep='')

mM <- Matrix(1:4 >= 4, 2,2)
mm <- as(mM, "matrix")
for(cl in getClass("Matrix")@subclasses) {
    clNam <- cl@subClass
    cat(clNam)
    if(isVirtualClass(clNam)) {
	cat(" - is virtual\n")
    }
    else {
	cat("; new(..):")
	m <- new(clNam)
	if(any(clNam == not.ok.classes)) {
	    cat(" in 'stop list' - no validity\n")
	}
	else {
	    cat("valid: ", validObject(m), "\n")

	    ## The show() method implicitly tests
	    ##	as( <obj> , "matrix")
	    if(all(clNam != no.show.classes))
		show(m)
	    else cat("	-- no show() yet \n")

	    if(FALSE && is(m, "dMatrix")) { ## (FIXME) ?
		cat("as(dge*, <class>): ")
		m2 <- as(mM, clNam)
		cat("valid:", validObject(m2), "\n")

		if(TRUE) { ## (FIXME) ?
		    cat("as(matrix, <class>): ")
		    m3 <- as(mm, clNam)
		    cat("valid:", validObject(m3), "\n")
		}
	    }
	}
    }
}
