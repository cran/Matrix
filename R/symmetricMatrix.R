#### symmetricMatrix : virtual class


### autogenerate coercions
###  as(*,  "symmetricMatrix")
###  ~~~~~~~~~~~~~~~~~~~~~~~~~

## This should work at package install time when package:Matrix does not exist!
if(FALSE)
local({
    allCl <- getClasses("package:Matrix") ## << fails at install time!!!!
    clss <- allCl[sapply(allCl, extends, class2 = "Matrix")]
    virt <- sapply(clss, isVirtualClass)
    ## Now ensure coercions for all  non-virtual "Matrix" inheriting classes:
    for(cl in clss[!virt]) {
        cld <- getClassDef(cl)
        if(extends(cld, "symmetricMatrix"))
            cat("\tsymmetric:\t", cl,"\n")
        else if(extends(cld, "triangularMatrix"))
            cat("\ttriangular:\t", cl,"\n")
        else if(extends(cld, "diagonalMatrix"))
            cat("\tdiagonal:\t", cl,"\n")
        else {
            cat("do ",cl,"\n")
##             setAs(cl, "symmetricMatrix",
##                   function(from) as(from, ".s.Matrix"))
        }
    }## for
})
