#### Superclass Methods for all sparse nonzero-pattern matrices

setAs("CsparseMatrix", "nsparseMatrix",
      function(from) .Call(Csparse_to_nz_pattern, from,
			   is(from, "triangularMatrix")))
setAs("CsparseMatrix", "nMatrix",
      function(from) .Call(Csparse_to_nz_pattern, from,
			   is(from, "triangularMatrix")))

setAs("nsparseMatrix", "dsparseMatrix", function(from) as(from, "dMatrix"))


setMethod("is.na", signature(x = "nsparseMatrix"), is.na_nsp)

if(getRversion() >= "3.1.0")
setMethod("anyNA", signature(x = "nsparseMatrix"), function(x) FALSE)


setMethod("image", "nsparseMatrix", function(x, ...) image(as(x,"dMatrix")))
