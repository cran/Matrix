sscCrosstab = 
    function(flist, permutations = TRUE)
{
    .Call("sscCrosstab", lapply(as(flist, "list"),
                                function(x) as(x, "factor")[drop = TRUE]),
          as(permutations, "logical")[1], PACKAGE = "Matrix")
}

