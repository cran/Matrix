lu.Matrix <- function(x, norm.comp = c(one = TRUE, infinity = TRUE))
{
    .Call("R_LapackPP_lu", x, norm.comp)
}
