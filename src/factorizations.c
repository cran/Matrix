#include "factorizations.h"

SEXP LU_validate(SEXP obj)
{
    return ScalarLogical(1);
}

SEXP Cholesky_validate(SEXP obj)
{
    return ScalarLogical(1);
}

SEXP SVD_validate(SEXP obj)
{
    return ScalarLogical(1);
}
