//   R : A Computer Language for Statistical Data Analysis
//   Copyright (C) 2000  the R Development Core Team

//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.

//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "orthonormal.h"

SEXP LaColOrthogonalMatDouble::asSEXP() const
{
    if (size(0) < 1 || size(1) < 1) return R_NilValue;

    SEXP val = PROTECT(allocMatrix(REALSXP, size(0), size(1)));
    F77_CALL(dlacpy)('A', size(0), size(1), &(*this)(0,0), gdim(0),
		     REAL(val), size(0));
    SEXP classes = PROTECT(allocVector(STRSXP, 2));
    STRING(classes)[0] = mkChar("ColOrthogonal");
    STRING(classes)[1] = mkChar("Matrix");
    setAttrib(val, R_ClassSymbol, classes);
    UNPROTECT(2);
    return val;
}

SEXP LaRowOrthogonalMatDouble::asSEXP() const
{
    if (size(0) < 1 || size(1) < 1) return R_NilValue;

    SEXP val = PROTECT(allocMatrix(REALSXP, size(0), size(1)));
    F77_CALL(dlacpy)('A', size(0), size(1), &(*this)(0,0), gdim(0),
		     REAL(val), size(0));
    SEXP classes = PROTECT(allocVector(STRSXP, 2));
    STRING(classes)[0] = mkChar("RowOrthogonal");
    STRING(classes)[1] = mkChar("Matrix");
    setAttrib(val, R_ClassSymbol, classes);
    UNPROTECT(2);
    return val;
}

SEXP LaColOrthonormalMatDouble::asSEXP() const
{
    if (size(0) < 1 || size(1) < 1) return R_NilValue;

    SEXP val = PROTECT(allocMatrix(REALSXP, size(0), size(1)));
    F77_CALL(dlacpy)('A', size(0), size(1), &(*this)(0,0), gdim(0),
		     REAL(val), size(0));
    SEXP classes = PROTECT(allocVector(STRSXP, 3));
    STRING(classes)[0] = mkChar("ColOrthonormal");
    STRING(classes)[1] = mkChar("ColOrthogonal");
    STRING(classes)[2] = mkChar("Matrix");
    setAttrib(val, R_ClassSymbol, classes);
    UNPROTECT(2);
    return val;
}

SEXP LaRowOrthonormalMatDouble::asSEXP() const
{
    if (size(0) < 1 || size(1) < 1) return R_NilValue;

    SEXP val = PROTECT(allocMatrix(REALSXP, size(0), size(1)));
    F77_CALL(dlacpy)('A', size(0), size(1), &(*this)(0,0), gdim(0),
		     REAL(val), size(0));
    SEXP classes = PROTECT(allocVector(STRSXP, 3));
    STRING(classes)[0] = mkChar("RowOrthonormal");
    STRING(classes)[1] = mkChar("RowOrthogonal");
    STRING(classes)[2] = mkChar("Matrix");
    setAttrib(val, R_ClassSymbol, classes);
    UNPROTECT(2);
    return val;
}

SEXP LaOrthogonalMatDouble::asSEXP() const
{
    int m = LaGenMatDouble::size(0), n = LaGenMatDouble::size(1),
	lda = LaGenMatDouble::gdim(0);

    if (m < 1 || n < 1) return R_NilValue;

    SEXP val = PROTECT(allocMatrix(REALSXP, m, n));
    F77_CALL(dlacpy)('A', m, n, &(*this)(0,0), lda, REAL(val), m);
    SEXP classes = PROTECT(allocVector(STRSXP, 6));
    STRING(classes)[0] = mkChar("Orthogonal");
    STRING(classes)[1] = mkChar("RowOrthonormal");
    STRING(classes)[2] = mkChar("RowOrthogonal");
    STRING(classes)[3] = mkChar("ColOrthonormal");
    STRING(classes)[4] = mkChar("ColOrthogonal");
    STRING(classes)[5] = mkChar("Matrix");
    setAttrib(val, R_ClassSymbol, classes);
    UNPROTECT(2);
    return val;
}
