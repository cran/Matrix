#ifndef MATRIX_FLAME_H
#define MATRIX_FLAME_H

#include "Rdefines.h"
#include "R_ext/Lapack.h"
#include "FLAME/FLAME.h"

FLA_Obj *R_to_FLA_copy(SEXP Ain);

FLA_Obj *R_to_FLA_inPlace(SEXP Ain);

SEXP R_FLA_Init();
SEXP R_FLA_Finalize();

#define RFLAME_CHOL_NB 104
#define RFLAME_QR_NB    96

SEXP lsq_Chol_flame(SEXP Xin, SEXP yin);

#endif
