#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <stddef.h> /* size_t */
#include <Rinternals.h>
#include <Rversion.h>

#if R_VERSION < R_Version(4, 5, 0)
int ANY_ATTRIB(SEXP);
void CLEAR_ATTRIB(SEXP);
#endif

void *Matrix_memset(void *, int, R_xlen_t, size_t);
void *Matrix_memcpy(void *, const void *, R_xlen_t, size_t);
char *Matrix_sprintf(const char *, ...);

int equal_character_vectors(SEXP, SEXP, int);

void conjugate(SEXP);
void zeroRe(SEXP);
void zeroIm(SEXP);
void naToOne(SEXP);

#endif /* MATRIX_UTILS_H */
