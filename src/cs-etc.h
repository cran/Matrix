#ifndef MATRIX_CS_ETC_H
#define MATRIX_CS_ETC_H

#include <Rinternals.h>
#include "SuiteSparse/CXSparse/Include/cs.h"

#define MCS_PATTERN 0
#define MCS_REAL    1
#define MCS_COMPLEX 2

#define MCS_XTYPE_GET(       ) Matrix_cs_xtype
#define MCS_XTYPE_SET(_VALUE_) Matrix_cs_xtype = _VALUE_

extern int Matrix_cs_xtype;

typedef struct Matrix_cs_sparse
{
	int nzmax;
	int m;
	int n;
	int *p;
	int *i;
	void *x; /* (double *) or (double _Complex *) */
	int nz;
	int xtype; /* Matrix-only */
} Matrix_cs;

typedef struct Matrix_cs_symbolic
{
	int *pinv;
	int *q;
	int *parent;
	int *cp;
	int *leftmost;
	int m2;
	double lnz;
	double unz;
} Matrix_css;

typedef struct Matrix_cs_numeric
{
	Matrix_cs *L; /* (cs_di *) or (cs_ci *) */
	Matrix_cs *U; /* (cs_di *) or (cs_ci *) */
	int *pinv;
	double *B;
} Matrix_csn;

typedef struct Matrix_cs_dmperm_results
{
	int *p;
	int *q;
	int *r;
	int *s;
	int nb;
	int rr[5];
	int cc[5];
} Matrix_csd;

Matrix_cs *M2CXS(SEXP, int);
SEXP CXS2M(Matrix_cs *, int, char);

/* Wrappers for the functions that we use at least once : */

Matrix_csd *Matrix_cs_dfree     (Matrix_csd *);
Matrix_csd *Matrix_cs_dmperm    (const Matrix_cs *, int);
int         Matrix_cs_dropzeros (Matrix_cs *);
void       *Matrix_cs_free      (void *);
int         Matrix_cs_happly    (const Matrix_cs *, int, double, void *);
int         Matrix_cs_ipvec     (const int *, const void *, void *, int);
int         Matrix_cs_lsolve    (const Matrix_cs *, void *);
Matrix_csn *Matrix_cs_lu        (const Matrix_cs *, const Matrix_css *, double);
int         Matrix_cs_lusol     (int, const Matrix_cs *, void *, double);
Matrix_csn *Matrix_cs_nfree     (Matrix_csn *);
Matrix_cs  *Matrix_cs_permute   (const Matrix_cs *, const int *, const int *, int);
int        *Matrix_cs_pinv      (const int *, int);
int         Matrix_cs_pvec      (const int *, const void *, void *, int);
Matrix_csn *Matrix_cs_qr        (const Matrix_cs *, const Matrix_css *);
int         Matrix_cs_qrsol     (int, const Matrix_cs *, void *);
Matrix_css *Matrix_cs_sfree     (Matrix_css *);
Matrix_cs  *Matrix_cs_spalloc   (int, int, int, int, int);
Matrix_cs  *Matrix_cs_speye     (int, int, int, int);
Matrix_cs  *Matrix_cs_spfree    (Matrix_cs *);
int         Matrix_cs_sprealloc (Matrix_cs *, int);
int         Matrix_cs_spsolve   (Matrix_cs *, const Matrix_cs *, int, int *, void *, const int *, int);
Matrix_css *Matrix_cs_sqr       (int, const Matrix_cs *, int);
Matrix_cs  *Matrix_cs_transpose (const Matrix_cs *, int);
int         Matrix_cs_usolve    (const Matrix_cs *, void *);

#endif /* MATRIX_CS_ETC_H */
