#ifndef MATRIX_MUTILS_H
#define MATRIX_MUTILS_H

#include <Rdefines.h>
#include <Rconfig.h>
#include "cblas.h"

/* short forms of some enum constants from cblas.h */
#define RMJ CblasRowMajor
#define CMJ CblasColMajor
#define NTR CblasNoTrans
#define TRN CblasTrans
#define CTR CblasConjTrans
#define UPP CblasUpper
#define LOW CblasLower
#define NUN CblasNonUnit
#define UNT CblasUnit
#define LFT CblasLeft
#define RGT CblasRight

char norm_type(char *typstr);
char rcond_type(char *typstr);
double get_double_by_name(SEXP obj, char *nm);
SEXP set_double_by_name(SEXP obj, double val, char *nm);
SEXP as_det_obj(double val, int log, int sign);
SEXP get_factorization(SEXP obj, char *nm);
SEXP set_factorization(SEXP obj, SEXP val, char *nm);
SEXP cscMatrix_set_Dim(SEXP x, int nrow);
int csc_unsorted_columns(int ncol, const int p[], const int i[]);
void csc_sort_columns(int ncol, const int p[], int i[], double x[]);
SEXP triple_as_SEXP(int nrow, int ncol, int nz,
		    const int Ti [], const int Tj [], const double Tx [],
		    char *Rclass);
SEXP csc_check_column_sorting(SEXP A);
void csc_components_transpose(int m, int n, int nnz,
			      const int xp[], const int xi[],
			      const double xx[],
			      int ap[], int ai[], double ax[]);
void ssc_symbolic_permute(int n, int upper, const int perm[],
			  int Ap[], int Ai[]);
double *nlme_symmetrize(double *a, const int nc);
void nlme_check_Lapack_error(int info, const char *laName);
SEXP nlme_replaceSlot(SEXP obj, SEXP names, SEXP value);
SEXP nlme_weight_matrix_list(SEXP MLin, SEXP wts, SEXP adjst, SEXP MLout);

				/* stored pointers to symbols */
				/* initialized in R_init_Matrix */
extern
#include "Syms.h"

/* zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}


#endif

