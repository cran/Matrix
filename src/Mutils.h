#ifndef MATRIX_MUTILS_H
#define MATRIX_MUTILS_H

#include <Rdefines.h>

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
void triplet_to_col(int nrow, int ncol, int nz,
		    const int Ti [], const int Tj [], const double Tx [],
		    int Ap [], int Ai [], double Ax []);
void ssc_symbolic_permute(int n, int upper, const int perm[],
			  int Ap[], int Ai[]);
double *nlme_symmetrize(double *a, const int nc);
void nlme_check_Lapack_error(int info, const char *laName);
SEXP nlme_replaceSlot(SEXP obj, SEXP names, SEXP value);
SEXP nlme_weight_matrix_list(SEXP MLin, SEXP wts, SEXP adjst, SEXP MLout);

				/* stored pointers to symbols */
				/* initialized in Matrix_init */
extern SEXP
    Matrix_DSym,
    Matrix_DIsqrtSym,
    Matrix_DimSym,
    Matrix_GpSym,
    Matrix_LiSym,
    Matrix_LpSym,
    Matrix_LxSym,
    Matrix_OmegaSym,
    Matrix_ParentSym,
    Matrix_RXXSym,
    Matrix_RZXSym,
    Matrix_XtXSym,
    Matrix_ZtXSym,
    Matrix_ZZxSym,
    Matrix_bVarSym,
    Matrix_cnamesSym,
    Matrix_devianceSym,
    Matrix_devCompSym,
    Matrix_diagSym,
    Matrix_iSym,
    Matrix_ipermSym,
    Matrix_jSym,
    Matrix_matSym,
    Matrix_ncSym,
    Matrix_pSym,
    Matrix_permSym,
    Matrix_statusSym,
    Matrix_uploSym,
    Matrix_xSym,
    Matrix_zSym;

#endif

