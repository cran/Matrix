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
SEXP get_factors(SEXP obj, char *nm);
SEXP set_factors(SEXP obj, SEXP val, char *nm);
SEXP dgCMatrix_set_Dim(SEXP x, int nrow);
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
SEXP Matrix_make_named(int TYP, char **names);
				/* stored pointers to symbols */
				/* initialized in R_init_Matrix */
extern
#include "Syms.h"

/* zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

/** 
 * Allocate an SEXP of given type and length, assign it as slot nm in
 * the object, and return the SEXP.  The validity of this function
 * depends on SET_SLOT not duplicating val when NAMED(val) == 0.  If
 * this behavior changes then ALLOC_SLOT must use SET_SLOT followed by
 * GET_SLOT to ensure that the value returned is indeed the SEXP in
 * the slot.
 * 
 * @param obj object in which to assign the slot
 * @param nm name of the slot, as an R name object
 * @param type type of SEXP to allocate
 * @param length length of SEXP to allocate
 * 
 * @return SEXP of given type and length assigned as slot nm in obj
 */
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length)
{
    SEXP val = allocVector(type, length);

    SET_SLOT(obj, nm, val);
    return val;
}

/** 
 * Expand the column pointers in the array mp into a full set of column indices
 * in the array mj.
 * 
 * @param ncol number of columns
 * @param mp column pointer vector of length ncol + 1
 * @param mj vector of length mp[ncol] - 1 to hold the result
 * 
 * @return mj
 */
static R_INLINE
int* expand_column_pointers(int ncol, const int mp[], int mj[])
{
    int j;
    for (j = 0; j < ncol; j++) {
	int j2 = mp[j+1], jj;
	for (jj = mp[j]; jj < j2; jj++) mj[jj] = j;
    }
    return mj;
}


/** 
 * Return the linear index of the (row, col) entry in a csc structure.
 * If the entry is not found and missing is 0 an error is signaled;
 * otherwise the missing value is returned.
 * 
 * @param p vector of column pointers
 * @param i vector of row indices
 * @param row row index
 * @param col column index
 * @param missing the value to return is the row, col entry does not
 * exist.  If this is zero and the row, col entry does not exist an
 * error is signaled.
 * 
 * @return index of element at (row, col) if it exists, otherwise missing
 */
static R_INLINE int
check_csc_index(const int p[], const int i[], int row, int col, int missing)
{
    int k, k2 = p[col + 1];
				/* linear search - perhaps replace by bsearch */
    for (k = p[col]; k < k2; k++) if (i[k] == row) return k;
    if (!missing)
	error("row %d and column %d not defined in rowind and colptr",
	      row, col);
    return missing;
}

SEXP alloc3Darray(SEXPTYPE mode, int nrow, int ncol, int nface);

/**
 * Calculate the zero-based index in a row-wise packed lower triangular matrix.
 * This is used for the arrays of blocked sparse matrices.
 *
 * @param i column number (zero-based)
 * @param k row number (zero-based)
 *
 * @return The index of the (k,i) element of a packed lower triangular matrix
 */
static R_INLINE
int Lind(int k, int i)
{
    if (k < i) error("Lind(k = %d, i = %d) must have k >= i", k, i);
    return (k * (k + 1))/2 + i;
}

/** 
 * Check for a complete match on matrix dimensions
 * 
 * @param xd dimensions of first matrix
 * @param yd dimensions of second matrix
 * 
 * @return 1 if dimensions match, otherwise 0
 */
static R_INLINE
int match_mat_dims(const int xd[], const int yd[])
{
    return xd[0] == yd[0] && xd[1] == yd[1];
}

double *expand_csc_column(double *dest, int m, int j,
			  const int Ap[], const int Ai[], const double Ax[]);

#endif

