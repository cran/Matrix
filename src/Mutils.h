#ifndef MATRIX_MUTILS_H
#define MATRIX_MUTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <Rdefines.h>
#include <Rconfig.h>
#include <R.h>  /* to include Rconfig.h */

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("Matrix", String)
#else
#define _(String) (String)
#endif

SEXP triangularMatrix_validate(SEXP obj);
SEXP symmetricMatrix_validate(SEXP obj);

/* enum constants from cblas.h and some short forms */
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};
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
char uplo_value(SEXP x);
char diag_value(SEXP x);

int csc_unsorted_columns(int ncol, const int p[], const int i[]);
void csc_sort_columns(int ncol, const int p[], int i[], double x[]);
SEXP triple_as_SEXP(int nrow, int ncol, int nz,
		    const int Ti [], const int Tj [], const double Tx [],
		    char *Rclass);
SEXP csc_check_column_sorting(SEXP A);
void csc_compTr(int m, int n, int nnz,
		const int xp[], const int xi[], const double xx[],
		int ap[], int ai[], double ax[]);
void ssc_symbolic_permute(int n, int upper, const int perm[],
			  int Ap[], int Ai[]);
SEXP Matrix_make_named(int TYP, char **names);
SEXP check_scalar_string(SEXP sP, char *vals, char *nm);
double *packed_getDiag(double *dest, SEXP x);
SEXP Matrix_getElement(SEXP list, char *nm);

#define PACKED_TO_FULL(TYPE)						\
TYPE *packed_to_full_ ## TYPE(TYPE *dest, const TYPE *src,		\
			     int n, enum CBLAS_UPLO uplo)
PACKED_TO_FULL(double);
PACKED_TO_FULL(int);
#undef PACKED_TO_FULL

#define FULL_TO_PACKED(TYPE)						\
TYPE *full_to_packed_ ## TYPE(TYPE *dest, const TYPE *src, int n,	\
			      enum CBLAS_UPLO uplo, enum CBLAS_DIAG diag)
FULL_TO_PACKED(double);
FULL_TO_PACKED(int);
#undef FULL_TO_PACKED


extern	 /* stored pointers to symbols initialized in R_init_Matrix */
#include "Syms.h"

/* zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

/* number of elements in one triangle of a square matrix of order n */
#define PACKED_LENGTH(n)   ((n) * ((n) + 1))/2

/* duplicate the slot with name given by sym from src to dest */
#define slot_dup(dest, src, sym)  SET_SLOT(dest, sym, duplicate(GET_SLOT(src, sym)))

#define uplo_P(_x_) CHAR(STRING_ELT(GET_SLOT(_x_, Matrix_uploSym), 0))
#define diag_P(_x_) CHAR(STRING_ELT(GET_SLOT(_x_, Matrix_diagSym), 0))


/**
 * Check for valid length of a packed triangular array and return the
 * corresponding number of columns
 *
 * @param len length of a packed triangular array
 *
 * @return number of columns
 */
static R_INLINE
int packed_ncol(int len)
{
    int disc = 8 * len + 1;	/* discriminant */
    int sqrtd = (int) sqrt((double) disc);

    if (len < 0 || disc != sqrtd * sqrtd)
	error(_("invalid 'len' = %d in packed_ncol"));
    return (sqrtd - 1)/2;
}

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
 * Expand compressed pointers in the array mp into a full set of indices
 * in the array mj.
 *
 * @param ncol number of columns (or rows)
 * @param mp column pointer vector of length ncol + 1
 * @param mj vector of length mp[ncol] to hold the result
 *
 * @return mj
 */
static R_INLINE
int* expand_cmprPt(int ncol, const int mp[], int mj[])
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

/**
 * Apply a permutation to an integer vector
 *
 * @param i vector of 0-based indices
 * @param n length of vector i
 * @param perm 0-based permutation vector of length max(i) + 1
 */
static R_INLINE void
int_permute(int i[], int n, const int perm[])
{
    int j;
    for (j = 0; j < n; j++) i[j] = perm[i[j]];
}

/**
 * Force index pairs to be in the upper triangle of a matrix
 *
 * @param i vector of 0-based row indices
 * @param j vector of 0-based column indices
 * @param nnz length of index vectors
 */
static R_INLINE void
make_upper_triangular(int i[], int j[], int nnz)
{
    int k;
    for (k = 0; k < nnz; k++) {
	if (i[k] > j[k]) {
	    int tmp = i[k];
	    i[k] = j[k];
	    j[k] = tmp;
	}
    }
}

void make_array_triangular(double *x, SEXP from);

SEXP Matrix_expand_pointers(SEXP pP);


/**
 * Elementwise increment dest by src
 *
 * @param dest vector to be incremented
 * @param src vector to be added to dest
 * @param n length of vectors
 *
 * @return dest
 */
static R_INLINE double*
vecIncrement(double dest[], const double src[], int n) {
    int i;
    for (i = 0; i < n; i++) dest[i] += src[i];
    return dest;
}

/**
 * Elementwise sum of src1 and src2 into dest
 *
 * @param dest vector to be incremented
 * @param src1 vector to be added
 * @param src1 second vector to be added
 * @param n length of vectors
 *
 * @return dest
 */
static R_INLINE double*
vecSum(double dest[], const double src1[], const double src2[],
       int n) {
    int i;
    for (i = 0; i < n; i++) dest[i] = src1[i] + src2[i];
    return dest;
}

SEXP alloc_real_classed_matrix(char *class, int nrow, int ncol);

#ifdef __cplusplus
}
#endif

#endif /* MATRIX_MUTILS_H_ */
