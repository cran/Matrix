#ifndef MATRIX_TRIPLET_TO_COL_H
#define MATRIX_TRIPLET_TO_COL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <R_ext/RS.h>
#include <R.h>  /* to include Rconfig.h */

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("Matrix", String)
#else
#define _(String) (String)
#endif

void triplet_to_col(int n_row, int n_col, int nz,
		    const int Ti[], const int Tj[], const double Tx[],
		    int Ap[], int Ai[], double Ax[]);


#ifdef __cplusplus
}
#endif

#endif /* MATRIX_TRIPLET_TO_COL_H_ */
