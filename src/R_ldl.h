/* This is a modified version of the ldl.h file released by Timothy A. Davis  */
/* in the LDL package and carrying the copyright shown below. The             */
/* modifications are to replace scratch arrays passed as arguments by         */
/* dynamically allocated arrays.                                              */
/* Douglas Bates (Nov., 2004)                                                 */

/* ========================================================================== */
/* === ldl.h:  include file for the LDL package ============================= */
/* ========================================================================== */

/* LDL Version 1.0 (Dec. 31, 2003), Copyright (c) 2003 by Timothy A Davis,
 * University of Florida.  All Rights Reserved.  See README for the License.
 */

#include <R_ext/RS.h>
#include <R_ext/Memory.h>

extern void R_ldl_symbolic(int n, const int Ap[], const int Ai[],
			   int Lp[], int Parent[], const int P[],
                           int Pinv[]);

extern int R_ldl_numeric(int n, const int Ap[], const int Ai[],
			 const double Ax[], const int Lp[], const int Parent[],
			 int Li[], double Lx[], double D[],
			 const int P[], const int Pinv[]);

extern void R_ldl_lsolve(int n, double X[], const int Lp[], const int Li[],
			 const double Lx[]);

extern void R_ldl_dsolve(int n, double X[], const double D[]);

extern void R_ldl_ltsolve(int n, double X[], const int Lp[], const int Li[],
			  const double Lx[]);

extern void R_ldl_perm(int n, double X[], const double B[], const int P[]);
extern void R_ldl_permt(int n, double X[], const double B[], const int P[]);

int R_ldl_valid_perm(int n, const int P[]);
int R_ldl_valid_matrix(int n, const int Ap[], const int Ai[]);
