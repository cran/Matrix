/* ========================================================================== */
/* === umfpack_global ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.5, Copyright (c) 2005 by Timothy A. Davis.  CISE Dept,   */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* prototypes for global variables, and basic operators for complex values  */

extern double (*umfpack_hypot) (double, double) ;
extern
int (*umfpack_divcomplex) (double, double, double, double, double *, double *) ;

double umf_hypot (double x, double y) ;
int umf_divcomplex (double, double, double, double, double *, double *) ;

