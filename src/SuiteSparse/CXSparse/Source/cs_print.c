// CXSparse/Source/cs_print: print a sparse matrix
// CXSparse, Copyright (c) 2006-2022, Timothy A. Davis. All Rights Reserved.
// SPDX-License-Identifier: LGPL-2.1+
#include "cs.h"
#include <R_ext/Print.h>

/* print a sparse matrix; use %g for integers to avoid differences with CS_INT */
CS_INT cs_print (const cs *A, CS_INT brief)
{
    CS_INT p, j, m, n, nzmax, nz, *Ap, *Ai ;
    CS_ENTRY *Ax ;
    if (!A) { Rprintf ("(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    Rprintf ("CXSparse Version %d.%d.%d, %s.  %s\n", CS_VER, CS_SUBVER,
        CS_SUBSUB, CS_DATE, CS_COPYRIGHT) ;
    if (nz < 0)
    {
        Rprintf ("%g-by-%g, nzmax: %g nnz: %g, 1-norm: %g\n", (double) m,
            (double) n, (double) nzmax, (double) (Ap [n]), cs_norm (A)) ;
        for (j = 0 ; j < n ; j++)
        {
            Rprintf ("    col %g : locations %g to %g\n", (double) j, 
                (double) (Ap [j]), (double) (Ap [j+1]-1)) ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                Rprintf ("      %g : ", (double) (Ai [p])) ;
#ifdef CS_COMPLEX
                Rprintf ("(%g, %g)\n",
                    Ax ? CS_REAL (Ax [p]) : 1, Ax ? CS_IMAG (Ax [p]) : 0) ;
#else
                Rprintf ("%g\n", Ax ? Ax [p] : 1) ;
#endif
                if (brief && p > 20) { Rprintf ("  ...\n") ; return (1) ; }
            }
        }
    }
    else
    {
        Rprintf ("triplet: %g-by-%g, nzmax: %g nnz: %g\n", (double) m,
            (double) n, (double) nzmax, (double) nz) ;
        for (p = 0 ; p < nz ; p++)
        {

            Rprintf ("    %g %g : ", (double) (Ai [p]), (double) (Ap [p])) ;
#ifdef CS_COMPLEX
            Rprintf ("(%g, %g)\n",
                Ax ? CS_REAL (Ax [p]) : 1, Ax ? CS_IMAG (Ax [p]) : 0) ;
#else
            Rprintf ("%g\n", Ax ? Ax [p] : 1) ;
#endif
            if (brief && p > 20) { Rprintf ("  ...\n") ; return (1) ; }
        }
    }
    return (1) ;
}
