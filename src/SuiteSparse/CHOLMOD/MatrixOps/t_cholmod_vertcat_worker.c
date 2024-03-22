//------------------------------------------------------------------------------
// CHOLMOD/MatrixOps/t_cholmod_vertcat_worker
//------------------------------------------------------------------------------

// CHOLMOD/MatrixOps Module.  Copyright (C) 2005-2023, Timothy A. Davis.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0+

//------------------------------------------------------------------------------

#include "cholmod_template.h"

static void TEMPLATE (cholmod_vertcat_worker)
(
    cholmod_sparse *C,  // output matrix
    cholmod_sparse *A,  // top matrix to concatenate
    cholmod_sparse *B   // bottom matrix to concatenate
)
{

    //--------------------------------------------------------------------------
    // get inputs
    //--------------------------------------------------------------------------

    Int *Ap  = A->p ;
    Int *Anz = A->nz ;
    Int *Ai  = A->i ;
#ifndef PATTERN
    Real *Ax  = A->x ;
#ifdef ZOMPLEX
    Real *Az  = A->z ;
#endif
#endif
    bool apacked = A->packed ;
    Int anrow = A->nrow ;

    Int *Bp  = B->p ;
    Int *Bnz = B->nz ;
    Int *Bi  = B->i ;
#ifndef PATTERN
    Real *Bx  = B->x ;
#ifdef ZOMPLEX
    Real *Bz  = B->z ;
#endif
#endif
    bool bpacked = B->packed ;

    Int *Cp = C->p ;
    Int *Ci = C->i ;
#ifndef PATTERN
    Real *Cx = C->x ;
#ifdef ZOMPLEX
    Real *Cz = C->z ;
#endif
#endif
    Int ncol = C->ncol ;

    //--------------------------------------------------------------------------
    // C = [A ; B]
    //--------------------------------------------------------------------------

    Int pc = 0 ;

    for (Int j = 0 ; j < ncol ; j++)
    {
        // append A(:,j) as the first part of C(:,j)
        Int p = Ap [j] ;
        Int pend = (apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
        Cp [j] = pc ;
        for ( ; p < pend ; p++)
        {
            Ci [pc] = Ai [p] ;
            ASSIGN (Cx, Cz, pc, Ax, Az, p) ;
            pc++ ;
        }

        // append B(:,j) as the second part of C(:,j)
        p = Bp [j] ;
        pend = (bpacked) ? (Bp [j+1]) : (p + Bnz [j]) ;
        for ( ; p < pend ; p++)
        {
            Ci [pc] = Bi [p] + anrow ;
            ASSIGN (Cx, Cz, pc, Bx, Bz, p) ;
            pc++ ;
        }
    }
    Cp [ncol] = pc ;
}

#undef PATTERN
#undef REAL
#undef COMPLEX
#undef ZOMPLEX

