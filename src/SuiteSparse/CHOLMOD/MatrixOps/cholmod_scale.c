//------------------------------------------------------------------------------
// CHOLMOD/MatrixOps/cholmod_scale: scale a sparse matrix
//------------------------------------------------------------------------------

// CHOLMOD/MatrixOps Module.  Copyright (C) 2005-2023, Timothy A. Davis.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0+

//------------------------------------------------------------------------------

// scale a matrix:  A = diag(s)*A, A*diag(s), s*A, or diag(s)*A*diag(s)
//
// A can be of any type (packed/unpacked, upper/lower/unsymmetric).
// The symmetry of A is ignored; all entries in the matrix are modified.
//
// If A is m-by-n unsymmetric but scaled symmtrically, the result is
// A = diag (s (1:m)) * A * diag (s (1:n)).
//
// Note: diag(s) should be interpretted as spdiags(s,0,n,n) where n=length(s).
//
// Row or column scaling of a symmetric matrix still results in a symmetric
// matrix, since entries are still ignored by other routines.
// For example, when row-scaling a symmetric matrix where just the upper
// triangular part is stored (and lower triangular entries ignored)
// A = diag(s)*triu(A) is performed, where the result A is also
// symmetric-upper.  This has the effect of modifying the implicit lower
// triangular part.  In MATLAB notation:
//
//      U = diag(s)*triu(A) ;
//      L = tril (U',-1)
//      A = L + U ;
//
// The scale parameter determines the kind of scaling to perform:
//
//      CHOLMOD_SCALAR: s[0]*A
//      CHOLMOD_ROW:    diag(s)*A
//      CHOLMOD_COL:    A*diag(s)
//      CHOLMOD_SYM:    diag(s)*A*diag(s)
//
// The size of S depends on the scale parameter:
//
//      CHOLMOD_SCALAR: size 1
//      CHOLMOD_ROW:    size nrow-by-1 or 1-by-nrow
//      CHOLMOD_COL:    size ncol-by-1 or 1-by-ncol
//      CHOLMOD_SYM:    size max(nrow,ncol)-by-1, or 1-by-max(nrow,ncol)
//
// workspace: none
//
// Real, complex, and zomplex matrices are supported, of any dtype.
// The xtype and dtype of A and S must match.

#include "cholmod_internal.h"

#ifndef NGPL
#ifndef NMATRIXOPS

//------------------------------------------------------------------------------
// t_cholmod_scale_worker
//------------------------------------------------------------------------------

#define DOUBLE
#define REAL
#include "t_cholmod_scale_worker.c"
#define COMPLEX
#include "t_cholmod_scale_worker.c"
#define ZOMPLEX
#include "t_cholmod_scale_worker.c"

#undef  DOUBLE
#define SINGLE
#define REAL
#include "t_cholmod_scale_worker.c"
#define COMPLEX
#include "t_cholmod_scale_worker.c"
#define ZOMPLEX
#include "t_cholmod_scale_worker.c"

//------------------------------------------------------------------------------
// cholmod_scale
//------------------------------------------------------------------------------

int CHOLMOD(scale)
(
    // input:
    cholmod_dense *S,   // scale factors (scalar or vector)
    int scale,          // type of scaling to compute
    // input/output:
    cholmod_sparse *A,  // matrix to scale
    cholmod_common *Common
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (S, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (S, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    if (A->xtype != S->xtype || A->dtype != S->dtype)
    {
        ERROR (CHOLMOD_INVALID, "xtype and dtype of A and S must match") ;
        return (FALSE) ;
    }

    Int ncol = A->ncol ;
    Int nrow = A->nrow ;
    Int sncol = S->ncol ;
    Int snrow = S->nrow ;
    bool ok ;

    if (scale == CHOLMOD_SCALAR)
    {
        ok = (snrow == 1 && sncol == 1) ;
    }
    else if (scale == CHOLMOD_ROW)
    {
        ok = (snrow == nrow && sncol == 1) || (snrow == 1 && sncol == nrow) ;
    }
    else if (scale == CHOLMOD_COL)
    {
        ok = (snrow == ncol && sncol == 1) || (snrow == 1 && sncol == ncol) ;
    }
    else if (scale == CHOLMOD_SYM)
    {
        Int nn = MAX (nrow, ncol) ;
        ok = (snrow == nn && sncol == 1) || (snrow == 1 && sncol == nn) ;
    }
    else
    {
        // scale invalid
        ERROR (CHOLMOD_INVALID, "invalid scaling option") ;
        return (FALSE) ;
    }
    if (!ok)
    {
        // S is wrong size
        ERROR (CHOLMOD_INVALID, "invalid scale factors") ;
        return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    //--------------------------------------------------------------------------
    // scale the matrix
    //--------------------------------------------------------------------------

    switch ((A->xtype + A->dtype) % 8)
    {
        case CHOLMOD_REAL    + CHOLMOD_SINGLE:
            rs_cholmod_scale_worker (S, scale, A) ;
            break ;

        case CHOLMOD_COMPLEX + CHOLMOD_SINGLE:
            cs_cholmod_scale_worker (S, scale, A) ;
            break ;

        case CHOLMOD_ZOMPLEX + CHOLMOD_SINGLE:
            zs_cholmod_scale_worker (S, scale, A) ;
            break ;

        case CHOLMOD_REAL    + CHOLMOD_DOUBLE:
            rd_cholmod_scale_worker (S, scale, A) ;
            break ;

        case CHOLMOD_COMPLEX + CHOLMOD_DOUBLE:
            cd_cholmod_scale_worker (S, scale, A) ;
            break ;

        case CHOLMOD_ZOMPLEX + CHOLMOD_DOUBLE:
            zd_cholmod_scale_worker (S, scale, A) ;
            break ;
    }

    ASSERT (CHOLMOD(dump_sparse) (A, "A scaled", Common) >= 0) ;
    return (TRUE) ;
}

#endif
#endif

