/* Based on UMF_triplet.c from UMFPACK which carries the following copyright  */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

#include <R_ext/RS.h>

void triplet_to_col
(
    int n_row,
    int n_col,
    int nz,
    const int Ti [ ],		/* size nz */
    const int Tj [ ],		/* size nz */
    const double Tx [ ],	/* size nz */
    int Ap [ ],			/* size n_col + 1 */
    int Ai [ ],			/* size nz */
    double Ax [ ]		/* size nz */
)
{
    int i, j, k, p, cp, p1, p2, pdest, pj;
    int *Rp = Calloc((n_row + 1), int), 
	*Rj = Calloc(nz, int),
	*W = Calloc((n_row > n_col) ? n_row : n_col, int),
	*RowCount = Calloc(n_row, int);
    double *Rx = Calloc(nz, double);

    /* count the entries in each row (including duplicates) */
    /* use W as workspace for row counts (including duplicates) */
    memset(W, 0, sizeof(int) * n_row);
    for (k = 0 ; k < nz ; k++) {
	i = Ti [k] ;
	j = Tj [k] ;
	if (i < 0 || i >= n_row || j < 0 || j >= n_col)
	    error("entry %d in input has row %d and column %d", k, i, j);
	W [i]++ ;
    }
    /* compute the row pointers */
    Rp [0] = 0 ;
    for (i = 0 ; i < n_row ; i++) {
	Rp [i+1] = Rp [i] + W [i] ;
	W [i] = Rp [i] ;
    }
    /* W is now equal to the row pointers */

    /* ---------------------------------------------------------------------- */
    /* construct the row form */
    /* ---------------------------------------------------------------------- */
    for (k = 0 ; k < nz ; k++) {
	p = W [Ti [k]]++ ;
	Rj [p] = Tj [k] ;
	Rx [p] = Tx [k] ;
    }
    /* Rp stays the same, but W [i] is advanced to the start of row i+1 */
    /* ---------------------------------------------------------------------- */
    /* sum up duplicates */
    /* ---------------------------------------------------------------------- */
    /* use W [j] to hold position in Ri/Rx/Rz of a_ij, for row i [ */
    for (j = 0 ; j < n_col ; j++) {
	W [j] = -1;
    }
    for (i = 0 ; i < n_row ; i++) {
	p1 = Rp [i] ;
	p2 = Rp [i+1] ;
	pdest = p1 ;
	/* At this point, W [j] < p1 holds true for all columns j, */
	/* because Ri/Rx/Rz is stored in row oriented order. */
	for (p = p1; p < p2; p++) {
	    j = Rj [p] ;
	    pj = W [j] ;
	    if (pj >= p1) {
		/* this column index, j, is already in row i, at position pj */
		/* sum the entry */
		Rx [pj] += Rx [p] ;
	    } else {
		/* keep the entry */
		/* also keep track in W[j] of position of a_ij for case above */
		W [j] = pdest ;
		/* no need to move the entry if pdest is equal to p */
		if (pdest != p)
		{
		    Rj [pdest] = j ;
		    Rx [pdest] = Rx [p] ;
		}
		pdest++ ;
	    }
	}
	RowCount [i] = pdest - p1 ;
    }
    /* done using W for position of a_ij ] */
    /* ---------------------------------------------------------------------- */
    /* count the entries in each column */
    /* ---------------------------------------------------------------------- */
    /* [ use W as work space for column counts of A */
    memset(W, 0, sizeof(int) * n_col);
    for (i = 0 ; i < n_row ; i++) {
	for (p = Rp [i] ; p < Rp [i] + RowCount [i] ; p++) {
	    j = Rj [p] ;
	    W [j]++ ;
	}
    }
    /* ---------------------------------------------------------------------- */
    /* create the column pointers */
    /* ---------------------------------------------------------------------- */
    Ap [0] = 0 ;
    for (j = 0 ; j < n_col ; j++) {
	Ap [j+1] = Ap [j] + W [j] ;
    }
    /* done using W as workspace for column counts of A ] */
    for (j = 0 ; j < n_col ; j++) {
	W [j] = Ap [j] ;
    }
    /* ---------------------------------------------------------------------- */
    /* construct the column form */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < n_row ; i++) {
	for (p = Rp [i] ; p < Rp [i] + RowCount [i] ; p++) {
	    cp = W [Rj [p]]++ ;
	    Ai [cp] = i ;
	    Ax [cp] = Rx [p] ;
	}
    }
    Free(Rp); Free(Rj); Free(W); Free(RowCount); Free(Rx);
}
