/*------ Definition of a template for [diln]Csparse_subassign(...) : *
 *                       --------     ~~~~~~~~~~~~~~~~~~~~~~
 * i.e., included several times from ./Csparse.c
 *                                   ~~~~~~~~~~~
 *
 _slot_kind : use the integer codes matching  x_slot_kind in ./Mutils.h
 *							       ~~~~~~~~
 */

#ifdef _DOUBLE_x

# define has_x_slot
# define Type_x double
# define STYP_x REAL
# define SXP_x  REALSXP
# define Zero_x 0.

# undef _DOUBLE_x

#elif defined _LOGICAL_x

// --> lgCMatrix etc
# define has_x_slot
# define Type_x int
# define STYP_x LOGICAL
# define SXP_x  LGLSXP
# define Zero_x FALSE

# undef _LOGICAL_x

#elif defined _PATTERN_x

// --> ngCMatrix etc
// # undef has_x_slot
# define Type_x int
# define STYP_x LOGICAL
# define SXP_x  LGLSXP
# define Zero_x FALSE

# undef _PATTERN_x

#elif defined _COMPLEX_x // not yet existing

// --> igCMatrix etc
# define has_x_slot
# define Type_x Rcomplex
# define STYP_x COMPLEX
# define SXP_x  CPLXSXP
# define Zero_x {0., 0.} // FIXME !

# undef _COMPLEX_x

#elif defined _INTEGER_x // not yet existing

// --> igCMatrix etc
# define has_x_slot
# define Type_x int
# define STYP_x INTEGER
# define SXP_x  INTSXP
# define Zero_x 0

# undef _INTEGER_x

#endif

//--------------- The code, included inside a switch() case, ------------------
{
    Type_x* xx = STYP_x(x);
#ifdef has_x_slot
    Type_x* rx = Calloc(nnz, Type_x); // to become  x slot
#endif
    ii = 0; // ii in  0..(n-1)
    for(int j=0; j < nc; ) { // look at 0-based column 'j' <=> 1-based R: x[, j+1]
	int nr_j = rp[j]; // cumulative number of non-zero entries in this column 'j'
	for(int i=0; i < nr; i++, ii++) {
	    // look at 0-based row 'i' -- C's x[ii]  == R's  x[i+1, j+1]
	    if(xx[ii] != Zero_x) {
		ri[nz] = i; // 0-based row number
#ifdef has_x_slot
		rx[nz] = xx[ii];
#endif
		nr_j++;
		if(++nz >= nnz && ii < n-1) {// increase nnz and grow both 'rx' and 'ri'
		    // current density ~= nz / ii == estim.final dens. ==> est. nnz = nz*n/ii
		    nnz = imax2(nnz+256, imax2(5*nnz/4, (nz * n) / ii));
		    ri = Realloc(ri, nnz, int);
#ifdef has_x_slot
		    rx = Realloc(rx, nnz, Type_x);
#endif
		}
	    }
	}
	rp[++j] = nr_j;
    }

    // final number of non zeros: nz, almost always *smaller* than 'nnz':
    nnz = nz;
#ifdef has_x_slot
    Memcpy( STYP_x(ALLOC_SLOT(ans, Matrix_xSym,   SXP_x, nnz)), rx, nnz);
    Free(rx);
#endif
}

// clean up remaining defines from header
#ifdef has_x_slot
# undef has_x_slot
#endif
#undef Type_x
#undef STYP_x
#undef SXP_x
#undef Zero_x
