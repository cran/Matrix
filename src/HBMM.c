#include "HBMM.h"

#if 0
#include "iohb.h"
#endif

#include "mmio.h"

#if 0
SEXP Matrix_writeHarwellBoeing(SEXP obj, SEXP file, SEXP typep)
{
    char *type = CHAR(asChar(typep)), *Type = strdup("RUA");
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym)),
	*ii = (int *) NULL, *pp = (int *) NULL;
    int M = dims[0], N = dims[1], nz = -1;
    double *xx = (double *) NULL;

    if (type[2] == 'C') {
	SEXP islot = GET_SLOT(obj, Matrix_iSym);
	nz = LENGTH(islot);
	ii = INTEGER(islot);
	pp = INTEGER(GET_SLOT(obj, Matrix_pSym));
    } else error("Only type 'C' allowed");

    if (type[0] == 'D') {
	xx = REAL(GET_SLOT(obj, Matrix_xSym));
    } else error("Only real matrices allowed");

    if (!isString(file))
	error("non-string values for file not presently accepted");

    if (type[1] == 'S') {
	if (*uplo_P(obj) != 'L')
	    error("Symmetric matrices must be stored in lower triangle");
	Type[1] = 'S';
    }

    writeHB_mat_double(CHAR(asChar(file)), M, N, nz, pp, ii, xx, 0,
		       (double *)NULL, (double *)NULL, (double *)NULL,
		       "", "", Type, (char*)NULL, (char*)NULL,
		       (char*)NULL, (char*)NULL, "RUA");

    Free(Type);
    return R_NilValue;
}
#endif

SEXP Matrix_writeMatrixMarket(SEXP obj, SEXP file, SEXP typep)
{
    char *type = CHAR(asChar(typep));
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym)),
	*ii = (int *) NULL, *jj = (int *) NULL;
    int M = dims[0], N = dims[1], nz = -1;
    MM_typecode matcode;
    double *xx = (double *) NULL;

    mm_set_matrix(&matcode);
    if (type[2] == 'C' || type[2] == 'T') {
	SEXP islot = GET_SLOT(obj, Matrix_iSym);
	nz = LENGTH(islot);
	ii = INTEGER(islot);
	mm_set_coordinate(&matcode);
    } else error("Only types 'C' and 'T' allowed");

    if (type[0] == 'D') {
	xx = REAL(GET_SLOT(obj, Matrix_xSym));
	mm_set_real(&matcode);
    } else error("Only real matrices allowed");

    if (type[1] == 'S') {
	if (*uplo_P(obj) != 'L')
	    error("Symmetric matrices must be stored in lower triangle");
	mm_set_symmetric(&matcode);
    }
    if (type[1] == 'G') mm_set_general(&matcode);

    if (type[2] == 'C')
	jj = expand_cmprPt(N, INTEGER(GET_SLOT(obj, Matrix_pSym)),
			   Calloc(nz, int));
    if (type[2] == 'T')
	jj = INTEGER(GET_SLOT(obj, Matrix_jSym));
    if (!jj) error("storage mode must be T or C");

    mm_write_mtx_crd(CHAR(asChar(file)), M, N, nz, ii, jj, xx,
		     matcode);

    if (type[2] == 'C') Free(jj);
    return R_NilValue;

}
