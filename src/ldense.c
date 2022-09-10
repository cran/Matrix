#include "ldense.h"

/* MJ: no longer needed ... prefer more general (un)?pack() */
#if 0

/* dense logical Matrices "ldenseMatrix" classes --- almost identical to
 * dense nonzero-pattern: "ndenseMatrix" ones
 */

/* this is very close to dspMatrix_as_dsy* () in ./dspMatrix.c : */
SEXP lspMatrix_as_lsyMatrix(SEXP from, SEXP kind)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS(
			   (asInteger(kind) == 1) ? "nsyMatrix" : "lsyMatrix")),
	uplo = GET_SLOT(from, Matrix_uploSym),
	dimP = GET_SLOT(from, Matrix_DimSym),
	dmnP = GET_SLOT(from, Matrix_DimNamesSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_DimNamesSym, duplicate(dmnP));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    idense_unpack(LOGICAL(ALLOC_SLOT(val, Matrix_xSym, LGLSXP, n*n)),
		  LOGICAL( GET_SLOT(from, Matrix_xSym)),
		  n,
		  *CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
		  NUN);
    UNPROTECT(1);
    return val;
}

// this is very close to dsyMatrix_as_lsp*() in ./dsyMatrix.c  -- keep synced !
SEXP lsyMatrix_as_lspMatrix(SEXP from, SEXP kind)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS(
			   (asInteger(kind) == 1) ? "nspMatrix" : "lspMatrix")),
	uplo = GET_SLOT(from, Matrix_uploSym),
	dimP = GET_SLOT(from, Matrix_DimSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    idense_pack(
	LOGICAL(ALLOC_SLOT(val, Matrix_xSym, LGLSXP, (n*(n+1))/2)),
	LOGICAL( GET_SLOT(from, Matrix_xSym)),
	n,
	*CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
	NUN);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(from, Matrix_DimNamesSym)));
    SET_SLOT(val, Matrix_factorSym,
	     duplicate(GET_SLOT(from, Matrix_factorSym)));
    UNPROTECT(1);
    return val;
}

// this is very close to dtpMatrix_as_dtr*() in ./dtpMatrix.c -- keep synced!
SEXP ltpMatrix_as_ltrMatrix(SEXP from, SEXP kind)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS(
			   (asInteger(kind) == 1) ? "ntrMatrix" : "ltrMatrix")),
	uplo = GET_SLOT(from, Matrix_uploSym),
	diag = GET_SLOT(from, Matrix_diagSym),
	dimP = GET_SLOT(from, Matrix_DimSym),
	dmnP = GET_SLOT(from, Matrix_DimNamesSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_DimNamesSym, duplicate(dmnP));
    SET_SLOT(val, Matrix_diagSym, duplicate(diag));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    idense_unpack(LOGICAL(ALLOC_SLOT(val, Matrix_xSym, LGLSXP, n*n)),
		  LOGICAL(GET_SLOT(from, Matrix_xSym)),
		  n,
		  *CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
		  *CHAR(STRING_ELT(diag, 0)) == 'N' ? NUN : UNT);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(from, Matrix_DimNamesSym)));
    UNPROTECT(1);
    return val;
}

/* this is very close to dtrMatrix_as_dtp* () in ./dtrMatrix.c : */
SEXP ltrMatrix_as_ltpMatrix(SEXP from, SEXP kind)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS(
			   (asInteger(kind) == 1) ? "ntpMatrix" : "ltpMatrix")),
	uplo = GET_SLOT(from, Matrix_uploSym),
	diag = GET_SLOT(from, Matrix_diagSym),
	dimP = GET_SLOT(from, Matrix_DimSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_diagSym, duplicate(diag));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    idense_pack(
	LOGICAL(ALLOC_SLOT(val, Matrix_xSym, LGLSXP, (n*(n+1))/2)),
	LOGICAL(GET_SLOT(from, Matrix_xSym)),
	n,
	*CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
	*CHAR(STRING_ELT(diag, 0)) == 'N' ? NUN : UNT);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(from, Matrix_DimNamesSym)));
    UNPROTECT(1);
    return val;
}

#endif

/* MJ: no longer needed ... prefer more general dense_as_general() */
#if 0

/* this is very close to dtrMatrix_as_dge*() :*/
SEXP ltrMatrix_as_lgeMatrix(SEXP from, SEXP kind)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS(
			   (asInteger(kind) == 1) ? "ngeMatrix" : "lgeMatrix"));
    slot_dup(val, from, Matrix_xSym);
    slot_dup(val, from, Matrix_DimSym);
    slot_dup(val, from, Matrix_DimNamesSym);
    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));

    idense_unpacked_make_triangular(LOGICAL(GET_SLOT(val, Matrix_xSym)), from);
    UNPROTECT(1);
    return val;
}

// this is somewhat close to dense_as_general() :
SEXP lsyMatrix_as_lgeMatrix(SEXP from, SEXP kind)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS(
			   (asInteger(kind) == 1) ? "ngeMatrix" : "lgeMatrix"));
    slot_dup(val, from, Matrix_xSym);
    slot_dup(val, from, Matrix_DimSym);
    set_symmetrized_DimNames(val, GET_SLOT(from, Matrix_DimNamesSym), -1);
    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));

    idense_unpacked_make_symmetric(LOGICAL(GET_SLOT(val, Matrix_xSym)), from);
    UNPROTECT(2);
    return val;
}

#endif /* MJ */
