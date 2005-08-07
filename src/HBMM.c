#include "HBMM.h"
#include "iohb.h"
#include "mmio.h"

SEXP Matrix_readHarwellBoeing(SEXP filename)
{
    char *fnm = CHAR(asChar(filename)), *Type = Calloc(4, char);
    int M, N, nz, Nrhs;
    SEXP ans = R_NilValue;

    readHB_info(fnm, &M, &N, &nz, &Type, &Nrhs);
    Rprintf("Filename: %s, M=%d, N=%d, nz=%d, Type=\"%s\", Nrhs=%d\n",
	    fnm, M, N, nz, Type, Nrhs);

    if (toupper(Type[0]) == 'R') { /* Real (double precision) matrix */
	char upT1 = toupper(Type[1]);
	int *dims;

	if (upT1 == 'S') {
	    ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsCMatrix")));
	    SET_SLOT(ans, Matrix_uploSym, mkString("L"));
	}
	if (upT1 == 'R' || upT1 == 'U')
	    ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix")));
	if (ans == R_NilValue) {Free(Type); return ans;}
	SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, N + 1));
	SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nz));
	SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nz));
	SET_SLOT(ans, Matrix_DimSym, allocVector(INTSXP, 2));
	dims = INTEGER(GET_SLOT(ans, Matrix_DimSym));
	dims[0] = M; dims[1] = N;
	readHB_mat_double(fnm, INTEGER(GET_SLOT(ans, Matrix_pSym)),
			  INTEGER(GET_SLOT(ans, Matrix_iSym)),
			  REAL(GET_SLOT(ans, Matrix_xSym)));
    }
    Free(Type);
    UNPROTECT(1);
    return ans;
}

SEXP Matrix_readMatrixMarket(SEXP filename)
{
    FILE *conn; 
    MM_typecode code;
    int *dims, M, N, i, nz;
    SEXP ans = R_NilValue;

    if (isString(filename)) {
	conn = fopen(CHAR(asChar(filename)), "r");
	if (conn == NULL) {
	    error("Unable to open file: %s", CHAR(asChar(filename)));
	    return R_NilValue;
	}
    } else {
	error("non-string values not presently accepted");
    }

    if ((nz = mm_read_banner(conn, &code))) {
	fclose(conn);
	error("mm_read_banner returned code %d", nz);
    }
    if (!mm_is_valid(code)) {
	fclose(conn);
	error("Invalid code: %s", mm_typecode_to_str(code));
    }

    if (mm_is_sparse(code)) mm_read_mtx_crd_size(conn, &M, &N, &nz);
    if (mm_is_dense(code)) mm_read_mtx_array_size(conn, &M, &N);
    
    if (mm_is_real(code)) {
	if (mm_is_sparse(code)) {
	    int *ii, *jj;
	    double *vv;

	    if (mm_is_general(code))
		ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgTMatrix")));
	    if (mm_is_symmetric(code)) {
		ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsTMatrix")));
		SET_SLOT(ans, Matrix_uploSym, mkString("L"));
	    }
	    if (ans == R_NilValue)
		error("Unrecognized matrix type: %s",
		       mm_typecode_to_str(code));
	    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nz));
	    ii = INTEGER(GET_SLOT(ans, Matrix_iSym));
	    SET_SLOT(ans, Matrix_jSym, allocVector(INTSXP, nz));
	    jj = INTEGER(GET_SLOT(ans, Matrix_jSym));
	    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nz));
	    vv = REAL(GET_SLOT(ans, Matrix_xSym));
	    for (i = 0; i < nz; i++) {
		if (fscanf(conn, "%d %d %lg", &ii[i], &jj[i], &vv[i]) != 3)
		    error("Premature end of file");
		ii[i] -= 1;	/* change indices to zero-based */
		jj[i] -= 1;
	    }
	}
	if (mm_is_dense(code)) {
	    int j;
	    double *vv;

	    if (mm_is_general(code))
		ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
	    /* Need to adjust reading loop for symmetric matrices */
/* 	    if (mm_is_symmetric(code)) { */
/* 		ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsyMatrix"))); */
/* 		SET_SLOT(ans, Matrix_uploSym, mkString("L")); */
/* 	    } */
	    if (ans == R_NilValue)
		error("Unrecognized matrix type: %s",
		       mm_typecode_to_str(code));
	    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, M * N));
	    vv = REAL(GET_SLOT(ans, Matrix_xSym));
	    for (j = 0; j < N; j++) {
		for (i = 0; i < M; i++) {
		    if (fscanf(conn, "%lg", &vv[i + j * M]) != 1)
			error("Premature end of file");
		}
	    }
	}
    } else {
	error("Only real matrices handled at this point");
    }
    SET_SLOT(ans, Matrix_DimSym, allocVector(INTSXP, 2));
    dims = INTEGER(GET_SLOT(ans, Matrix_DimSym));
    dims[0] = M; dims[1] = N;

    fclose(conn);
    UNPROTECT(1);
    return ans;
}

SEXP Matrix_writeHarwellBoeing(SEXP obj, SEXP file, SEXP typep)
{
    char *type = CHAR(asChar(typep)), Type[4] = "RUA";
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym)), *ii, *pp;
    int M = dims[0], N = dims[1], nz;
    double *xx;

    if (type[2] == 'C' || type[2] == 'T') {
	SEXP islot = GET_SLOT(obj, Matrix_iSym);
	nz = LENGTH(islot);
	ii = INTEGER(islot);
	if (type[2] == 'T') {	/* create column pointers */
	    int *i1 = Calloc(nz, int);
	    double *x1 = Calloc(nz, double);
	    
	    pp = Calloc(N + 1, int);
	    triplet_to_col(M, N, nz, ii,
			   INTEGER(GET_SLOT(obj, Matrix_jSym)), xx,
			   pp, i1, x1);
	    nz = pp[N];
	    xx = x1;
	    ii = i1;
	} else pp = INTEGER(GET_SLOT(obj, Matrix_pSym));
    } else error("Only types 'C' and 'T' allowed");

    if (type[0] == 'D') {
	xx = REAL(GET_SLOT(obj, Matrix_xSym));
    } else error("Only real matrices allowed");

    if (!isString(file)) 
	error("non-string values for file not presently accepted");

    if (type[1] == 'S') {
	if (CHAR(STRING_ELT(GET_SLOT(obj, Matrix_uploSym), 0))[0] != 'L')
	    error("Symmetric matrices must be stored in lower triangle");
	Type[1] = 'S';
    }

    writeHB_mat_double(CHAR(asChar(file)), M, N, nz, pp, ii, xx, 0,
		       (double *)NULL, (double *)NULL, (double *)NULL, 
		       "", "", Type, (char*)NULL, (char*)NULL,
		       (char*)NULL, (char*)NULL, "RUA"); 

    if (type[2] == 'T') {Free(ii); Free(pp); Free(xx);}
    return R_NilValue;
}

SEXP Matrix_writeMatrixMarket(SEXP obj, SEXP file, SEXP typep)
{
    char *type = CHAR(asChar(typep));
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym)), *ii,
	*jj = (int *) NULL;
    int M = dims[0], N = dims[1], nz;
    MM_typecode matcode;
    double *xx;

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

    if (!isString(file)) 
	error("non-string values for file not currently allowed");

    if (type[1] == 'S') {
	if (CHAR(STRING_ELT(GET_SLOT(obj, Matrix_uploSym), 0))[0] != 'L')
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

    mm_write_mtx_crd(CHAR(STRING_ELT(file, 0)), M, N, nz, ii, jj, xx,
		     matcode);

    if (type[2] == 'C') Free(jj);
    return R_NilValue;
    
}
