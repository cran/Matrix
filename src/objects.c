#include "Mdefines.h"
#include "objects.h"

SEXP newObject(const char *what)
{
	SEXP class = PROTECT(R_do_MAKE_CLASS(what)), obj = R_do_new_object(class);
	UNPROTECT(1);
	return obj;
}

char typeToKind(SEXPTYPE type)
{
	switch (type) {
	case LGLSXP:
		return 'l';
	case INTSXP:
		return 'i';
	case REALSXP:
		return 'd';
	case CPLXSXP:
		return 'z';
	default:
		error(_("unexpected type \"%s\" in '%s'"), type2char(type), __func__);
		return '\0';
	}
}

SEXPTYPE kindToType(char kind)
{
	switch (kind) {
	case 'n':
	case 'l':
		return LGLSXP;
	case 'i':
		return INTSXP;
	case 'd':
		return REALSXP;
	case 'z':
		return CPLXSXP;
	default:
		error(_("unexpected kind \"%c\" in '%s'"), kind, __func__);
		return NILSXP;
	}
}

size_t kindToSize(char kind)
{
	switch (kind) {
	case 'n':
	case 'l':
	case 'i':
		return sizeof(int);
	case 'd':
		return sizeof(double);
	case 'z':
		return sizeof(Rcomplex);
	default:
		error(_("unexpected kind \"%c\" in '%s'"), kind, __func__);
		return 0;
	}
}

const char *Matrix_nonvirtual(SEXP obj, int strict)
{
	if (!IS_S4_OBJECT(obj))
		return "";
	static const char *valid[] = { VALID_NONVIRTUAL, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		return "";
	if (!strict)
		ivalid += VALID_NONVIRTUAL_SHIFT(ivalid, 1);
	return valid[ivalid];
}

char Matrix_kind(SEXP obj)
{
	if (IS_S4_OBJECT(obj)) {
		static const char *valid[] = { VALID_NONVIRTUAL, "" };
		int ivalid = R_check_class_etc(obj, valid);
		if (ivalid < 0)
			return '\0';
		ivalid += VALID_NONVIRTUAL_SHIFT(ivalid, 1);
		const char *cl = valid[ivalid];
		return (cl[2] == 'd') ? 'n' : cl[0];
	} else {
		switch (TYPEOF(obj)) {
		case LGLSXP:
			return 'l';
		case INTSXP:
			return 'i';
		case REALSXP:
			return 'd';
		case CPLXSXP:
			return 'z';
		default:
			return '\0';
		}
	}
}

char Matrix_shape(SEXP obj)
{
	if (!IS_S4_OBJECT(obj))
		return '\0';
	static const char *valid[] = { VALID_NONVIRTUAL, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		return '\0';
	ivalid += VALID_NONVIRTUAL_SHIFT(ivalid, 1);
	const char *cl = valid[ivalid];
	return (cl[2] == 'd' || cl[3] != 'M') ? 'g' : cl[1];
}

char Matrix_repr(SEXP obj)
{
	if (!IS_S4_OBJECT(obj))
		return '\0';
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		return '\0';
	ivalid += VALID_NONVIRTUAL_SHIFT(ivalid, 1);
	const char *cl = valid[ivalid];
	switch (cl[2]) {
	case 'e':
	case 'y':
	case 'r':
		return 'u'; /* unpackedMatrix */
	case 'p':
		return 'p'; /* packedMatrix */
	case 'C':
	case 'R':
	case 'T':
		return cl[2]; /* [CRT]sparseMatrix */
	case 'i':
		return 'd'; /* diagonalMatrix */
	case 'd':
		return 'i'; /* indMatrix */
	default:
		return '\0';
	}
}

SEXP R_Matrix_nonvirtual(SEXP obj, SEXP strict)
{
	return mkString(Matrix_nonvirtual(obj, asLogical(strict)));
}

#define RETURN_AS_STRSXP(_C_) \
do { \
	char c = _C_; \
	if (!c) \
		return mkString(""); \
	else { \
		char s[] = { c, '\0' }; \
		return mkString(s); \
	} \
} while (0)

SEXP R_Matrix_kind(SEXP obj)
{
	RETURN_AS_STRSXP(Matrix_kind (obj));
}

SEXP R_Matrix_shape(SEXP obj)
{
	RETURN_AS_STRSXP(Matrix_shape(obj));
}

SEXP R_Matrix_repr(SEXP obj)
{
	RETURN_AS_STRSXP(Matrix_repr (obj));
}
