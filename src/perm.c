#include "Mdefines.h"
#include "perm.h"

int isPerm(const int *p, int n, int off)
{
	int res = 1;
	if (n <= 0)
		return res;
	int i, j;
	char *work;
	Matrix_Calloc(work, n, char);
	for (i = 0; i < n; ++i) {
		if (p[i] == NA_INTEGER || (j = p[i] - off) < 0 || j >= n || work[j]) {
			res = 0;
			break;
		}
		work[j] = 1;
	}
	Matrix_Free(work, n);
	return res;
}

int signPerm(const int *p, int n, int off)
{
	if (!isPerm(p, n, off))
		error(_("attempt to get sign of non-permutation"));
	int sign = 1;
	if (n <= 0)
		return sign;
	int i, pos = 0;
	char *work;
	Matrix_Calloc(work, n, char);
	while (pos < n) {
		work[pos] = 1;
		i = p[pos] - off;
		while (!work[i]) { /* transposition */
			sign = -sign;
			work[i] = 1;
			i = p[i] - off;
		}
		while (pos < n && work[pos])
			++pos;
	}
	Matrix_Free(work, n);
	return sign;
}

void invertPerm(const int *p, int *ip, int n, int off, int ioff)
{
	if (!isPerm(p, n, off))
		error(_("attempt to invert non-permutation"));
	int j;
	for (j = 0; j < n; ++j)
		ip[p[j] - off] = j + ioff;
	return;
}

void asPerm(const int *p, int *ip, int m, int n, int off, int ioff)
{
	int i, j, tmp;
	for (i = 0; i < n; ++i)
		ip[i] = i + ioff;
	for (i = 0; i < m; ++i) {
		j = p[i] - off;
		if (j < 0 || j >= n)
			error(_("invalid transposition vector"));
		if (j != i) {
			tmp = ip[j];
			ip[j] = ip[i];
			ip[i] = tmp;
		}
	}
	return;
}

SEXP R_isPerm(SEXP p, SEXP off)
{
	if (TYPEOF(p) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "p", "integer");
	if (TYPEOF(off) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "off", "integer");
	if (XLENGTH(off) != 1)
		error(_("'%s' does not have length %d"), "off", 1);
	int off_ = INTEGER(off)[0];
	if (off_ == NA_INTEGER)
		error(_("'%s' is NA"), "off");
	R_xlen_t n_ = XLENGTH(p);
	if (n_ > INT_MAX)
		return ScalarLogical(0);
	return ScalarLogical(isPerm(INTEGER(p), (int) n_, off_));
}

SEXP R_signPerm(SEXP p, SEXP off)
{
	if (TYPEOF(p) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "p", "integer");
	if (TYPEOF(off) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "off", "integer");
	if (XLENGTH(off) != 1)
		error(_("'%s' does not have length %d"), "off", 1);
	int off_ = INTEGER(off)[0];
	if (off_ == NA_INTEGER)
		error(_("'%s' is NA"), "off");
	R_xlen_t n_ = XLENGTH(p);
	if (n_ > INT_MAX)
		error(_("attempt to get sign of non-permutation"));
	return ScalarInteger(signPerm(INTEGER(p), (int) n_, off_));
}

SEXP R_invertPerm(SEXP p, SEXP off, SEXP ioff)
{
	if (TYPEOF(p) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "p", "integer");
	if (TYPEOF(off) != INTSXP || TYPEOF(ioff) != INTSXP)
		error(_("'%s' or '%s' is not of type \"%s\""), "off", "ioff", "integer");
	if (XLENGTH(off) != 1 || XLENGTH(ioff) != 1)
		error(_("'%s' or '%s' does not have length %d"), "off", "ioff", 1);
	int off_ = INTEGER(off)[0], ioff_ = INTEGER(ioff)[0];
	if (off_ == NA_INTEGER || ioff_ == NA_INTEGER)
		error(_("'%s' or '%s' is NA"), "off", "ioff");
	R_xlen_t n_ = XLENGTH(p);
	if (n_ > INT_MAX)
		error(_("attempt to invert non-permutation"));
	SEXP ip = PROTECT(allocVector(INTSXP, n_));
	invertPerm(INTEGER(p), INTEGER(ip), (int) n_, off_, ioff_);
	UNPROTECT(1);
	return ip;
}

SEXP R_asPerm(SEXP p, SEXP off, SEXP ioff, SEXP n)
{
	if (TYPEOF(p) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "p", "integer");
	R_xlen_t m_ = XLENGTH(p);
	if (m_ > INT_MAX)
		error(_("'%s' has length exceeding %s"), "p", "2^31-1");
	if (TYPEOF(off) != INTSXP || TYPEOF(ioff) != INTSXP)
		error(_("'%s' or '%s' is not of type \"%s\""), "off", "ioff", "integer");
	if (XLENGTH(off) != 1 || XLENGTH(ioff) != 1)
		error(_("'%s' or '%s' does not have length %d"), "off", "ioff", 1);
	int off_ = INTEGER(off)[0], ioff_ = INTEGER(ioff)[0];
	if (off_ == NA_INTEGER || ioff_ == NA_INTEGER)
		error(_("'%s' or '%s' is NA"), "off", "ioff");
	if (TYPEOF(n) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "n", "integer");
	if (XLENGTH(n) != 1)
		error(_("'%s' does not have length %d"), "n", 1);
	int n_ = INTEGER(n)[0];
	if (n_ == NA_INTEGER || n_ < m_)
		error(_("'%s' is NA or less than %s"), "n", "length(p)");
	SEXP ip = PROTECT(allocVector(INTSXP, n_));
	asPerm(INTEGER(p), INTEGER(ip), (int) m_, n_, off_, ioff_);
	UNPROTECT(1);
	return ip;
}
