#include "Metis_utils.h"

void ssc_metis_order(int n, int nnz,
		     const int Tp [], const int Ti [],
		     idxtype* perm, idxtype* iperm)
{
    int i, j, ip, nnodes, num_flag = 0, options_flag = 0;
    idxtype *xadj, *adj;

    nnodes = 2 * (nnz - n);
    xadj = (idxtype *) R_alloc(n + 1, sizeof(idxtype));
    adj = (idxtype *) R_alloc(nnodes, sizeof(idxtype));
				/* temporarily use perm to store lengths */
    for (j = 0; j < n; j++) perm[j] = 0;
    for (j = 0; j < n; j++) {
	for (ip = Tp[j]; ip < Tp[j+1]; ip++) {
	    i = Ti[ip];
	    if (i != j) {
		perm[i]++;
		perm[j]++;
	    }
	}
    }
    xadj[0] = 0;
    for (j = 0; j < n; j++) xadj[j+1] = xadj[j] + perm[j];
				/* temporarily use perm to store pointers */
    for (j = 0; j < n; j++) perm[j] = xadj[j];
    for (j = 0; j < n; j++) {
	for (ip = Tp[j]; ip < Tp[j+1]; ip++) {
	    i = Ti[ip];
	    if (i != j) {
		adj[perm[i]] = j;
		adj[perm[j]] = i;
		perm[i]++;
		perm[j]++;
	    }
	}
    }
    METIS_NodeND(&n, xadj, adj, &num_flag, &options_flag, perm, iperm);
}
