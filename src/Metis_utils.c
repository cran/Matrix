#include "Metis_utils.h"

void ssc_metis_order(int n, const int Tp [], const int Ti [],
		     idxtype* perm, idxtype* iperm)
{
    int  j, num_flag = 0, options_flag = 0;
    idxtype *xadj, *adj;

    xadj = (idxtype *) R_alloc(n + 1, sizeof(idxtype));
    adj = (idxtype *) R_alloc(2 * (Tp[n] - n), sizeof(idxtype));
				/* temporarily use perm to store lengths */
    memset(perm, 0, sizeof(idxtype) * n);
    for (j = 0; j < n; j++) {
	int ip, p2 = Tp[j+1];
	for (ip = Tp[j]; ip < p2; ip++) {
	    int i = Ti[ip];
	    if (i != j) {
		perm[i]++;
		perm[j]++;
	    }
	}
    }
    xadj[0] = 0;
    for (j = 0; j < n; j++) xadj[j+1] = xadj[j] + perm[j];
				/* temporarily use perm to store pointers */
    Memcpy(perm, xadj, n);
    for (j = 0; j < n; j++) {
	int ip, p2 = Tp[j+1];
	for (ip = Tp[j]; ip < p2; ip++) {
	    int i = Ti[ip];
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
