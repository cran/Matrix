#include <Rdefines.h>
#include "metis.h"
#include "Metis_utils.h"

void ssc_metis_order(int n, const int Tp [], const int Ti [],
		     int Perm[], int iPerm[])
{
    int  j, num_flag = 0, options_flag = 0;
    idxtype
	*perm = Calloc(n, idxtype), /* in case idxtype != int */
	*iperm = Calloc(n, idxtype),
	*xadj = Calloc(n+1, idxtype),
	*adj = Calloc(2 * (Tp[n] - n), idxtype);

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
    for (j = 0; j < n; j++) {
	Perm[j] = (int) perm[j];
	iPerm[j] = (int) iperm[j];
    }
    Free(iperm); Free(perm); Free(xadj); Free(adj);
}
