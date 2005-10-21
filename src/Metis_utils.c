#include <Rdefines.h>
#include "Metis/metis.h"
#include "Metis_utils.h"



/** 
 * Establish a fill-reducing permutation for the sparse symmetric
 * matrix of order n represented by the column pointers Tp and row
 * indices Ti.
 * 
 * @param n  order of the sparse symmetric matrix
 * @param Tp  column pointers (total length n + 1)
 * @param Ti  row indices (total length Tp[n])
 * @param Perm array of length n to hold the permutation
 * @param iPerm array of length n to hold the inverse permutation
 * 
 */
void ssc_metis_order(int n, const int Tp [], const int Ti [],
		     int Perm[], int iPerm[])
{
    int  j, num_flag = 0, options_flag = 0;
    idxtype
	*perm = Calloc(n, idxtype), /* in case idxtype != int */
	*iperm = Calloc(n, idxtype),
	*xadj = Calloc(n+1, idxtype),
	*adj = Calloc(2 * (Tp[n] - n), idxtype);

				/* check row indices for correct range */
    for (j = 0; j < Tp[n]; j++)
      if (Ti[j] < 0 || Ti[j] >= n)
	error(_("row index Ti[%d] = %d is out of range [0,%d]"),
	      j, Ti[j], n - 1);
				/* temporarily use perm to store lengths */
    AZERO(perm, n);
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
