/* Subsequent to changes made for SuiteSparse by Timothy A. Davis,           */
/* which are documented in the file                                          */
/* ../../../../../inst/doc/SuiteSparse/CHOLMOD/SuiteSparse_metis/README.txt, */
/* the METIS library sources, which include this file, have been patched     */
/* for R package Matrix by its authors to resolve warnings issued by GCC     */
/* and Clang with options -Wall and -Wextra.  See the files ssget.sh and     */
/* *.patch below ../../../../../inst/scripts for details.                    */

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * omp.c
 *
 * This file contains "fake" implementations of OpenMP's runtime libraries
 *
 */

#include "GKlib.h"

#ifdef GK_NOOPENMP  /* remove those for now */
#if !defined(_OPENMP)
void omp_set_num_threads(int num_threads) { return; }
int omp_get_num_threads(void) { return 1; }
int omp_get_max_threads(void) { return 1; }
int omp_get_thread_num(void) { return 0; }
int omp_get_num_procs(void) { return 1; }
int omp_in_parallel(void) { return 0; }
void omp_set_dynamic(int num_threads) { return; }
int omp_get_dynamic(void) { return 0; }
void omp_set_nested(int nested) { return; }
int omp_get_nested(void) { return 0; }
#endif
#endif


