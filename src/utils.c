/* Utilities needed by some of the libraries */

/* TAUCS utilities */
#include "R.h"
#include "taucs/taucs.h"
				/* timers */
double taucs_wtime() { return 0.0; }
double taucs_ctime() { return 0.0; }
				/* memory allocation */
#undef malloc
#undef calloc
#undef realloc
#undef free

void* taucs_malloc_stub (size_t size)               { return malloc(size); }
void* taucs_calloc_stub (size_t nmemb, size_t size) { return calloc(nmemb,size); }
void* taucs_realloc_stub(void* ptr, size_t size)    { return realloc(ptr,size); }
void  taucs_free_stub   (void* ptr)                 { free(ptr); }

double taucs_allocation_amount()   { return 0.0; }
int    taucs_allocation_count()    { return 0; }
int    taucs_allocation_attempts() { return 0; }
void   taucs_allocation_assert_clean() {}
void   taucs_allocation_mark_clean() {}
void   taucs_allocation_induce_failure(int i) {}
				/* logging */
int
taucs_printf(char *fmt, ...)
{
    return 0;
}
				/* arithmetic constants */
double taucs_get_nan() { return R_NaN; }
double taucs_dzero_const     =  0.0;
double taucs_done_const      =  1.0;
double taucs_dminusone_const = -1.0;
