/* Subsequent to changes made for SuiteSparse by Timothy A. Davis,           */
/* which are documented in the file                                          */
/* ../../../../../inst/doc/SuiteSparse/CHOLMOD/SuiteSparse_metis/README.txt, */
/* the METIS library sources, which include this file, have been patched     */
/* for R package Matrix by its authors to resolve warnings issued by GCC     */
/* and Clang with options -Wall and -Wextra.  See the files ssget.sh and     */
/* *.patch below ../../../../../inst/scripts for details.                    */

/*!
\file gk_arch.h
\brief This file contains various architecture-specific declerations

\date   Started 3/27/2007
\author George
\version\verbatim $Id: gk_arch.h 10711 2011-08-31 22:23:04Z karypis $ \endverbatim
*/

#ifndef _GK_ARCH_H_
#define _GK_ARCH_H_

/*************************************************************************
* Architecture-specific differences in header files
**************************************************************************/

// stdint.h, inttypes.h: added for SuiteSparse, Dec 2022
#include <stdint.h>
#include <inttypes.h>

#if 0
#ifdef LINUX
#if !defined(__USE_XOPEN)
#define __USE_XOPEN
#endif
#if !defined(_XOPEN_SOURCE)
#define _XOPEN_SOURCE 600
#endif
#if !defined(__USE_XOPEN2K)
#define __USE_XOPEN2K
#endif
#endif


#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif


#ifdef __MSC__ 
  #include "ms_stdint.h"
  #include "ms_inttypes.h"
  #include "ms_stat.h"
#else
#ifndef SUNOS
  #include <stdint.h>
#endif
  #include <inttypes.h>
  #include <sys/types.h>
  #include <sys/resource.h>
  #include <sys/time.h>
#endif
#endif


/*************************************************************************
* Architecture-specific modifications
**************************************************************************/
// revised for SuiteSparse, Jan 2023
#if defined ( NO_SSIZE_T )
// #ifdef WIN32
typedef ptrdiff_t ssize_t;
#else
// POSIX: ssize_t is defined in sys/types.h
#include <sys/types.h>
#endif


#ifdef SUNOS
#define PTRDIFF_MAX  INT64_MAX
#endif

#if 0
// rint and INFINITY disabled for SuiteSparse, Dec 2022
#ifdef __MSC__
/* MSC does not have rint() function */
#define rint(x) ((int)((x)+0.5))  

/* MSC does not have INFINITY defined */
#ifndef INFINITY
#define INFINITY FLT_MAX
#endif
#endif
#endif

#endif
