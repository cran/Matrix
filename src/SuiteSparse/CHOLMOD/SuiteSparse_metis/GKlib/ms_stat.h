/* Subsequent to changes made for SuiteSparse by Timothy A. Davis,           */
/* which are documented in the file                                          */
/* ../../../../../inst/doc/SuiteSparse/CHOLMOD/SuiteSparse_metis/README.txt, */
/* the METIS library sources, which include this file, have been patched     */
/* for R package Matrix by its authors to resolve warnings issued by GCC     */
/* and Clang with options -Wall and -Wextra.  See the files ssget.sh and     */
/* *.patch below ../../../../../inst/scripts for details.                    */

#ifndef _MSC_VER // [
#error "Use this header only with Microsoft Visual C++ compilers!"
#endif // _MSC_VER ]

#ifndef _MS_STAT_H_
#define _MS_STAT_H_

#if _MSC_VER > 1000
#pragma once
#endif

#include <sys/stat.h>
/* Test macros for file types.  */

#define __S_ISTYPE(mode, mask)  (((mode) & S_IFMT) == (mask))

#define S_ISDIR(mode)    __S_ISTYPE((mode), S_IFDIR)
#define S_ISCHR(mode)    __S_ISTYPE((mode), S_IFCHR)
#define S_ISBLK(mode)    __S_ISTYPE((mode), S_IFBLK)
#define S_ISREG(mode)    __S_ISTYPE((mode), S_IFREG)

#endif 
