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
 * stdheaders.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: stdheaders.h 5993 2009-01-07 02:09:57Z karypis $
 */

#ifndef _LIBMETIS_STDHEADERS_H_
#define _LIBMETIS_STDHEADERS_H_

#include <stdio.h>
#ifdef __STDC__
#include <stdlib.h>
#else
#include <malloc.h>
#endif
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>

#endif
