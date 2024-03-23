//------------------------------------------------------------------------------
// CHOLMOD/Utility/cholmod_l_malloc: malloc (int64 version)
//------------------------------------------------------------------------------

// CHOLMOD/Utility Module. Copyright (C) 2023, Timothy A. Davis, All Rights
// Reserved.
// SPDX-License-Identifier: LGPL-2.1+

//------------------------------------------------------------------------------

#define CHOLMOD_ALLOC_FUNCTION      cholmod_l_malloc
#define SUITESPARSE_ALLOC_FUNCTION  SuiteSparse_malloc
#define CHOLMOD_INT64
#include "t_cholmod_malloc.c"

