// -*- c++ -*-
//              LAPACK++ 1.1 Linear Algebra Package 1.1
//               University of Tennessee, Knoxvilee, TN.
//            Oak Ridge National Laboratory, Oak Ridge, TN.
//        Authors: J. J. Dongarra, E. Greaser, R. Pozo, D. Walker
//                 (C) 1992-1996 All Rights Reserved
//
//                             NOTICE
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Institutions (University of Tennessee, and Oak Ridge National
// Laboratory) nor the Authors make any representations about the suitability 
// of this software for any purpose.  This software is provided ``as is'' 
// without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.
//
// Modifications Copyright (C) 2000-2000 the R Development Core Team

#ifndef LA_F_NAMES_H
#define LA_F_NAMES_H

#if 1

#define LA_BAND_FACT_DOUBLE_H              "bfd.h"
#define LA_BAND_MAT_DOUBLE_H               "bmd.h"
#define LA_GEN_FACT_DOUBLE_H               "gfd.h"
#define LA_EXCEPTION_H                     "laexcp.h"
#define LA_GEN_MAT_COMPLEX_H               "gmc.h"
#define LA_GEN_MAT_DOUBLE_H                "gmd.h"
#define LA_GEN_MAT_FLOAT_H                 "gmf.h"
#define LA_GEN_MAT_INT_H                   "gmi.h"
#define LA_GEN_MAT_LONG_INT_H              "gmli.h"
#define LA_GENERATE_MAT_DOUBLE_H           "genmd.h"
#define LA_INDEX_H                         "laindex.h"
#define LA_SOLVE_DOUBLE_H                  "laslv.h"
#define LA_LOWER_TRIANG_MAT_DOUBLE_H       "ltgmd.h"
#define LA_SPD_BAND_MAT_DOUBLE_H           "spdbmd.h"
#define LA_SPD_MAT_DOUBLE_H                "spdmd.h"
#define LA_SPD_TRIDIAG_MAT_DOUBLE_H        "spdtrmd.h"
#define LA_SYMM_BAND_FACT_DOUBLE_H         "sybfd.h"
#define LA_SYMM_BAND_MAT_DOUBLE_H          "sybmd.h"
#define LA_SYMM_FACT_DOUBLE_H              "syfd.h"
#define LA_SYMM_MAT_DOUBLE_H               "symd.h"
#define LA_SYMM_TRIDIAG_MAT_DOUBLE_H       "sytrmd.h"
#define LA_TRIDIAG_FACT_DOUBLE_H           "trfd.h"
#define LA_TRIDIAG_MAT_DOUBLE_H            "trmd.h"
#define LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H  "ultgmd.h"
#define LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H  "uutgmd.h"
#define LA_UTIL_H                          "lautil.h"
#define LA_UPPER_TRIANG_MAT_DOUBLE_H       "utgmd.h"
#define LA_VECTOR_COMPLEX_H                "lavc.h"
#define LA_VECTOR_DOUBLE_H                 "lavd.h"
#define LA_VECTOR_LONG_INT_H               "lavli.h"
#define LA_COL_VECTOR_DOUBLE_H             "lacvd.h"
#define LA_ROW_VECTOR_DOUBLE_H             "larvd.h"
#define LA_VECTOR_INT_H                    "lavi.h"
#define VECTOR_COMPLEX_H                   "vc.h"
#define VECTOR_DOUBLE_H                    "vd.h"
#define VECTOR_INT_H                       "vi.h"
#define VECTOR_LONG_INT_H                  "vli.h"

#define VECTOR_FLOAT_H                      "vf.h"

#else   
/* old-style long names */

#define LA_BAND_FACT_DOUBLE_H              "LaBandFactDouble.h"
#define LA_BAND_MAT_DOUBLE_H               "LaBandMatDouble.h"
#define LA_GEN_FACT_DOUBLE_H               "LaGenFactDouble.h"
#define LA_GEN_MAT_DOUBLE_H                "LaGenMatDouble.h"
#define LA_GEN_MAT_INT_H                   "LaGenMatInt.h"
#define LA_LOWER_TRIANG_MAT_DOUBLE_H       "LaLowerTriangMatDouble.h"
#define LA_SPD_BAND_MAT_DOUBLE_H           "LaSpdBandMatDouble.h"
#define LA_SPD_TRIDIAG_MAT_DOUBLE_H        "LaSpdTridiagMatDouble.h"
#define LA_SYMM_BAND_FACT_DOUBLE_H         "LaSymmBandFactDouble.h"
#define LA_SYMM_FACT_DOUBLE_H              "LaSymmFactDouble.h"
#define LA_SYMM_MAT_DOUBLE_H               "LaSymmMatDouble.h"
#define LA_TRIDIAG_FACT_DOUBLE_H           "LaTridiagFactDouble.h"
#define LA_TRIDIAG_MAT_DOUBLE_H            "LaTridiagMatDouble.h"
#define LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H  "LaUnitLowerTriangMatDouble.h"
#define LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H  "LaUnitUpperTriangMatDouble.h"
#define LA_UPPER_TRIANG_MAT_DOUBLE_H       "LaUpperTriangMatDouble.h"
#define LA_VECTOR_DOUBLE_H                 "LaVectorDouble.h"
#define LA_VECTOR_INT_H                    "LaVectorInt.h"
#define LA_UTIL_H                           "lautil"
#define VECTOR_COMPLEX_H                    "VectorComplex.h"
#define VECTOR_DOUBLE_H                    "VectorDouble.h"
#define VECTOR_INT_H                       "VectorInt.h"
#define VECTOR_FLOAT_H                      "VectorFloat.h"

#endif 
#endif 
