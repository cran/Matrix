#ifndef R_MATRIX_VERSION_H
#define R_MATRIX_VERSION_H

/* Users wanting to do version comparison will include Rversion.h then do, */
/* e.g., R_MATRIX_PACKAGE_VERSION <op> R_version(major, minor, patch) :    */

/* (version)_{10} = (major minor patch)_{256} */
#define R_MATRIX_PACKAGE_VERSION 67077
#define R_MATRIX_PACKAGE_MAJOR 1
#define R_MATRIX_PACKAGE_MINOR 6
#define R_MATRIX_PACKAGE_PATCH 5

#define R_MATRIX_ABI_VERSION 1

/* (version)_{10} = (major minor patch)_{256} */
#define R_MATRIX_SUITESPARSE_VERSION 330241
#define R_MATRIX_SUITESPARSE_MAJOR 5
#define R_MATRIX_SUITESPARSE_MINOR 10
#define R_MATRIX_SUITESPARSE_PATCH 1

#endif /* R_MATRIX_VERSION_H */
