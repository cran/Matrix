AC_INIT([src/SVD.cc])
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`${R_HOME}/bin/R CMD config CC`
AC_PROG_CC
F77=`${R_HOME}/bin/R CMD config F77`
AC_PROG_F77
ACX_BLAS([], AC_MSG_ERROR([could not find your BLAS library], 1))
ACX_LAPACK([], AC_MSG_ERROR([could not find your LAPACK library], 1))
AC_OUTPUT([src/Makevars])
