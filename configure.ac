AC_PREREQ(2.50)
AC_INIT([src/Makevars.in])
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`${R_HOME}/bin/R CMD config CC`
CFLAGS=`${R_HOME}/bin/R CMD config CFLAGS`
AC_PROG_CC
F77=`${R_HOME}/bin/R CMD config F77`
if test -z "${F77}"; then
  AC_MSG_ERROR([R was configured without a FORTRAN compiler], 1)
fi
AC_PROG_F77
FLIBS=`${R_HOME}/bin/R CMD config FLIBS`
ACX_BLAS([], AC_MSG_ERROR([could not find your BLAS library], 1))
ACX_LAPACK([], AC_MSG_ERROR([could not find your LAPACK library], 1))
AC_OUTPUT([src/Makevars])
