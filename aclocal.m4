AC_DEFUN(OCTAVE_BLAS_LIBS, [
AC_ARG_WITH(fastblas,
  [  --with-fastblas         use specified fast BLAS],
  with_fastblas=$withval,
  with_fastblas=yes)
if test "$with_fastblas" = "no"; then
  BLAS_LIBS=" "
elif test "$with_fastblas" != "yes"; then
  # user specified a BLAS library to try on the command line
  AC_CHECK_LIB($with_fastblas, $dgemm_func, 
	       BLAS_LIBS="-l$with_fastblas", , $FLIBS)
fi

if test "x$BLAS_LIBS" = x; then
  # Checks for ATLAS BLAS library:
  AC_CHECK_LIB(atlas, ATL_xerbla, BLAS_LIBS="-latlas")
  if test "x$BLAS_LIBS" != x; then
    # check for other atlas libs:
    AC_CHECK_LIB(cblas, cblas_dgemm,BLAS_LIBS="-lcblas $BLAS_LIBS",,$BLAS_LIBS)
    AC_CHECK_LIB(f77blas, $dgemm_func, 
		 BLAS_LIBS="-lf77blas $BLAS_LIBS", , $BLAS_LIBS $FLIBS)
  fi
fi

if test "x$BLAS_LIBS" = x; then
  # BLAS in Alpha CXML library?
  AC_CHECK_LIB(cxml, $dgemm_func, BLAS_LIBS="-lcxml", , $FLIBS)
fi

if test "x$BLAS_LIBS" = x; then
  # BLAS in Alpha DXML library? (now called CXML, see above)
  AC_CHECK_LIB(dxml, $dgemm_func, BLAS_LIBS="-ldxml", , $FLIBS)
fi

if test "x$BLAS_LIBS" = x; then
  # Check for BLAS in Sun Performance library:
  AC_CHECK_LIB(sunmath, acosp, BLAS_LIBS="-lsunmath")
  AC_CHECK_LIB(sunperf, $dgemm_func, BLAS_LIBS="-xlic_lib=sunperf $BLAS_LIBS",
               , $BLAS_LIBS)
fi

if test "x$BLAS_LIBS" = x; then
  # Check for BLAS in SCSL and SGIMATH libraries (prefer SCSL):
  AC_CHECK_LIB(scs, $dgemm_func,
               BLAS_LIBS="-lscs", 
	       AC_CHECK_LIB(complib.sgimath, $dgemm_func,
			    BLAS_LIBS="-lcomplib.sgimath", , $FLIBS), $FLIBS)
fi

if test "x$BLAS_LIBS" = x; then
  # Checks for BLAS in IBM ESSL library.  We must also link
  # with -lblas in this case (ESSL does not include the full BLAS):
  AC_CHECK_LIB(blas, zherk, 
	       AC_CHECK_LIB(essl, $dgemm_func, 
			    BLAS_LIBS="-lessl -lblas", , $FLIBS), , $FLIBS)
fi

if test "x$BLAS_LIBS" = x; then
  # Finally, check for the generic BLAS library:
  AC_CHECK_LIB(blas, $dgemm_func, BLAS_LIBS="-lblas", , $FLIBS)
fi

if test "$with_fastblas" = "no"; then
  # Unset BLAS_LIBS so that we know below that nothing was found.
  BLAS_LIBS=""
fi

AC_SUBST(BLAS_LIBS)
])
