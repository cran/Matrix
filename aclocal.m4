dnl @synopsis ACX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the BLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/).
dnl On success, it sets the BLAS_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with BLAS, you should link with:
dnl
dnl 	$BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl Many libraries are searched for, from ATLAS to CXML to ESSL.
dnl The user may also use --with-blas=<lib> in order to use some
dnl specific BLAS library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the BLAS library.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a BLAS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_BLAS.
dnl
dnl This macro requires autoconf 2.50 or later.
dnl
dnl @version $Id: aclocal.m4,v 1.3 2002/04/24 18:37:44 bates Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_BLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_blas_ok=no

AC_ARG_WITH(blas,
	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) acx_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
AC_F77_FUNC(sgemm)
AC_F77_FUNC(dgemm)

acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $acx_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
	AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($acx_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_CHECK_FUNC($sgemm, [acx_blas_ok=yes])
	LIBS="$save_LIBS"
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $sgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[acx_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(cxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(dxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $sgemm,
        			[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 acx_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(scs, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

AC_SUBST(BLAS_LIBS)

LIBS="$acx_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        acx_blas_ok=no
        $2
fi
])dnl ACX_BLAS
dnl @synopsis ACX_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the LAPACK
dnl linear-algebra interface (see http://www.netlib.org/lapack/).
dnl On success, it sets the LAPACK_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with LAPACK, you should link with:
dnl
dnl 	$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-lapack=<lib> in order to use some
dnl specific LAPACK library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the LAPACK and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a LAPACK
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_LAPACK.
dnl
dnl @version $Id: aclocal.m4,v 1.3 2002/04/24 18:37:44 bates Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_LAPACK], [
AC_REQUIRE([ACX_BLAS])
acx_lapack_ok=no

AC_ARG_WITH(lapack,
	[AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
case $with_lapack in
	yes | "") ;;
	no) acx_lapack_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
	*) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
AC_F77_FUNC(cheev)

# We cannot use LAPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
	acx_lapack_ok=noblas
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
	AC_MSG_CHECKING([for $cheev in $LAPACK_LIBS])
	AC_TRY_LINK_FUNC($cheev, [acx_lapack_ok=yes], [LAPACK_LIBS=""])
	AC_MSG_RESULT($acx_lapack_ok)
	LIBS="$save_LIBS"
	if test acx_lapack_ok = no; then
		LAPACK_LIBS=""
	fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
	AC_CHECK_FUNC($cheev, [acx_lapack_ok=yes])
	LIBS="$save_LIBS"
fi

# Generic LAPACK library?
if test $acx_lapack_ok = no; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_CHECK_LIB(lapack, $cheev,
		[acx_lapack_ok=yes; LAPACK_LIBS="-llapack"], [], [$FLIBS])
	LIBS="$save_LIBS"
fi

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        acx_lapack_ok=no
        $2
fi
])dnl ACX_LAPACK
#! /bin/sh
# Guess values for system-dependent variables and create Makefiles.
# Generated by Autoconf 2.52.
#
# Copyright 1992, 1993, 1994, 1995, 1996, 1998, 1999, 2000, 2001
# Free Software Foundation, Inc.
# This configure script is free software; the Free Software Foundation
# gives unlimited permission to copy, distribute and modify it.

# Avoid depending upon Character Ranges.
as_cr_letters='abcdefghijklmnopqrstuvwxyz'
as_cr_LETTERS='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
as_cr_Letters=$as_cr_letters$as_cr_LETTERS
as_cr_digits='0123456789'
as_cr_alnum=$as_cr_Letters$as_cr_digits

# Sed expression to map a string onto a valid variable name.
as_tr_sh="sed y%*+%pp%;s%[^_$as_cr_alnum]%_%g"

# Sed expression to map a string onto a valid CPP name.
as_tr_cpp="sed y%*$as_cr_letters%P$as_cr_LETTERS%;s%[^_$as_cr_alnum]%_%g"

# Be Bourne compatible
if test -n "${ZSH_VERSION+set}" && (emulate sh) >/dev/null 2>&1; then
  emulate sh
  NULLCMD=:
elif test -n "${BASH_VERSION+set}" && (set -o posix) >/dev/null 2>&1; then
  set -o posix
fi

# Name of the executable.
as_me=`echo "$0" |sed 's,.*[\\/],,'`

if expr a : '\(a\)' >/dev/null 2>&1; then
  as_expr=expr
else
  as_expr=false
fi

rm -f conf$$ conf$$.exe conf$$.file
echo >conf$$.file
if ln -s conf$$.file conf$$ 2>/dev/null; then
  # We could just check for DJGPP; but this test a) works b) is more generic
  # and c) will remain valid once DJGPP supports symlinks (DJGPP 2.04).
  if test -f conf$$.exe; then
    # Don't use ln at all; we don't have any links
    as_ln_s='cp -p'
  else
    as_ln_s='ln -s'
  fi
elif ln conf$$.file conf$$ 2>/dev/null; then
  as_ln_s=ln
else
  as_ln_s='cp -p'
fi
rm -f conf$$ conf$$.exe conf$$.file

as_executable_p="test -f"

# Support unset when possible.
if (FOO=FOO; unset FOO) >/dev/null 2>&1; then
  as_unset=unset
else
  as_unset=false
fi

# NLS nuisances.
$as_unset LANG || test "${LANG+set}" != set || { LANG=C; export LANG; }
$as_unset LC_ALL || test "${LC_ALL+set}" != set || { LC_ALL=C; export LC_ALL; }
$as_unset LC_TIME || test "${LC_TIME+set}" != set || { LC_TIME=C; export LC_TIME; }
$as_unset LC_CTYPE || test "${LC_CTYPE+set}" != set || { LC_CTYPE=C; export LC_CTYPE; }
$as_unset LANGUAGE || test "${LANGUAGE+set}" != set || { LANGUAGE=C; export LANGUAGE; }
$as_unset LC_COLLATE || test "${LC_COLLATE+set}" != set || { LC_COLLATE=C; export LC_COLLATE; }
$as_unset LC_NUMERIC || test "${LC_NUMERIC+set}" != set || { LC_NUMERIC=C; export LC_NUMERIC; }
$as_unset LC_MESSAGES || test "${LC_MESSAGES+set}" != set || { LC_MESSAGES=C; export LC_MESSAGES; }

# IFS
# We need space, tab and new line, in precisely that order.
as_nl='
'
IFS=" 	$as_nl"

# CDPATH.
$as_unset CDPATH || test "${CDPATH+set}" != set || { CDPATH=:; export CDPATH; }

# Name of the host.
# hostname on some systems (SVR3.2, Linux) returns a bogus exit status,
# so uname gets run too.
ac_hostname=`(hostname || uname -n) 2>/dev/null | sed 1q`

exec 6>&1

#
# Initializations.
#
ac_default_prefix=/usr/local
cross_compiling=no
subdirs=
MFLAGS= MAKEFLAGS=
SHELL=${CONFIG_SHELL-/bin/sh}

# Maximum number of lines to put in a shell here document.
# This variable seems obsolete.  It should probably be removed, and
# only ac_max_sed_lines should be used.
: ${ac_max_here_lines=38}

ac_unique_file="src/SVD.cc"

# Initialize some variables set by options.
ac_init_help=
ac_init_version=false
# The variables have the same names as the options, with
# dashes changed to underlines.
cache_file=/dev/null
exec_prefix=NONE
no_create=
no_recursion=
prefix=NONE
program_prefix=NONE
program_suffix=NONE
program_transform_name=s,x,x,
silent=
site=
srcdir=
verbose=
x_includes=NONE
x_libraries=NONE

# Installation directory options.
# These are left unexpanded so users can "make install exec_prefix=/foo"
# and all the variables that are supposed to be based on exec_prefix
# by default will actually change.
# Use braces instead of parens because sh, perl, etc. also accept them.
bindir='${exec_prefix}/bin'
sbindir='${exec_prefix}/sbin'
libexecdir='${exec_prefix}/libexec'
datadir='${prefix}/share'
sysconfdir='${prefix}/etc'
sharedstatedir='${prefix}/com'
localstatedir='${prefix}/var'
libdir='${exec_prefix}/lib'
includedir='${prefix}/include'
oldincludedir='/usr/include'
infodir='${prefix}/info'
mandir='${prefix}/man'

# Identity of this package.
PACKAGE_NAME=
PACKAGE_TARNAME=
PACKAGE_VERSION=
PACKAGE_STRING=
PACKAGE_BUGREPORT=

ac_prev=
for ac_option
do
  # If the previous option needs an argument, assign it.
  if test -n "$ac_prev"; then
    eval "$ac_prev=\$ac_option"
    ac_prev=
    continue
  fi

  ac_optarg=`expr "x$ac_option" : 'x[^=]*=\(.*\)'`

  # Accept the important Cygnus configure options, so we can diagnose typos.

  case $ac_option in

  -bindir | --bindir | --bindi | --bind | --bin | --bi)
    ac_prev=bindir ;;
  -bindir=* | --bindir=* | --bindi=* | --bind=* | --bin=* | --bi=*)
    bindir=$ac_optarg ;;

  -build | --build | --buil | --bui | --bu)
    ac_prev=build_alias ;;
  -build=* | --build=* | --buil=* | --bui=* | --bu=*)
    build_alias=$ac_optarg ;;

  -cache-file | --cache-file | --cache-fil | --cache-fi \
  | --cache-f | --cache- | --cache | --cach | --cac | --ca | --c)
    ac_prev=cache_file ;;
  -cache-file=* | --cache-file=* | --cache-fil=* | --cache-fi=* \
  | --cache-f=* | --cache-=* | --cache=* | --cach=* | --cac=* | --ca=* | --c=*)
    cache_file=$ac_optarg ;;

  --config-cache | -C)
    cache_file=config.cache ;;

  -datadir | --datadir | --datadi | --datad | --data | --dat | --da)
    ac_prev=datadir ;;
  -datadir=* | --datadir=* | --datadi=* | --datad=* | --data=* | --dat=* \
  | --da=*)
    datadir=$ac_optarg ;;

  -disable-* | --disable-*)
    ac_feature=`expr "x$ac_option" : 'x-*disable-\(.*\)'`
    # Reject names that are not valid shell variable names.
    expr "x$ac_feature" : ".*[^-_$as_cr_alnum]" >/dev/null &&
      { echo "$as_me: error: invalid feature name: $ac_feature" >&2
   { (exit 1); exit 1; }; }
    ac_feature=`echo $ac_feature | sed 's/-/_/g'`
    eval "enable_$ac_feature=no" ;;

  -enable-* | --enable-*)
    ac_feature=`expr "x$ac_option" : 'x-*enable-\([^=]*\)'`
    # Reject names that are not valid shell variable names.
    expr "x$ac_feature" : ".*[^-_$as_cr_alnum]" >/dev/null &&
      { echo "$as_me: error: invalid feature name: $ac_feature" >&2
   { (exit 1); exit 1; }; }
    ac_feature=`echo $ac_feature | sed 's/-/_/g'`
    case $ac_option in
      *=*) ac_optarg=`echo "$ac_optarg" | sed "s/'/'\\\\\\\\''/g"`;;
      *) ac_optarg=yes ;;
    esac
    eval "enable_$ac_feature='$ac_optarg'" ;;

  -exec-prefix | --exec_prefix | --exec-prefix | --exec-prefi \
  | --exec-pref | --exec-pre | --exec-pr | --exec-p | --exec- \
  | --exec | --exe | --ex)
    ac_prev=exec_prefix ;;
  -exec-prefix=* | --exec_prefix=* | --exec-prefix=* | --exec-prefi=* \
  | --exec-pref=* | --exec-pre=* | --exec-pr=* | --exec-p=* | --exec-=* \
  | --exec=* | --exe=* | --ex=*)
    exec_prefix=$ac_optarg ;;

  -gas | --gas | --ga | --g)
    # Obsolete; use --with-gas.
    with_gas=yes ;;

  -help | --help | --hel | --he | -h)
    ac_init_help=long ;;
  -help=r* | --help=r* | --hel=r* | --he=r* | -hr*)
    ac_init_help=recursive ;;
  -help=s* | --help=s* | --hel=s* | --he=s* | -hs*)
    ac_init_help=short ;;

  -host | --host | --hos | --ho)
    ac_prev=host_alias ;;
  -host=* | --host=* | --hos=* | --ho=*)
    host_alias=$ac_optarg ;;

  -includedir | --includedir | --includedi | --included | --include \
  | --includ | --inclu | --incl | --inc)
    ac_prev=includedir ;;
  -includedir=* | --includedir=* | --includedi=* | --included=* | --include=* \
  | --includ=* | --inclu=* | --incl=* | --inc=*)
    includedir=$ac_optarg ;;

  -infodir | --infodir | --infodi | --infod | --info | --inf)
    ac_prev=infodir ;;
  -infodir=* | --infodir=* | --infodi=* | --infod=* | --info=* | --inf=*)
    infodir=$ac_optarg ;;

  -libdir | --libdir | --libdi | --libd)
    ac_prev=libdir ;;
  -libdir=* | --libdir=* | --libdi=* | --libd=*)
    libdir=$ac_optarg ;;

  -libexecdir | --libexecdir | --libexecdi | --libexecd | --libexec \
  | --libexe | --libex | --libe)
    ac_prev=libexecdir ;;
  -libexecdir=* | --libexecdir=* | --libexecdi=* | --libexecd=* | --libexec=* \
  | --libexe=* | --libex=* | --libe=*)
    libexecdir=$ac_optarg ;;

  -localstatedir | --localstatedir | --localstatedi | --localstated \
  | --localstate | --localstat | --localsta | --localst \
  | --locals | --local | --loca | --loc | --lo)
    ac_prev=localstatedir ;;
  -localstatedir=* | --localstatedir=* | --localstatedi=* | --localstated=* \
  | --localstate=* | --localstat=* | --localsta=* | --localst=* \
  | --locals=* | --local=* | --loca=* | --loc=* | --lo=*)
    localstatedir=$ac_optarg ;;

  -mandir | --mandir | --mandi | --mand | --man | --ma | --m)
    ac_prev=mandir ;;
  -mandir=* | --mandir=* | --mandi=* | --mand=* | --man=* | --ma=* | --m=*)
    mandir=$ac_optarg ;;

  -nfp | --nfp | --nf)
    # Obsolete; use --without-fp.
    with_fp=no ;;

  -no-create | --no-create | --no-creat | --no-crea | --no-cre \
  | --no-cr | --no-c)
    no_create=yes ;;

  -no-recursion | --no-recursion | --no-recursio | --no-recursi \
  | --no-recurs | --no-recur | --no-recu | --no-rec | --no-re | --no-r)
    no_recursion=yes ;;

  -oldincludedir | --oldincludedir | --oldincludedi | --oldincluded \
  | --oldinclude | --oldinclud | --oldinclu | --oldincl | --oldinc \
  | --oldin | --oldi | --old | --ol | --o)
    ac_prev=oldincludedir ;;
  -oldincludedir=* | --oldincludedir=* | --oldincludedi=* | --oldincluded=* \
  | --oldinclude=* | --oldinclud=* | --oldinclu=* | --oldincl=* | --oldinc=* \
  | --oldin=* | --oldi=* | --old=* | --ol=* | --o=*)
    oldincludedir=$ac_optarg ;;

  -prefix | --prefix | --prefi | --pref | --pre | --pr | --p)
    ac_prev=prefix ;;
  -prefix=* | --prefix=* | --prefi=* | --pref=* | --pre=* | --pr=* | --p=*)
    prefix=$ac_optarg ;;

  -program-prefix | --program-prefix | --program-prefi | --program-pref \
  | --program-pre | --program-pr | --program-p)
    ac_prev=program_prefix ;;
  -program-prefix=* | --program-prefix=* | --program-prefi=* \
  | --program-pref=* | --program-pre=* | --program-pr=* | --program-p=*)
    program_prefix=$ac_optarg ;;

  -program-suffix | --program-suffix | --program-suffi | --program-suff \
  | --program-suf | --program-su | --program-s)
    ac_prev=program_suffix ;;
  -program-suffix=* | --program-suffix=* | --program-suffi=* \
  | --program-suff=* | --program-suf=* | --program-su=* | --program-s=*)
    program_suffix=$ac_optarg ;;

  -program-transform-name | --program-transform-name \
  | --program-transform-nam | --program-transform-na \
  | --program-transform-n | --program-transform- \
  | --program-transform | --program-transfor \
  | --program-transfo | --program-transf \
  | --program-trans | --program-tran \
  | --progr-tra | --program-tr | --program-t)
    ac_prev=program_transform_name ;;
  -program-transform-name=* | --program-transform-name=* \
  | --program-transform-nam=* | --program-transform-na=* \
  | --program-transform-n=* | --program-transform-=* \
  | --program-transform=* | --program-transfor=* \
  | --program-transfo=* | --program-transf=* \
  | --program-trans=* | --program-tran=* \
  | --progr-tra=* | --program-tr=* | --program-t=*)
    program_transform_name=$ac_optarg ;;

  -q | -quiet | --quiet | --quie | --qui | --qu | --q \
  | -silent | --silent | --silen | --sile | --sil)
    silent=yes ;;

  -sbindir | --sbindir | --sbindi | --sbind | --sbin | --sbi | --sb)
    ac_prev=sbindir ;;
  -sbindir=* | --sbindir=* | --sbindi=* | --sbind=* | --sbin=* \
  | --sbi=* | --sb=*)
    sbindir=$ac_optarg ;;

  -sharedstatedir | --sharedstatedir | --sharedstatedi \
  | --sharedstated | --sharedstate | --sharedstat | --sharedsta \
  | --sharedst | --shareds | --shared | --share | --shar \
  | --sha | --sh)
    ac_prev=sharedstatedir ;;
  -sharedstatedir=* | --sharedstatedir=* | --sharedstatedi=* \
  | --sharedstated=* | --sharedstate=* | --sharedstat=* | --sharedsta=* \
  | --sharedst=* | --shareds=* | --shared=* | --share=* | --shar=* \
  | --sha=* | --sh=*)
    sharedstatedir=$ac_optarg ;;

  -site | --site | --sit)
    ac_prev=site ;;
  -site=* | --site=* | --sit=*)
    site=$ac_optarg ;;

  -srcdir | --srcdir | --srcdi | --srcd | --src | --sr)
    ac_prev=srcdir ;;
  -srcdir=* | --srcdir=* | --srcdi=* | --srcd=* | --src=* | --sr=*)
    srcdir=$ac_optarg ;;

  -sysconfdir | --sysconfdir | --sysconfdi | --sysconfd | --sysconf \
  | --syscon | --sysco | --sysc | --sys | --sy)
    ac_prev=sysconfdir ;;
  -sysconfdir=* | --sysconfdir=* | --sysconfdi=* | --sysconfd=* | --sysconf=* \
  | --syscon=* | --sysco=* | --sysc=* | --sys=* | --sy=*)
    sysconfdir=$ac_optarg ;;

  -target | --target | --targe | --targ | --tar | --ta | --t)
    ac_prev=target_alias ;;
  -target=* | --target=* | --targe=* | --targ=* | --tar=* | --ta=* | --t=*)
    target_alias=$ac_optarg ;;

  -v | -verbose | --verbose | --verbos | --verbo | --verb)
    verbose=yes ;;

  -version | --version | --versio | --versi | --vers | -V)
    ac_init_version=: ;;

  -with-* | --with-*)
    ac_package=`expr "x$ac_option" : 'x-*with-\([^=]*\)'`
    # Reject names that are not valid shell variable names.
    expr "x$ac_package" : ".*[^-_$as_cr_alnum]" >/dev/null &&
      { echo "$as_me: error: invalid package name: $ac_package" >&2
   { (exit 1); exit 1; }; }
    ac_package=`echo $ac_package| sed 's/-/_/g'`
    case $ac_option in
      *=*) ac_optarg=`echo "$ac_optarg" | sed "s/'/'\\\\\\\\''/g"`;;
      *) ac_optarg=yes ;;
    esac
    eval "with_$ac_package='$ac_optarg'" ;;

  -without-* | --without-*)
    ac_package=`expr "x$ac_option" : 'x-*without-\(.*\)'`
    # Reject names that are not valid shell variable names.
    expr "x$ac_package" : ".*[^-_$as_cr_alnum]" >/dev/null &&
      { echo "$as_me: error: invalid package name: $ac_package" >&2
   { (exit 1); exit 1; }; }
    ac_package=`echo $ac_package | sed 's/-/_/g'`
    eval "with_$ac_package=no" ;;

  --x)
    # Obsolete; use --with-x.
    with_x=yes ;;

  -x-includes | --x-includes | --x-include | --x-includ | --x-inclu \
  | --x-incl | --x-inc | --x-in | --x-i)
    ac_prev=x_includes ;;
  -x-includes=* | --x-includes=* | --x-include=* | --x-includ=* | --x-inclu=* \
  | --x-incl=* | --x-inc=* | --x-in=* | --x-i=*)
    x_includes=$ac_optarg ;;

  -x-libraries | --x-libraries | --x-librarie | --x-librari \
  | --x-librar | --x-libra | --x-libr | --x-lib | --x-li | --x-l)
    ac_prev=x_libraries ;;
  -x-libraries=* | --x-libraries=* | --x-librarie=* | --x-librari=* \
  | --x-librar=* | --x-libra=* | --x-libr=* | --x-lib=* | --x-li=* | --x-l=*)
    x_libraries=$ac_optarg ;;

  -*) { echo "$as_me: error: unrecognized option: $ac_option
Try \`$0 --help' for more information." >&2
   { (exit 1); exit 1; }; }
    ;;

  *=*)
    ac_envvar=`expr "x$ac_option" : 'x\([^=]*\)='`
    # Reject names that are not valid shell variable names.
    expr "x$ac_envvar" : ".*[^_$as_cr_alnum]" >/dev/null &&
      { echo "$as_me: error: invalid variable name: $ac_envvar" >&2
   { (exit 1); exit 1; }; }
    ac_optarg=`echo "$ac_optarg" | sed "s/'/'\\\\\\\\''/g"`
    eval "$ac_envvar='$ac_optarg'"
    export $ac_envvar ;;

  *)
    # FIXME: should be removed in autoconf 3.0.
    echo "$as_me: WARNING: you should use --build, --host, --target" >&2
    expr "x$ac_option" : ".*[^-._$as_cr_alnum]" >/dev/null &&
      echo "$as_me: WARNING: invalid host type: $ac_option" >&2
    : ${build_alias=$ac_option} ${host_alias=$ac_option} ${target_alias=$ac_option}
    ;;

  esac
done

if test -n "$ac_prev"; then
  ac_option=--`echo $ac_prev | sed 's/_/-/g'`
  { echo "$as_me: error: missing argument to $ac_option" >&2
   { (exit 1); exit 1; }; }
fi

# Be sure to have absolute paths.
for ac_var in exec_prefix prefix
do
  eval ac_val=$`echo $ac_var`
  case $ac_val in
    [\\/$]* | ?:[\\/]* | NONE | '' ) ;;
    *)  { echo "$as_me: error: expected an absolute path for --$ac_var: $ac_val" >&2
   { (exit 1); exit 1; }; };;
  esac
done

# Be sure to have absolute paths.
for ac_var in bindir sbindir libexecdir datadir sysconfdir sharedstatedir \
              localstatedir libdir includedir oldincludedir infodir mandir
do
  eval ac_val=$`echo $ac_var`
  case $ac_val in
    [\\/$]* | ?:[\\/]* ) ;;
    *)  { echo "$as_me: error: expected an absolute path for --$ac_var: $ac_val" >&2
   { (exit 1); exit 1; }; };;
  esac
done

# There might be people who depend on the old broken behavior: `$host'
# used to hold the argument of --host etc.
build=$build_alias
host=$host_alias
target=$target_alias

# FIXME: should be removed in autoconf 3.0.
if test "x$host_alias" != x; then
  if test "x$build_alias" = x; then
    cross_compiling=maybe
    echo "$as_me: WARNING: If you wanted to set the --build type, don't use --host.
    If a cross compiler is detected then cross compile mode will be used." >&2
  elif test "x$build_alias" != "x$host_alias"; then
    cross_compiling=yes
  fi
fi

ac_tool_prefix=
test -n "$host_alias" && ac_tool_prefix=$host_alias-

test "$silent" = yes && exec 6>/dev/null

# Find the source files, if location was not specified.
if test -z "$srcdir"; then
  ac_srcdir_defaulted=yes
  # Try the directory containing this script, then its parent.
  ac_prog=$0
  ac_confdir=`echo "$ac_prog" | sed 's%[\\/][^\\/][^\\/]*$%%'`
  test "x$ac_confdir" = "x$ac_prog" && ac_confdir=.
  srcdir=$ac_confdir
  if test ! -r $srcdir/$ac_unique_file; then
    srcdir=..
  fi
else
  ac_srcdir_defaulted=no
fi
if test ! -r $srcdir/$ac_unique_file; then
  if test "$ac_srcdir_defaulted" = yes; then
    { echo "$as_me: error: cannot find sources in $ac_confdir or .." >&2
   { (exit 1); exit 1; }; }
  else
    { echo "$as_me: error: cannot find sources in $srcdir" >&2
   { (exit 1); exit 1; }; }
  fi
fi
srcdir=`echo "$srcdir" | sed 's%\([^\\/]\)[\\/]*$%\1%'`
ac_env_build_alias_set=${build_alias+set}
ac_env_build_alias_value=$build_alias
ac_cv_env_build_alias_set=${build_alias+set}
ac_cv_env_build_alias_value=$build_alias
ac_env_host_alias_set=${host_alias+set}
ac_env_host_alias_value=$host_alias
ac_cv_env_host_alias_set=${host_alias+set}
ac_cv_env_host_alias_value=$host_alias
ac_env_target_alias_set=${target_alias+set}
ac_env_target_alias_value=$target_alias
ac_cv_env_target_alias_set=${target_alias+set}
ac_cv_env_target_alias_value=$target_alias
ac_env_CC_set=${CC+set}
ac_env_CC_value=$CC
ac_cv_env_CC_set=${CC+set}
ac_cv_env_CC_value=$CC
ac_env_CFLAGS_set=${CFLAGS+set}
ac_env_CFLAGS_value=$CFLAGS
ac_cv_env_CFLAGS_set=${CFLAGS+set}
ac_cv_env_CFLAGS_value=$CFLAGS
ac_env_LDFLAGS_set=${LDFLAGS+set}
ac_env_LDFLAGS_value=$LDFLAGS
ac_cv_env_LDFLAGS_set=${LDFLAGS+set}
ac_cv_env_LDFLAGS_value=$LDFLAGS
ac_env_CPPFLAGS_set=${CPPFLAGS+set}
ac_env_CPPFLAGS_value=$CPPFLAGS
ac_cv_env_CPPFLAGS_set=${CPPFLAGS+set}
ac_cv_env_CPPFLAGS_value=$CPPFLAGS
ac_env_F77_set=${F77+set}
ac_env_F77_value=$F77
ac_cv_env_F77_set=${F77+set}
ac_cv_env_F77_value=$F77
ac_env_FFLAGS_set=${FFLAGS+set}
ac_env_FFLAGS_value=$FFLAGS
ac_cv_env_FFLAGS_set=${FFLAGS+set}
ac_cv_env_FFLAGS_value=$FFLAGS

#
# Report the --help message.
#
if test "$ac_init_help" = "long"; then
  # Omit some internal or obsolete options to make the list less imposing.
  # This message is too long to be a string in the A/UX 3.1 sh.
  cat <<EOF
\`configure' configures this package to adapt to many kinds of systems.

Usage: $0 [OPTION]... [VAR=VALUE]...

To assign environment variables (e.g., CC, CFLAGS...), specify them as
VAR=VALUE.  See below for descriptions of some of the useful variables.

Defaults for the options are specified in brackets.

Configuration:
  -h, --help              display this help and exit
      --help=short        display options specific to this package
      --help=recursive    display the short help of all the included packages
  -V, --version           display version information and exit
  -q, --quiet, --silent   do not print \`checking...' messages
      --cache-file=FILE   cache test results in FILE [disabled]
  -C, --config-cache      alias for \`--cache-file=config.cache'
  -n, --no-create         do not create output files
      --srcdir=DIR        find the sources in DIR [configure dir or \`..']

EOF

  cat <<EOF
Installation directories:
  --prefix=PREFIX         install architecture-independent files in PREFIX
                          [$ac_default_prefix]
  --exec-prefix=EPREFIX   install architecture-dependent files in EPREFIX
                          [PREFIX]

By default, \`make install' will install all the files in
\`$ac_default_prefix/bin', \`$ac_default_prefix/lib' etc.  You can specify
an installation prefix other than \`$ac_default_prefix' using \`--prefix',
for instance \`--prefix=\$HOME'.

For better control, use the options below.

Fine tuning of the installation directories:
  --bindir=DIR           user executables [EPREFIX/bin]
  --sbindir=DIR          system admin executables [EPREFIX/sbin]
  --libexecdir=DIR       program executables [EPREFIX/libexec]
  --datadir=DIR          read-only architecture-independent data [PREFIX/share]
  --sysconfdir=DIR       read-only single-machine data [PREFIX/etc]
  --sharedstatedir=DIR   modifiable architecture-independent data [PREFIX/com]
  --localstatedir=DIR    modifiable single-machine data [PREFIX/var]
  --libdir=DIR           object code libraries [EPREFIX/lib]
  --includedir=DIR       C header files [PREFIX/include]
  --oldincludedir=DIR    C header files for non-gcc [/usr/include]
  --infodir=DIR          info documentation [PREFIX/info]
  --mandir=DIR           man documentation [PREFIX/man]
EOF

  cat <<\EOF
EOF
fi

if test -n "$ac_init_help"; then

  cat <<\EOF

Optional Packages:
  --with-PACKAGE[=ARG]    use PACKAGE [ARG=yes]
  --without-PACKAGE       do not use PACKAGE (same as --with-PACKAGE=no)
  --with-blas=<lib>       use BLAS library <lib>
  --with-lapack=<lib>     use LAPACK library <lib>

Some influential environment variables:
  CC          C compiler command
  CFLAGS      C compiler flags
  LDFLAGS     linker flags, e.g. -L<lib dir> if you have libraries in a
              nonstandard directory <lib dir>
  CPPFLAGS    C/C++ preprocessor flags, e.g. -I<include dir> if you have
              headers in a nonstandard directory <include dir>
  F77         Fortran 77 compiler command
  FFLAGS      Fortran 77 compiler flags

Use these variables to override the choices made by `configure' or to help
it to find libraries and programs with nonstandard names/locations.

EOF
fi

if test "$ac_init_help" = "recursive"; then
  # If there are subdirs, report their specific --help.
  ac_popdir=`pwd`
  for ac_subdir in : $ac_subdirs_all; do test "x$ac_subdir" = x: && continue
    cd $ac_subdir
    # A "../" for each directory in /$ac_subdir.
    ac_dots=`echo $ac_subdir |
             sed 's,^\./,,;s,[^/]$,&/,;s,[^/]*/,../,g'`

    case $srcdir in
    .) # No --srcdir option.  We are building in place.
      ac_sub_srcdir=$srcdir ;;
    [\\/]* | ?:[\\/]* ) # Absolute path.
      ac_sub_srcdir=$srcdir/$ac_subdir ;;
    *) # Relative path.
      ac_sub_srcdir=$ac_dots$srcdir/$ac_subdir ;;
    esac

    # Check for guested configure; otherwise get Cygnus style configure.
    if test -f $ac_sub_srcdir/configure.gnu; then
      echo
      $SHELL $ac_sub_srcdir/configure.gnu  --help=recursive
    elif test -f $ac_sub_srcdir/configure; then
      echo
      $SHELL $ac_sub_srcdir/configure  --help=recursive
    elif test -f $ac_sub_srcdir/configure.ac ||
           test -f $ac_sub_srcdir/configure.in; then
      echo
      $ac_configure --help
    else
      echo "$as_me: WARNING: no configuration information is in $ac_subdir" >&2
    fi
    cd $ac_popdir
  done
fi

test -n "$ac_init_help" && exit 0
if $ac_init_version; then
  cat <<\EOF

Copyright 1992, 1993, 1994, 1995, 1996, 1998, 1999, 2000, 2001
Free Software Foundation, Inc.
This configure script is free software; the Free Software Foundation
gives unlimited permission to copy, distribute and modify it.
EOF
  exit 0
fi
exec 5>config.log
cat >&5 <<EOF
This file contains any messages produced by compilers while
running configure, to aid debugging if configure makes a mistake.

It was created by $as_me, which was
generated by GNU Autoconf 2.52.  Invocation command line was

  $ $0 $@

EOF
{
cat <<_ASUNAME
## ---------- ##
## Platform.  ##
## ---------- ##

hostname = `(hostname || uname -n) 2>/dev/null | sed 1q`
uname -m = `(uname -m) 2>/dev/null || echo unknown`
uname -r = `(uname -r) 2>/dev/null || echo unknown`
uname -s = `(uname -s) 2>/dev/null || echo unknown`
uname -v = `(uname -v) 2>/dev/null || echo unknown`

/usr/bin/uname -p = `(/usr/bin/uname -p) 2>/dev/null || echo unknown`
/bin/uname -X     = `(/bin/uname -X) 2>/dev/null     || echo unknown`

/bin/arch              = `(/bin/arch) 2>/dev/null              || echo unknown`
/usr/bin/arch -k       = `(/usr/bin/arch -k) 2>/dev/null       || echo unknown`
/usr/convex/getsysinfo = `(/usr/convex/getsysinfo) 2>/dev/null || echo unknown`
hostinfo               = `(hostinfo) 2>/dev/null               || echo unknown`
/bin/machine           = `(/bin/machine) 2>/dev/null           || echo unknown`
/usr/bin/oslevel       = `(/usr/bin/oslevel) 2>/dev/null       || echo unknown`
/bin/universe          = `(/bin/universe) 2>/dev/null          || echo unknown`

PATH = $PATH

_ASUNAME
} >&5

cat >&5 <<EOF
## ------------ ##
## Core tests.  ##
## ------------ ##

EOF

# Keep a trace of the command line.
# Strip out --no-create and --no-recursion so they do not pile up.
# Also quote any args containing shell meta-characters.
ac_configure_args=
ac_sep=
for ac_arg
do
  case $ac_arg in
  -no-create | --no-create | --no-creat | --no-crea | --no-cre \
  | --no-cr | --no-c) ;;
  -no-recursion | --no-recursion | --no-recursio | --no-recursi \
  | --no-recurs | --no-recur | --no-recu | --no-rec | --no-re | --no-r) ;;
  *" "*|*"	"*|*[\[\]\~\#\$\^\&\*\(\)\{\}\\\|\;\<\>\?\"\']*)
    ac_arg=`echo "$ac_arg" | sed "s/'/'\\\\\\\\''/g"`
    ac_configure_args="$ac_configure_args$ac_sep'$ac_arg'"
    ac_sep=" " ;;
  *) ac_configure_args="$ac_configure_args$ac_sep$ac_arg"
     ac_sep=" " ;;
  esac
  # Get rid of the leading space.
done

# When interrupted or exit'd, cleanup temporary files, and complete
# config.log.  We remove comments because anyway the quotes in there
# would cause problems or look ugly.
trap 'exit_status=$?
  # Save into config.log some information that might help in debugging.
  echo >&5
  echo "## ----------------- ##" >&5
  echo "## Cache variables.  ##" >&5
  echo "## ----------------- ##" >&5
  echo >&5
  # The following way of writing the cache mishandles newlines in values,
{
  (set) 2>&1 |
    case `(ac_space='"'"' '"'"'; set | grep ac_space) 2>&1` in
    *ac_space=\ *)
      sed -n \
        "s/'"'"'/'"'"'\\\\'"'"''"'"'/g;
    	  s/^\\([_$as_cr_alnum]*_cv_[_$as_cr_alnum]*\\)=\\(.*\\)/\\1='"'"'\\2'"'"'/p"
      ;;
    *)
      sed -n \
        "s/^\\([_$as_cr_alnum]*_cv_[_$as_cr_alnum]*\\)=\\(.*\\)/\\1=\\2/p"
      ;;
    esac;
} >&5
  sed "/^$/d" confdefs.h >conftest.log
  if test -s conftest.log; then
    echo >&5
    echo "## ------------ ##" >&5
    echo "## confdefs.h.  ##" >&5
    echo "## ------------ ##" >&5
    echo >&5
    cat conftest.log >&5
  fi
  (echo; echo) >&5
  test "$ac_signal" != 0 &&
    echo "$as_me: caught signal $ac_signal" >&5
  echo "$as_me: exit $exit_status" >&5
  rm -rf conftest* confdefs* core core.* *.core conf$$* $ac_clean_files &&
    exit $exit_status
     ' 0
for ac_signal in 1 2 13 15; do
  trap 'ac_signal='$ac_signal'; { (exit 1); exit 1; }' $ac_signal
done
ac_signal=0

# confdefs.h avoids OS command line length limits that DEFS can exceed.
rm -rf conftest* confdefs.h
# AIX cpp loses on an empty file, so make sure it contains at least a newline.
echo >confdefs.h

# Let the site file select an alternate cache file if it wants to.
# Prefer explicitly selected file to automatically selected ones.
if test -z "$CONFIG_SITE"; then
  if test "x$prefix" != xNONE; then
    CONFIG_SITE="$prefix/share/config.site $prefix/etc/config.site"
  else
    CONFIG_SITE="$ac_default_prefix/share/config.site $ac_default_prefix/etc/config.site"
  fi
fi
for ac_site_file in $CONFIG_SITE; do
  if test -r "$ac_site_file"; then
    { echo "$as_me:830: loading site script $ac_site_file" >&5
echo "$as_me: loading site script $ac_site_file" >&6;}
    cat "$ac_site_file" >&5
    . "$ac_site_file"
  fi
done

if test -r "$cache_file"; then
  # Some versions of bash will fail to source /dev/null (special
  # files actually), so we avoid doing that.
  if test -f "$cache_file"; then
    { echo "$as_me:841: loading cache $cache_file" >&5
echo "$as_me: loading cache $cache_file" >&6;}
    case $cache_file in
      [\\/]* | ?:[\\/]* ) . $cache_file;;
      *)                      . ./$cache_file;;
    esac
  fi
else
  { echo "$as_me:849: creating cache $cache_file" >&5
echo "$as_me: creating cache $cache_file" >&6;}
  >$cache_file
fi

# Check that the precious variables saved in the cache have kept the same
# value.
ac_cache_corrupted=false
for ac_var in `(set) 2>&1 |
               sed -n 's/^ac_env_\([a-zA-Z_0-9]*\)_set=.*/\1/p'`; do
  eval ac_old_set=\$ac_cv_env_${ac_var}_set
  eval ac_new_set=\$ac_env_${ac_var}_set
  eval ac_old_val="\$ac_cv_env_${ac_var}_value"
  eval ac_new_val="\$ac_env_${ac_var}_value"
  case $ac_old_set,$ac_new_set in
    set,)
      { echo "$as_me:865: error: \`$ac_var' was set to \`$ac_old_val' in the previous run" >&5
echo "$as_me: error: \`$ac_var' was set to \`$ac_old_val' in the previous run" >&2;}
      ac_cache_corrupted=: ;;
    ,set)
      { echo "$as_me:869: error: \`$ac_var' was not set in the previous run" >&5
echo "$as_me: error: \`$ac_var' was not set in the previous run" >&2;}
      ac_cache_corrupted=: ;;
    ,);;
    *)
      if test "x$ac_old_val" != "x$ac_new_val"; then
        { echo "$as_me:875: error: \`$ac_var' has changed since the previous run:" >&5
echo "$as_me: error: \`$ac_var' has changed since the previous run:" >&2;}
        { echo "$as_me:877:   former value:  $ac_old_val" >&5
echo "$as_me:   former value:  $ac_old_val" >&2;}
        { echo "$as_me:879:   current value: $ac_new_val" >&5
echo "$as_me:   current value: $ac_new_val" >&2;}
        ac_cache_corrupted=:
      fi;;
  esac
  # Pass precious variables to config.status.  It doesn't matter if
  # we pass some twice (in addition to the command line arguments).
  if test "$ac_new_set" = set; then
    case $ac_new_val in
    *" "*|*"	"*|*[\[\]\~\#\$\^\&\*\(\)\{\}\\\|\;\<\>\?\"\']*)
      ac_arg=$ac_var=`echo "$ac_new_val" | sed "s/'/'\\\\\\\\''/g"`
      ac_configure_args="$ac_configure_args '$ac_arg'"
      ;;
    *) ac_configure_args="$ac_configure_args $ac_var=$ac_new_val"
       ;;
    esac
  fi
done
if $ac_cache_corrupted; then
  { echo "$as_me:898: error: changes in the environment can compromise the build" >&5
echo "$as_me: error: changes in the environment can compromise the build" >&2;}
  { { echo "$as_me:900: error: run \`make distclean' and/or \`rm $cache_file' and start over" >&5
echo "$as_me: error: run \`make distclean' and/or \`rm $cache_file' and start over" >&2;}
   { (exit 1); exit 1; }; }
fi

ac_ext=c
ac_cpp='$CPP $CPPFLAGS'
ac_compile='$CC -c $CFLAGS $CPPFLAGS conftest.$ac_ext >&5'
ac_link='$CC -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_c_compiler_gnu

case `echo "testing\c"; echo 1,2,3`,`echo -n testing; echo 1,2,3` in
  *c*,-n*) ECHO_N= ECHO_C='
' ECHO_T='	' ;;
  *c*,*  ) ECHO_N=-n ECHO_C= ECHO_T= ;;
  *)      ECHO_N= ECHO_C='\c' ECHO_T= ;;
esac
echo "#! $SHELL" >conftest.sh
echo  "exit 0"   >>conftest.sh
chmod +x conftest.sh
if { (echo "$as_me:920: PATH=\".;.\"; conftest.sh") >&5
  (PATH=".;."; conftest.sh) 2>&5
  ac_status=$?
  echo "$as_me:923: \$? = $ac_status" >&5
  (exit $ac_status); }; then
  ac_path_separator=';'
else
  ac_path_separator=:
fi
PATH_SEPARATOR="$ac_path_separator"
rm -f conftest.sh

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`${R_HOME}/bin/R CMD config CC`
ac_ext=c
ac_cpp='$CPP $CPPFLAGS'
ac_compile='$CC -c $CFLAGS $CPPFLAGS conftest.$ac_ext >&5'
ac_link='$CC -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_c_compiler_gnu
if test -n "$ac_tool_prefix"; then
  # Extract the first word of "${ac_tool_prefix}gcc", so it can be a program name with args.
set dummy ${ac_tool_prefix}gcc; ac_word=$2
echo "$as_me:946: checking for $ac_word" >&5
echo $ECHO_N "checking for $ac_word... $ECHO_C" >&6
if test "${ac_cv_prog_CC+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  if test -n "$CC"; then
  ac_cv_prog_CC="$CC" # Let the user override the test.
else
  ac_save_IFS=$IFS; IFS=$ac_path_separator
ac_dummy="$PATH"
for ac_dir in $ac_dummy; do
  IFS=$ac_save_IFS
  test -z "$ac_dir" && ac_dir=.
  $as_executable_p "$ac_dir/$ac_word" || continue
ac_cv_prog_CC="${ac_tool_prefix}gcc"
echo "$as_me:961: found $ac_dir/$ac_word" >&5
break
done

fi
fi
CC=$ac_cv_prog_CC
if test -n "$CC"; then
  echo "$as_me:969: result: $CC" >&5
echo "${ECHO_T}$CC" >&6
else
  echo "$as_me:972: result: no" >&5
echo "${ECHO_T}no" >&6
fi

fi
if test -z "$ac_cv_prog_CC"; then
  ac_ct_CC=$CC
  # Extract the first word of "gcc", so it can be a program name with args.
set dummy gcc; ac_word=$2
echo "$as_me:981: checking for $ac_word" >&5
echo $ECHO_N "checking for $ac_word... $ECHO_C" >&6
if test "${ac_cv_prog_ac_ct_CC+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  if test -n "$ac_ct_CC"; then
  ac_cv_prog_ac_ct_CC="$ac_ct_CC" # Let the user override the test.
else
  ac_save_IFS=$IFS; IFS=$ac_path_separator
ac_dummy="$PATH"
for ac_dir in $ac_dummy; do
  IFS=$ac_save_IFS
  test -z "$ac_dir" && ac_dir=.
  $as_executable_p "$ac_dir/$ac_word" || continue
ac_cv_prog_ac_ct_CC="gcc"
echo "$as_me:996: found $ac_dir/$ac_word" >&5
break
done

fi
fi
ac_ct_CC=$ac_cv_prog_ac_ct_CC
if test -n "$ac_ct_CC"; then
  echo "$as_me:1004: result: $ac_ct_CC" >&5
echo "${ECHO_T}$ac_ct_CC" >&6
else
  echo "$as_me:1007: result: no" >&5
echo "${ECHO_T}no" >&6
fi

  CC=$ac_ct_CC
else
  CC="$ac_cv_prog_CC"
fi

if test -z "$CC"; then
  if test -n "$ac_tool_prefix"; then
  # Extract the first word of "${ac_tool_prefix}cc", so it can be a program name with args.
set dummy ${ac_tool_prefix}cc; ac_word=$2
echo "$as_me:1020: checking for $ac_word" >&5
echo $ECHO_N "checking for $ac_word... $ECHO_C" >&6
if test "${ac_cv_prog_CC+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  if test -n "$CC"; then
  ac_cv_prog_CC="$CC" # Let the user override the test.
else
  ac_save_IFS=$IFS; IFS=$ac_path_separator
ac_dummy="$PATH"
for ac_dir in $ac_dummy; do
  IFS=$ac_save_IFS
  test -z "$ac_dir" && ac_dir=.
  $as_executable_p "$ac_dir/$ac_word" || continue
ac_cv_prog_CC="${ac_tool_prefix}cc"
echo "$as_me:1035: found $ac_dir/$ac_word" >&5
break
done

fi
fi
CC=$ac_cv_prog_CC
if test -n "$CC"; then
  echo "$as_me:1043: result: $CC" >&5
echo "${ECHO_T}$CC" >&6
else
  echo "$as_me:1046: result: no" >&5
echo "${ECHO_T}no" >&6
fi

fi
if test -z "$ac_cv_prog_CC"; then
  ac_ct_CC=$CC
  # Extract the first word of "cc", so it can be a program name with args.
set dummy cc; ac_word=$2
echo "$as_me:1055: checking for $ac_word" >&5
echo $ECHO_N "checking for $ac_word... $ECHO_C" >&6
if test "${ac_cv_prog_ac_ct_CC+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  if test -n "$ac_ct_CC"; then
  ac_cv_prog_ac_ct_CC="$ac_ct_CC" # Let the user override the test.
else
  ac_save_IFS=$IFS; IFS=$ac_path_separator
ac_dummy="$PATH"
for ac_dir in $ac_dummy; do
  IFS=$ac_save_IFS
  test -z "$ac_dir" && ac_dir=.
  $as_executable_p "$ac_dir/$ac_word" || continue
ac_cv_prog_ac_ct_CC="cc"
echo "$as_me:1070: found $ac_dir/$ac_word" >&5
break
done

fi
fi
ac_ct_CC=$ac_cv_prog_ac_ct_CC
if test -n "$ac_ct_CC"; then
  echo "$as_me:1078: result: $ac_ct_CC" >&5
echo "${ECHO_T}$ac_ct_CC" >&6
else
  echo "$as_me:1081: result: no" >&5
echo "${ECHO_T}no" >&6
fi

  CC=$ac_ct_CC
else
  CC="$ac_cv_prog_CC"
fi

fi
if test -z "$CC"; then
  # Extract the first word of "cc", so it can be a program name with args.
set dummy cc; ac_word=$2
echo "$as_me:1094: checking for $ac_word" >&5
echo $ECHO_N "checking for $ac_word... $ECHO_C" >&6
if test "${ac_cv_prog_CC+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  if test -n "$CC"; then
  ac_cv_prog_CC="$CC" # Let the user override the test.
else
  ac_prog_rejected=no
  ac_save_IFS=$IFS; IFS=$ac_path_separator
ac_dummy="$PATH"
for ac_dir in $ac_dummy; do
  IFS=$ac_save_IFS
  test -z "$ac_dir" && ac_dir=.
  $as_executable_p "$ac_dir/$ac_word" || continue
if test "$ac_dir/$ac_word" = "/usr/ucb/cc"; then
  ac_prog_rejected=yes
  continue
fi
ac_cv_prog_CC="cc"
echo "$as_me:1114: found $ac_dir/$ac_word" >&5
break
done

if test $ac_prog_rejected = yes; then
  # We found a bogon in the path, so make sure we never use it.
  set dummy $ac_cv_prog_CC
  shift
  if test $# != 0; then
    # We chose a different compiler from the bogus one.
    # However, it has the same basename, so the bogon will be chosen
    # first if we set CC to just the basename; use the full file name.
    shift
    set dummy "$ac_dir/$ac_word" ${1+"$@"}
    shift
    ac_cv_prog_CC="$@"
  fi
fi
fi
fi
CC=$ac_cv_prog_CC
if test -n "$CC"; then
  echo "$as_me:1136: result: $CC" >&5
echo "${ECHO_T}$CC" >&6
else
  echo "$as_me:1139: result: no" >&5
echo "${ECHO_T}no" >&6
fi

fi
if test -z "$CC"; then
  if test -n "$ac_tool_prefix"; then
  for ac_prog in cl
  do
    # Extract the first word of "$ac_tool_prefix$ac_prog", so it can be a program name with args.
set dummy $ac_tool_prefix$ac_prog; ac_word=$2
echo "$as_me:1150: checking for $ac_word" >&5
echo $ECHO_N "checking for $ac_word... $ECHO_C" >&6
if test "${ac_cv_prog_CC+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  if test -n "$CC"; then
  ac_cv_prog_CC="$CC" # Let the user override the test.
else
  ac_save_IFS=$IFS; IFS=$ac_path_separator
ac_dummy="$PATH"
for ac_dir in $ac_dummy; do
  IFS=$ac_save_IFS
  test -z "$ac_dir" && ac_dir=.
  $as_executable_p "$ac_dir/$ac_word" || continue
ac_cv_prog_CC="$ac_tool_prefix$ac_prog"
echo "$as_me:1165: found $ac_dir/$ac_word" >&5
break
done

fi
fi
CC=$ac_cv_prog_CC
if test -n "$CC"; then
  echo "$as_me:1173: result: $CC" >&5
echo "${ECHO_T}$CC" >&6
else
  echo "$as_me:1176: result: no" >&5
echo "${ECHO_T}no" >&6
fi

    test -n "$CC" && break
  done
fi
if test -z "$CC"; then
  ac_ct_CC=$CC
  for ac_prog in cl
do
  # Extract the first word of "$ac_prog", so it can be a program name with args.
set dummy $ac_prog; ac_word=$2
echo "$as_me:1189: checking for $ac_word" >&5
echo $ECHO_N "checking for $ac_word... $ECHO_C" >&6
if test "${ac_cv_prog_ac_ct_CC+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  if test -n "$ac_ct_CC"; then
  ac_cv_prog_ac_ct_CC="$ac_ct_CC" # Let the user override the test.
else
  ac_save_IFS=$IFS; IFS=$ac_path_separator
ac_dummy="$PATH"
for ac_dir in $ac_dummy; do
  IFS=$ac_save_IFS
  test -z "$ac_dir" && ac_dir=.
  $as_executable_p "$ac_dir/$ac_word" || continue
ac_cv_prog_ac_ct_CC="$ac_prog"
echo "$as_me:1204: found $ac_dir/$ac_word" >&5
break
done

fi
fi
ac_ct_CC=$ac_cv_prog_ac_ct_CC
if test -n "$ac_ct_CC"; then
  echo "$as_me:1212: result: $ac_ct_CC" >&5
echo "${ECHO_T}$ac_ct_CC" >&6
else
  echo "$as_me:1215: result: no" >&5
echo "${ECHO_T}no" >&6
fi

  test -n "$ac_ct_CC" && break
done

  CC=$ac_ct_CC
fi

fi

test -z "$CC" && { { echo "$as_me:1227: error: no acceptable cc found in \$PATH" >&5
echo "$as_me: error: no acceptable cc found in \$PATH" >&2;}
   { (exit 1); exit 1; }; }

# Provide some information about the compiler.
echo "$as_me:1232:" \
     "checking for C compiler version" >&5
ac_compiler=`set X $ac_compile; echo $2`
{ (eval echo "$as_me:1235: \"$ac_compiler --version </dev/null >&5\"") >&5
  (eval $ac_compiler --version </dev/null >&5) 2>&5
  ac_status=$?
  echo "$as_me:1238: \$? = $ac_status" >&5
  (exit $ac_status); }
{ (eval echo "$as_me:1240: \"$ac_compiler -v </dev/null >&5\"") >&5
  (eval $ac_compiler -v </dev/null >&5) 2>&5
  ac_status=$?
  echo "$as_me:1243: \$? = $ac_status" >&5
  (exit $ac_status); }
{ (eval echo "$as_me:1245: \"$ac_compiler -V </dev/null >&5\"") >&5
  (eval $ac_compiler -V </dev/null >&5) 2>&5
  ac_status=$?
  echo "$as_me:1248: \$? = $ac_status" >&5
  (exit $ac_status); }

cat >conftest.$ac_ext <<_ACEOF
#line 1252 "configure"
#include "confdefs.h"

int
main ()
{

  ;
  return 0;
}
_ACEOF
ac_clean_files_save=$ac_clean_files
ac_clean_files="$ac_clean_files a.out a.exe"
# Try to create an executable without -o first, disregard a.out.
# It will help us diagnose broken compilers, and finding out an intuition
# of exeext.
echo "$as_me:1268: checking for C compiler default output" >&5
echo $ECHO_N "checking for C compiler default output... $ECHO_C" >&6
ac_link_default=`echo "$ac_link" | sed 's/ -o *conftest[^ ]*//'`
if { (eval echo "$as_me:1271: \"$ac_link_default\"") >&5
  (eval $ac_link_default) 2>&5
  ac_status=$?
  echo "$as_me:1274: \$? = $ac_status" >&5
  (exit $ac_status); }; then
  # Find the output, starting from the most likely.  This scheme is
# not robust to junk in `.', hence go to wildcards (a.*) only as a last
# resort.
for ac_file in `ls a.exe conftest.exe 2>/dev/null;
                ls a.out conftest 2>/dev/null;
                ls a.* conftest.* 2>/dev/null`; do
  case $ac_file in
    *.$ac_ext | *.o | *.obj | *.xcoff | *.tds | *.d | *.pdb ) ;;
    a.out ) # We found the default executable, but exeext='' is most
            # certainly right.
            break;;
    *.* ) ac_cv_exeext=`expr "$ac_file" : '[^.]*\(\..*\)'`
          # FIXME: I believe we export ac_cv_exeext for Libtool --akim.
          export ac_cv_exeext
          break;;
    * ) break;;
  esac
done
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
{ { echo "$as_me:1297: error: C compiler cannot create executables" >&5
echo "$as_me: error: C compiler cannot create executables" >&2;}
   { (exit 77); exit 77; }; }
fi

ac_exeext=$ac_cv_exeext
echo "$as_me:1303: result: $ac_file" >&5
echo "${ECHO_T}$ac_file" >&6

# Check the compiler produces executables we can run.  If not, either
# the compiler is broken, or we cross compile.
echo "$as_me:1308: checking whether the C compiler works" >&5
echo $ECHO_N "checking whether the C compiler works... $ECHO_C" >&6
# FIXME: These cross compiler hacks should be removed for Autoconf 3.0
# If not cross compiling, check that we can run a simple program.
if test "$cross_compiling" != yes; then
  if { ac_try='./$ac_file'
  { (eval echo "$as_me:1314: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:1317: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
    cross_compiling=no
  else
    if test "$cross_compiling" = maybe; then
	cross_compiling=yes
    else
	{ { echo "$as_me:1324: error: cannot run C compiled programs.
If you meant to cross compile, use \`--host'." >&5
echo "$as_me: error: cannot run C compiled programs.
If you meant to cross compile, use \`--host'." >&2;}
   { (exit 1); exit 1; }; }
    fi
  fi
fi
echo "$as_me:1332: result: yes" >&5
echo "${ECHO_T}yes" >&6

rm -f a.out a.exe conftest$ac_cv_exeext
ac_clean_files=$ac_clean_files_save
# Check the compiler produces executables we can run.  If not, either
# the compiler is broken, or we cross compile.
echo "$as_me:1339: checking whether we are cross compiling" >&5
echo $ECHO_N "checking whether we are cross compiling... $ECHO_C" >&6
echo "$as_me:1341: result: $cross_compiling" >&5
echo "${ECHO_T}$cross_compiling" >&6

echo "$as_me:1344: checking for executable suffix" >&5
echo $ECHO_N "checking for executable suffix... $ECHO_C" >&6
if { (eval echo "$as_me:1346: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:1349: \$? = $ac_status" >&5
  (exit $ac_status); }; then
  # If both `conftest.exe' and `conftest' are `present' (well, observable)
# catch `conftest.exe'.  For instance with Cygwin, `ls conftest' will
# work properly (i.e., refer to `conftest.exe'), while it won't with
# `rm'.
for ac_file in `(ls conftest.exe; ls conftest; ls conftest.*) 2>/dev/null`; do
  case $ac_file in
    *.$ac_ext | *.o | *.obj | *.xcoff | *.tds | *.d | *.pdb ) ;;
    *.* ) ac_cv_exeext=`expr "$ac_file" : '[^.]*\(\..*\)'`
          export ac_cv_exeext
          break;;
    * ) break;;
  esac
done
else
  { { echo "$as_me:1365: error: cannot compute EXEEXT: cannot compile and link" >&5
echo "$as_me: error: cannot compute EXEEXT: cannot compile and link" >&2;}
   { (exit 1); exit 1; }; }
fi

rm -f conftest$ac_cv_exeext
echo "$as_me:1371: result: $ac_cv_exeext" >&5
echo "${ECHO_T}$ac_cv_exeext" >&6

rm -f conftest.$ac_ext
EXEEXT=$ac_cv_exeext
ac_exeext=$EXEEXT
echo "$as_me:1377: checking for object suffix" >&5
echo $ECHO_N "checking for object suffix... $ECHO_C" >&6
if test "${ac_cv_objext+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  cat >conftest.$ac_ext <<_ACEOF
#line 1383 "configure"
#include "confdefs.h"

int
main ()
{

  ;
  return 0;
}
_ACEOF
rm -f conftest.o conftest.obj
if { (eval echo "$as_me:1395: \"$ac_compile\"") >&5
  (eval $ac_compile) 2>&5
  ac_status=$?
  echo "$as_me:1398: \$? = $ac_status" >&5
  (exit $ac_status); }; then
  for ac_file in `(ls conftest.o conftest.obj; ls conftest.*) 2>/dev/null`; do
  case $ac_file in
    *.$ac_ext | *.xcoff | *.tds | *.d | *.pdb ) ;;
    *) ac_cv_objext=`expr "$ac_file" : '.*\.\(.*\)'`
       break;;
  esac
done
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
{ { echo "$as_me:1410: error: cannot compute OBJEXT: cannot compile" >&5
echo "$as_me: error: cannot compute OBJEXT: cannot compile" >&2;}
   { (exit 1); exit 1; }; }
fi

rm -f conftest.$ac_cv_objext conftest.$ac_ext
fi
echo "$as_me:1417: result: $ac_cv_objext" >&5
echo "${ECHO_T}$ac_cv_objext" >&6
OBJEXT=$ac_cv_objext
ac_objext=$OBJEXT
echo "$as_me:1421: checking whether we are using the GNU C compiler" >&5
echo $ECHO_N "checking whether we are using the GNU C compiler... $ECHO_C" >&6
if test "${ac_cv_c_compiler_gnu+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  cat >conftest.$ac_ext <<_ACEOF
#line 1427 "configure"
#include "confdefs.h"

int
main ()
{
#ifndef __GNUC__
       choke me
#endif

  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext
if { (eval echo "$as_me:1442: \"$ac_compile\"") >&5
  (eval $ac_compile) 2>&5
  ac_status=$?
  echo "$as_me:1445: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest.$ac_objext'
  { (eval echo "$as_me:1448: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:1451: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  ac_compiler_gnu=yes
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
ac_compiler_gnu=no
fi
rm -f conftest.$ac_objext conftest.$ac_ext
ac_cv_c_compiler_gnu=$ac_compiler_gnu

fi
echo "$as_me:1463: result: $ac_cv_c_compiler_gnu" >&5
echo "${ECHO_T}$ac_cv_c_compiler_gnu" >&6
GCC=`test $ac_compiler_gnu = yes && echo yes`
ac_test_CFLAGS=${CFLAGS+set}
ac_save_CFLAGS=$CFLAGS
CFLAGS="-g"
echo "$as_me:1469: checking whether $CC accepts -g" >&5
echo $ECHO_N "checking whether $CC accepts -g... $ECHO_C" >&6
if test "${ac_cv_prog_cc_g+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  cat >conftest.$ac_ext <<_ACEOF
#line 1475 "configure"
#include "confdefs.h"

int
main ()
{

  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext
if { (eval echo "$as_me:1487: \"$ac_compile\"") >&5
  (eval $ac_compile) 2>&5
  ac_status=$?
  echo "$as_me:1490: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest.$ac_objext'
  { (eval echo "$as_me:1493: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:1496: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  ac_cv_prog_cc_g=yes
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
ac_cv_prog_cc_g=no
fi
rm -f conftest.$ac_objext conftest.$ac_ext
fi
echo "$as_me:1506: result: $ac_cv_prog_cc_g" >&5
echo "${ECHO_T}$ac_cv_prog_cc_g" >&6
if test "$ac_test_CFLAGS" = set; then
  CFLAGS=$ac_save_CFLAGS
elif test $ac_cv_prog_cc_g = yes; then
  if test "$GCC" = yes; then
    CFLAGS="-g -O2"
  else
    CFLAGS="-g"
  fi
else
  if test "$GCC" = yes; then
    CFLAGS="-O2"
  else
    CFLAGS=
  fi
fi
# Some people use a C++ compiler to compile C.  Since we use `exit',
# in C++ we need to declare it.  In case someone uses the same compiler
# for both compiling C and C++ we need to have the C++ compiler decide
# the declaration of exit, since it's the most demanding environment.
cat >conftest.$ac_ext <<_ACEOF
#ifndef __cplusplus
  choke me
#endif
_ACEOF
rm -f conftest.$ac_objext
if { (eval echo "$as_me:1533: \"$ac_compile\"") >&5
  (eval $ac_compile) 2>&5
  ac_status=$?
  echo "$as_me:1536: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest.$ac_objext'
  { (eval echo "$as_me:1539: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:1542: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  for ac_declaration in \
   ''\
   '#include <stdlib.h>' \
   'extern "C" void std::exit (int) throw (); using std::exit;' \
   'extern "C" void std::exit (int); using std::exit;' \
   'extern "C" void exit (int) throw ();' \
   'extern "C" void exit (int);' \
   'void exit (int);'
do
  cat >conftest.$ac_ext <<_ACEOF
#line 1554 "configure"
#include "confdefs.h"
#include <stdlib.h>
$ac_declaration
int
main ()
{
exit (42);
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext
if { (eval echo "$as_me:1567: \"$ac_compile\"") >&5
  (eval $ac_compile) 2>&5
  ac_status=$?
  echo "$as_me:1570: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest.$ac_objext'
  { (eval echo "$as_me:1573: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:1576: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  :
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
continue
fi
rm -f conftest.$ac_objext conftest.$ac_ext
  cat >conftest.$ac_ext <<_ACEOF
#line 1586 "configure"
#include "confdefs.h"
$ac_declaration
int
main ()
{
exit (42);
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext
if { (eval echo "$as_me:1598: \"$ac_compile\"") >&5
  (eval $ac_compile) 2>&5
  ac_status=$?
  echo "$as_me:1601: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest.$ac_objext'
  { (eval echo "$as_me:1604: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:1607: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  break
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
fi
rm -f conftest.$ac_objext conftest.$ac_ext
done
rm -f conftest*
if test -n "$ac_declaration"; then
  echo '#ifdef __cplusplus' >>confdefs.h
  echo $ac_declaration      >>confdefs.h
  echo '#endif'             >>confdefs.h
fi

else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
fi
rm -f conftest.$ac_objext conftest.$ac_ext
ac_ext=c
ac_cpp='$CPP $CPPFLAGS'
ac_compile='$CC -c $CFLAGS $CPPFLAGS conftest.$ac_ext >&5'
ac_link='$CC -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_c_compiler_gnu

F77=`${R_HOME}/bin/R CMD config F77`
ac_ext=f
ac_compile='$F77 -c $FFLAGS conftest.$ac_ext >&5'
ac_link='$F77 -o conftest$ac_exeext $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_f77_compiler_gnu
if test -n "$ac_tool_prefix"; then
  for ac_prog in g77 f77 xlf cf77 cft77 frt pgf77 fl32 af77 fort77 f90 xlf90 pgf90 epcf90 f95 fort xlf95 lf95 g95 fc
  do
    # Extract the first word of "$ac_tool_prefix$ac_prog", so it can be a program name with args.
set dummy $ac_tool_prefix$ac_prog; ac_word=$2
echo "$as_me:1644: checking for $ac_word" >&5
echo $ECHO_N "checking for $ac_word... $ECHO_C" >&6
if test "${ac_cv_prog_F77+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  if test -n "$F77"; then
  ac_cv_prog_F77="$F77" # Let the user override the test.
else
  ac_save_IFS=$IFS; IFS=$ac_path_separator
ac_dummy="$PATH"
for ac_dir in $ac_dummy; do
  IFS=$ac_save_IFS
  test -z "$ac_dir" && ac_dir=.
  $as_executable_p "$ac_dir/$ac_word" || continue
ac_cv_prog_F77="$ac_tool_prefix$ac_prog"
echo "$as_me:1659: found $ac_dir/$ac_word" >&5
break
done

fi
fi
F77=$ac_cv_prog_F77
if test -n "$F77"; then
  echo "$as_me:1667: result: $F77" >&5
echo "${ECHO_T}$F77" >&6
else
  echo "$as_me:1670: result: no" >&5
echo "${ECHO_T}no" >&6
fi

    test -n "$F77" && break
  done
fi
if test -z "$F77"; then
  ac_ct_F77=$F77
  for ac_prog in g77 f77 xlf cf77 cft77 frt pgf77 fl32 af77 fort77 f90 xlf90 pgf90 epcf90 f95 fort xlf95 lf95 g95 fc
do
  # Extract the first word of "$ac_prog", so it can be a program name with args.
set dummy $ac_prog; ac_word=$2
echo "$as_me:1683: checking for $ac_word" >&5
echo $ECHO_N "checking for $ac_word... $ECHO_C" >&6
if test "${ac_cv_prog_ac_ct_F77+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  if test -n "$ac_ct_F77"; then
  ac_cv_prog_ac_ct_F77="$ac_ct_F77" # Let the user override the test.
else
  ac_save_IFS=$IFS; IFS=$ac_path_separator
ac_dummy="$PATH"
for ac_dir in $ac_dummy; do
  IFS=$ac_save_IFS
  test -z "$ac_dir" && ac_dir=.
  $as_executable_p "$ac_dir/$ac_word" || continue
ac_cv_prog_ac_ct_F77="$ac_prog"
echo "$as_me:1698: found $ac_dir/$ac_word" >&5
break
done

fi
fi
ac_ct_F77=$ac_cv_prog_ac_ct_F77
if test -n "$ac_ct_F77"; then
  echo "$as_me:1706: result: $ac_ct_F77" >&5
echo "${ECHO_T}$ac_ct_F77" >&6
else
  echo "$as_me:1709: result: no" >&5
echo "${ECHO_T}no" >&6
fi

  test -n "$ac_ct_F77" && break
done

  F77=$ac_ct_F77
fi

# Provide some information about the compiler.
echo "$as_me:1720:" \
     "checking for Fortran 77 compiler version" >&5
ac_compiler=`set X $ac_compile; echo $2`
{ (eval echo "$as_me:1723: \"$ac_compiler --version </dev/null >&5\"") >&5
  (eval $ac_compiler --version </dev/null >&5) 2>&5
  ac_status=$?
  echo "$as_me:1726: \$? = $ac_status" >&5
  (exit $ac_status); }
{ (eval echo "$as_me:1728: \"$ac_compiler -v </dev/null >&5\"") >&5
  (eval $ac_compiler -v </dev/null >&5) 2>&5
  ac_status=$?
  echo "$as_me:1731: \$? = $ac_status" >&5
  (exit $ac_status); }
{ (eval echo "$as_me:1733: \"$ac_compiler -V </dev/null >&5\"") >&5
  (eval $ac_compiler -V </dev/null >&5) 2>&5
  ac_status=$?
  echo "$as_me:1736: \$? = $ac_status" >&5
  (exit $ac_status); }

# If we don't use `.F' as extension, the preprocessor is not run on the
# input file.
ac_save_ext=$ac_ext
ac_ext=F
echo "$as_me:1743: checking whether we are using the GNU Fortran 77 compiler" >&5
echo $ECHO_N "checking whether we are using the GNU Fortran 77 compiler... $ECHO_C" >&6
if test "${ac_cv_f77_compiler_gnu+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  cat >conftest.$ac_ext <<_ACEOF
      program main
#ifndef __GNUC__
       choke me
#endif

      end
_ACEOF
rm -f conftest.$ac_objext
if { (eval echo "$as_me:1757: \"$ac_compile\"") >&5
  (eval $ac_compile) 2>&5
  ac_status=$?
  echo "$as_me:1760: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest.$ac_objext'
  { (eval echo "$as_me:1763: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:1766: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  ac_compiler_gnu=yes
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
ac_compiler_gnu=no
fi
rm -f conftest.$ac_objext conftest.$ac_ext
ac_cv_f77_compiler_gnu=$ac_compiler_gnu

fi
echo "$as_me:1778: result: $ac_cv_f77_compiler_gnu" >&5
echo "${ECHO_T}$ac_cv_f77_compiler_gnu" >&6
ac_ext=$ac_save_ext
G77=`test $ac_compiler_gnu = yes && echo yes`
ac_test_FFLAGS=${FFLAGS+set}
ac_save_FFLAGS=$FFLAGS
FFLAGS=
echo "$as_me:1785: checking whether $F77 accepts -g" >&5
echo $ECHO_N "checking whether $F77 accepts -g... $ECHO_C" >&6
if test "${ac_cv_prog_f77_g+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  FFLAGS=-g
cat >conftest.$ac_ext <<_ACEOF
      program main

      end
_ACEOF
rm -f conftest.$ac_objext
if { (eval echo "$as_me:1797: \"$ac_compile\"") >&5
  (eval $ac_compile) 2>&5
  ac_status=$?
  echo "$as_me:1800: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest.$ac_objext'
  { (eval echo "$as_me:1803: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:1806: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  ac_cv_prog_f77_g=yes
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
ac_cv_prog_f77_g=no
fi
rm -f conftest.$ac_objext conftest.$ac_ext

fi
echo "$as_me:1817: result: $ac_cv_prog_f77_g" >&5
echo "${ECHO_T}$ac_cv_prog_f77_g" >&6
if test "$ac_test_FFLAGS" = set; then
  FFLAGS=$ac_save_FFLAGS
elif test $ac_cv_prog_f77_g = yes; then
  if test "$G77" = yes; then
    FFLAGS="-g -O2"
  else
    FFLAGS="-g"
  fi
else
  if test "$G77" = yes; then
    FFLAGS="-O2"
  else
    FFLAGS=
  fi
fi
ac_ext=c
ac_cpp='$CPP $CPPFLAGS'
ac_compile='$CC -c $CFLAGS $CPPFLAGS conftest.$ac_ext >&5'
ac_link='$CC -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_c_compiler_gnu

ac_ext=f
ac_compile='$F77 -c $FFLAGS conftest.$ac_ext >&5'
ac_link='$F77 -o conftest$ac_exeext $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_f77_compiler_gnu
echo "$as_me:1844: checking how to get verbose linking output from $F77" >&5
echo $ECHO_N "checking how to get verbose linking output from $F77... $ECHO_C" >&6
if test "${ac_cv_prog_f77_v+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else

cat >conftest.$ac_ext <<_ACEOF
      program main

      end
_ACEOF
rm -f conftest.$ac_objext
if { (eval echo "$as_me:1856: \"$ac_compile\"") >&5
  (eval $ac_compile) 2>&5
  ac_status=$?
  echo "$as_me:1859: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest.$ac_objext'
  { (eval echo "$as_me:1862: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:1865: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  ac_cv_prog_f77_v=
# Try some options frequently used verbose output
for ac_verb in -v -verbose --verbose -V -\#\#\#; do
  ac_ext=f
ac_compile='$F77 -c $FFLAGS conftest.$ac_ext >&5'
ac_link='$F77 -o conftest$ac_exeext $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_f77_compiler_gnu

cat >conftest.$ac_ext <<_ACEOF
      program main

      end
_ACEOF

# Compile and link our simple test program by passing a flag (argument
# 1 to this macro) to the Fortran 77 compiler in order to get
# "verbose" output that we can then parse for the Fortran 77 linker
# flags.
ac_save_FFLAGS=$FFLAGS
FFLAGS="$FFLAGS $ac_verb"
(eval echo $as_me:1887: \"$ac_link\") >&5
ac_f77_v_output=`eval $ac_link 5>&1 2>&1 | grep -v 'Driving:'`
echo "$ac_f77_v_output" >&5
FFLAGS=$ac_save_FFLAGS

rm -f conftest*
ac_ext=f
ac_compile='$F77 -c $FFLAGS conftest.$ac_ext >&5'
ac_link='$F77 -o conftest$ac_exeext $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_f77_compiler_gnu

# If we are using xlf then replace all the commas with spaces.
if echo $ac_f77_v_output | grep xlfentry >/dev/null 2>&1; then
  ac_f77_v_output=`echo $ac_f77_v_output | sed 's/,/ /g'`
fi

# If we are using Cray Fortran then delete quotes.
# Use "\"" instead of '"' for font-lock-mode.
# FIXME: a more general fix for quoted arguments with spaces?
if echo $ac_f77_v_output | grep cft90 >/dev/null 2>&1; then
  ac_f77_v_output=`echo $ac_f77_v_output | sed "s/\"//g"`
fi
  # look for -l* and *.a constructs in the output
  for ac_arg in $ac_f77_v_output; do
     case $ac_arg in
        [\\/]*.a | ?:[\\/]*.a | -[lLRu]*)
          ac_cv_prog_f77_v=$ac_verb
          break 2 ;;
     esac
  done
done
if test -z "$ac_cv_prog_f77_v"; then
   { echo "$as_me:1919: WARNING: cannot determine how to obtain linking information from $F77" >&5
echo "$as_me: WARNING: cannot determine how to obtain linking information from $F77" >&2;}
fi
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
{ echo "$as_me:1925: WARNING: compilation failed" >&5
echo "$as_me: WARNING: compilation failed" >&2;}
fi
rm -f conftest.$ac_objext conftest.$ac_ext

fi
echo "$as_me:1931: result: $ac_cv_prog_f77_v" >&5
echo "${ECHO_T}$ac_cv_prog_f77_v" >&6
echo "$as_me:1933: checking for Fortran 77 libraries" >&5
echo $ECHO_N "checking for Fortran 77 libraries... $ECHO_C" >&6
if test "${ac_cv_flibs+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  if test "x$FLIBS" != "x"; then
  ac_cv_flibs="$FLIBS" # Let the user override the test.
else

ac_ext=f
ac_compile='$F77 -c $FFLAGS conftest.$ac_ext >&5'
ac_link='$F77 -o conftest$ac_exeext $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_f77_compiler_gnu

cat >conftest.$ac_ext <<_ACEOF
      program main

      end
_ACEOF

# Compile and link our simple test program by passing a flag (argument
# 1 to this macro) to the Fortran 77 compiler in order to get
# "verbose" output that we can then parse for the Fortran 77 linker
# flags.
ac_save_FFLAGS=$FFLAGS
FFLAGS="$FFLAGS $ac_cv_prog_f77_v"
(eval echo $as_me:1959: \"$ac_link\") >&5
ac_f77_v_output=`eval $ac_link 5>&1 2>&1 | grep -v 'Driving:'`
echo "$ac_f77_v_output" >&5
FFLAGS=$ac_save_FFLAGS

rm -f conftest*
ac_ext=f
ac_compile='$F77 -c $FFLAGS conftest.$ac_ext >&5'
ac_link='$F77 -o conftest$ac_exeext $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_f77_compiler_gnu

# If we are using xlf then replace all the commas with spaces.
if echo $ac_f77_v_output | grep xlfentry >/dev/null 2>&1; then
  ac_f77_v_output=`echo $ac_f77_v_output | sed 's/,/ /g'`
fi

# If we are using Cray Fortran then delete quotes.
# Use "\"" instead of '"' for font-lock-mode.
# FIXME: a more general fix for quoted arguments with spaces?
if echo $ac_f77_v_output | grep cft90 >/dev/null 2>&1; then
  ac_f77_v_output=`echo $ac_f77_v_output | sed "s/\"//g"`
fi

ac_cv_flibs=

# Save positional arguments (if any)
ac_save_positional="$@"

set X $ac_f77_v_output
while test $# != 1; do
  shift
  ac_arg=$1
  case $ac_arg in
        [\\/]*.a | ?:[\\/]*.a)

  ac_exists=false
  for ac_i in $ac_cv_flibs; do
    if test x"$ac_arg" = x"$ac_i"; then
      ac_exists=true
      break
    fi
  done

  if test x"$ac_exists" = xtrue; then
  :
else
  ac_cv_flibs="$ac_cv_flibs $ac_arg"
fi

          ;;
        -bI:*)

  ac_exists=false
  for ac_i in $ac_cv_flibs; do
    if test x"$ac_arg" = x"$ac_i"; then
      ac_exists=true
      break
    fi
  done

  if test x"$ac_exists" = xtrue; then
  :
else
  if test "$ac_compiler_gnu" = yes; then
  for ac_link_opt in $ac_arg; do
    ac_cv_flibs="$ac_cv_flibs -Xlinker $ac_link_opt"
  done
else
  ac_cv_flibs="$ac_cv_flibs $ac_arg"
fi
fi

          ;;
          # Ignore these flags.
        -lang* | -lcrt0.o | -lc | -lgcc | -LANG:=*)
          ;;
        -lkernel32)
          test x"$CYGWIN" != xyes && ac_cv_flibs="$ac_cv_flibs $ac_arg"
          ;;
        -[LRuY])
          # These flags, when seen by themselves, take an argument.
          # We remove the space between option and argument and re-iterate
          # unless we find an empty arg or a new option (starting with -)
	  case $2 in
             "" | -*);;
             *)
		ac_arg="$ac_arg$2"
		shift; shift
		set X $ac_arg "$@"
		;;
	  esac
          ;;
        -YP,*)
          for ac_j in `echo $ac_arg | sed -e 's/-YP,/-L/;s/:/ -L/g'`; do

  ac_exists=false
  for ac_i in $ac_cv_flibs; do
    if test x"$ac_j" = x"$ac_i"; then
      ac_exists=true
      break
    fi
  done

  if test x"$ac_exists" = xtrue; then
  :
else
  ac_arg="$ac_arg $ac_j"
                             ac_cv_flibs="$ac_cv_flibs $ac_j"
fi

          done
          ;;
        -[lLR]*)

  ac_exists=false
  for ac_i in $ac_cv_flibs; do
    if test x"$ac_arg" = x"$ac_i"; then
      ac_exists=true
      break
    fi
  done

  if test x"$ac_exists" = xtrue; then
  :
else
  ac_cv_flibs="$ac_cv_flibs $ac_arg"
fi

          ;;
          # Ignore everything else.
  esac
done
# restore positional arguments
set X $ac_save_positional; shift

# We only consider "LD_RUN_PATH" on Solaris systems.  If this is seen,
# then we insist that the "run path" must be an absolute path (i.e. it
# must begin with a "/").
case `(uname -sr) 2>/dev/null` in
   "SunOS 5"*)
      ac_ld_run_path=`echo $ac_f77_v_output |
                        sed -n 's,^.*LD_RUN_PATH *= *\(/[^ ]*\).*$,-R\1,p'`
      test "x$ac_ld_run_path" != x &&
        if test "$ac_compiler_gnu" = yes; then
  for ac_link_opt in $ac_ld_run_path; do
    ac_cv_flibs="$ac_cv_flibs -Xlinker $ac_link_opt"
  done
else
  ac_cv_flibs="$ac_cv_flibs $ac_ld_run_path"
fi
      ;;
esac
fi # test "x$FLIBS" = "x"

fi
echo "$as_me:2114: result: $ac_cv_flibs" >&5
echo "${ECHO_T}$ac_cv_flibs" >&6
FLIBS="$ac_cv_flibs"

ac_ext=c
ac_cpp='$CPP $CPPFLAGS'
ac_compile='$CC -c $CFLAGS $CPPFLAGS conftest.$ac_ext >&5'
ac_link='$CC -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_c_compiler_gnu

echo "$as_me:2124: checking for dummy main to link with Fortran 77 libraries" >&5
echo $ECHO_N "checking for dummy main to link with Fortran 77 libraries... $ECHO_C" >&6
if test "${ac_cv_f77_dummy_main+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_ext=c
ac_cpp='$CPP $CPPFLAGS'
ac_compile='$CC -c $CFLAGS $CPPFLAGS conftest.$ac_ext >&5'
ac_link='$CC -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_c_compiler_gnu
 ac_f77_dm_save_LIBS=$LIBS
 LIBS="$LIBS $FLIBS"

 # First, try linking without a dummy main:
 cat >conftest.$ac_ext <<_ACEOF
#line 2139 "configure"
#include "confdefs.h"

#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{

  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:2157: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:2160: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:2163: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:2166: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  ac_cv_f77_dummy_main=none
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
ac_cv_f77_dummy_main=unknown
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext

 if test $ac_cv_f77_dummy_main = unknown; then
   for ac_func in MAIN__ MAIN_ __main MAIN _MAIN __MAIN main_ main__ _main; do
     cat >conftest.$ac_ext <<_ACEOF
#line 2179 "configure"
#include "confdefs.h"
#define F77_DUMMY_MAIN $ac_func
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{

  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:2197: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:2200: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:2203: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:2206: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  ac_cv_f77_dummy_main=$ac_func; break
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
   done
 fi
 rm -f conftest*
 LIBS=$ac_f77_dm_save_LIBS
 ac_ext=c
ac_cpp='$CPP $CPPFLAGS'
ac_compile='$CC -c $CFLAGS $CPPFLAGS conftest.$ac_ext >&5'
ac_link='$CC -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_c_compiler_gnu

fi
echo "$as_me:2225: result: $ac_cv_f77_dummy_main" >&5
echo "${ECHO_T}$ac_cv_f77_dummy_main" >&6
F77_DUMMY_MAIN=$ac_cv_f77_dummy_main
if test "$F77_DUMMY_MAIN" != unknown; then
  if test $F77_DUMMY_MAIN != none; then

cat >>confdefs.h <<EOF
#define F77_DUMMY_MAIN $F77_DUMMY_MAIN
EOF

fi
else
  { { echo "$as_me:2237: error: Linking to Fortran libraries from C fails." >&5
echo "$as_me: error: Linking to Fortran libraries from C fails." >&2;}
   { (exit 1); exit 1; }; }
fi

echo "$as_me:2242: checking for Fortran 77 name-mangling scheme" >&5
echo $ECHO_N "checking for Fortran 77 name-mangling scheme... $ECHO_C" >&6
if test "${ac_cv_f77_mangling+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_ext=f
ac_compile='$F77 -c $FFLAGS conftest.$ac_ext >&5'
ac_link='$F77 -o conftest$ac_exeext $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_f77_compiler_gnu
cat >conftest.$ac_ext <<_ACEOF
      subroutine foobar()
      return
      end
      subroutine foo_bar()
      return
      end
_ACEOF
rm -f conftest.$ac_objext
if { (eval echo "$as_me:2260: \"$ac_compile\"") >&5
  (eval $ac_compile) 2>&5
  ac_status=$?
  echo "$as_me:2263: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest.$ac_objext'
  { (eval echo "$as_me:2266: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:2269: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  mv conftest.$ac_objext cf77_test.$ac_objext

  ac_ext=c
ac_cpp='$CPP $CPPFLAGS'
ac_compile='$CC -c $CFLAGS $CPPFLAGS conftest.$ac_ext >&5'
ac_link='$CC -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_c_compiler_gnu

  ac_save_LIBS=$LIBS
  LIBS="cf77_test.$ac_objext $LIBS $FLIBS"

  ac_success=no
  for ac_foobar in foobar FOOBAR; do
    for ac_underscore in "" "_"; do
      ac_func="$ac_foobar$ac_underscore"
      cat >conftest.$ac_ext <<_ACEOF
#line 2287 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $ac_func ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$ac_func ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:2312: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:2315: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:2318: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:2321: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  ac_success=yes; break 2
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
    done
  done

  if test "$ac_success" = "yes"; then
     case $ac_foobar in
        foobar)
           ac_case=lower
           ac_foo_bar=foo_bar
           ;;
        FOOBAR)
           ac_case=upper
           ac_foo_bar=FOO_BAR
           ;;
     esac

     ac_success_extra=no
     for ac_extra in "" "_"; do
        ac_func="$ac_foo_bar$ac_underscore$ac_extra"
        cat >conftest.$ac_ext <<_ACEOF
#line 2348 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $ac_func ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$ac_func ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:2373: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:2376: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:2379: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:2382: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  ac_success_extra=yes; break
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
     done

     if test "$ac_success_extra" = "yes"; then
	ac_cv_f77_mangling="$ac_case case"
        if test -z "$ac_underscore"; then
           ac_cv_f77_mangling="$ac_cv_f77_mangling, no underscore"
	else
           ac_cv_f77_mangling="$ac_cv_f77_mangling, underscore"
        fi
        if test -z "$ac_extra"; then
           ac_cv_f77_mangling="$ac_cv_f77_mangling, no extra underscore"
	else
           ac_cv_f77_mangling="$ac_cv_f77_mangling, extra underscore"
        fi
      else
	ac_cv_f77_mangling="unknown"
      fi
  else
     ac_cv_f77_mangling="unknown"
  fi

  LIBS=$ac_save_LIBS
  ac_ext=f
ac_compile='$F77 -c $FFLAGS conftest.$ac_ext >&5'
ac_link='$F77 -o conftest$ac_exeext $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_f77_compiler_gnu
  rm -f cf77_test* conftest*
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
fi
rm -f conftest.$ac_objext conftest.$ac_ext
ac_ext=c
ac_cpp='$CPP $CPPFLAGS'
ac_compile='$CC -c $CFLAGS $CPPFLAGS conftest.$ac_ext >&5'
ac_link='$CC -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=$ac_cv_c_compiler_gnu

fi
echo "$as_me:2429: result: $ac_cv_f77_mangling" >&5
echo "${ECHO_T}$ac_cv_f77_mangling" >&6

acx_blas_ok=no

# Check whether --with-blas or --without-blas was given.
if test "${with_blas+set}" = set; then
  withval="$with_blas"

fi;
case $with_blas in
	yes | "") ;;
	no) acx_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
case $ac_cv_f77_mangling in
  upper*) ac_val="SGEMM" ;;
  lower*) ac_val="sgemm" ;;
  *)      ac_val="unknown" ;;
esac
case $ac_cv_f77_mangling in *," underscore"*) ac_val="$ac_val"_ ;; esac

sgemm="$ac_val"

case $ac_cv_f77_mangling in
  upper*) ac_val="DGEMM" ;;
  lower*) ac_val="dgemm" ;;
  *)      ac_val="unknown" ;;
esac
case $ac_cv_f77_mangling in *," underscore"*) ac_val="$ac_val"_ ;; esac

dgemm="$ac_val"

acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $acx_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	echo "$as_me:2472: checking for $sgemm in $BLAS_LIBS" >&5
echo $ECHO_N "checking for $sgemm in $BLAS_LIBS... $ECHO_C" >&6
	cat >conftest.$ac_ext <<_ACEOF
#line 2475 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $sgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$sgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:2500: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:2503: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:2506: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:2509: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  acx_blas_ok=yes
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
BLAS_LIBS=""
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
	echo "$as_me:2518: result: $acx_blas_ok" >&5
echo "${ECHO_T}$acx_blas_ok" >&6
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	as_ac_var=`echo "ac_cv_func_$sgemm" | $as_tr_sh`
echo "$as_me:2528: checking for $sgemm" >&5
echo $ECHO_N "checking for $sgemm... $ECHO_C" >&6
if eval "test \"\${$as_ac_var+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  cat >conftest.$ac_ext <<_ACEOF
#line 2534 "configure"
#include "confdefs.h"
/* System header to define __stub macros and hopefully few prototypes,
    which can conflict with char $sgemm (); below.  */
#include <assert.h>
/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $sgemm ();
char (*f) ();

#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
/* The GNU C library defines this for functions which it implements
    to always fail with ENOSYS.  Some functions are actually named
    something starting with __ and the normal name is an alias.  */
#if defined (__stub_$sgemm) || defined (__stub___$sgemm)
choke me
#else
f = $sgemm;
#endif

  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:2571: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:2574: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:2577: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:2580: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_var=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_var=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
fi
echo "$as_me:2590: result: `eval echo '${'$as_ac_var'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_var'}'`" >&6
if test `eval echo '${'$as_ac_var'}'` = yes; then
  acx_blas_ok=yes
fi

	LIBS="$save_LIBS"
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_blas_ok = no; then
	echo "$as_me:2601: checking for ATL_xerbla in -latlas" >&5
echo $ECHO_N "checking for ATL_xerbla in -latlas... $ECHO_C" >&6
if test "${ac_cv_lib_atlas_ATL_xerbla+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-latlas  $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 2609 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char ATL_xerbla ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
ATL_xerbla ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:2634: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:2637: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:2640: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:2643: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  ac_cv_lib_atlas_ATL_xerbla=yes
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
ac_cv_lib_atlas_ATL_xerbla=no
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:2654: result: $ac_cv_lib_atlas_ATL_xerbla" >&5
echo "${ECHO_T}$ac_cv_lib_atlas_ATL_xerbla" >&6
if test $ac_cv_lib_atlas_ATL_xerbla = yes; then
  as_ac_Lib=`echo "ac_cv_lib_f77blas_$sgemm" | $as_tr_sh`
echo "$as_me:2658: checking for $sgemm in -lf77blas" >&5
echo $ECHO_N "checking for $sgemm in -lf77blas... $ECHO_C" >&6
if eval "test \"\${$as_ac_Lib+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lf77blas -latlas $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 2666 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $sgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$sgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:2691: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:2694: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:2697: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:2700: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_Lib=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_Lib=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:2711: result: `eval echo '${'$as_ac_Lib'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_Lib'}'`" >&6
if test `eval echo '${'$as_ac_Lib'}'` = yes; then
  echo "$as_me:2714: checking for cblas_dgemm in -lcblas" >&5
echo $ECHO_N "checking for cblas_dgemm in -lcblas... $ECHO_C" >&6
if test "${ac_cv_lib_cblas_cblas_dgemm+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lcblas -lf77blas -latlas $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 2722 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char cblas_dgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
cblas_dgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:2747: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:2750: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:2753: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:2756: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  ac_cv_lib_cblas_cblas_dgemm=yes
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
ac_cv_lib_cblas_cblas_dgemm=no
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:2767: result: $ac_cv_lib_cblas_cblas_dgemm" >&5
echo "${ECHO_T}$ac_cv_lib_cblas_cblas_dgemm" >&6
if test $ac_cv_lib_cblas_cblas_dgemm = yes; then
  acx_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"
fi

fi

fi

fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	as_ac_Lib=`echo "ac_cv_lib_blas_$sgemm" | $as_tr_sh`
echo "$as_me:2783: checking for $sgemm in -lblas" >&5
echo $ECHO_N "checking for $sgemm in -lblas... $ECHO_C" >&6
if eval "test \"\${$as_ac_Lib+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lblas  $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 2791 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $sgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$sgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:2816: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:2819: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:2822: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:2825: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_Lib=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_Lib=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:2836: result: `eval echo '${'$as_ac_Lib'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_Lib'}'`" >&6
if test `eval echo '${'$as_ac_Lib'}'` = yes; then
  as_ac_Lib=`echo "ac_cv_lib_dgemm_$dgemm" | $as_tr_sh`
echo "$as_me:2840: checking for $dgemm in -ldgemm" >&5
echo $ECHO_N "checking for $dgemm in -ldgemm... $ECHO_C" >&6
if eval "test \"\${$as_ac_Lib+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-ldgemm -lblas $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 2848 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $dgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$dgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:2873: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:2876: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:2879: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:2882: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_Lib=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_Lib=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:2893: result: `eval echo '${'$as_ac_Lib'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_Lib'}'`" >&6
if test `eval echo '${'$as_ac_Lib'}'` = yes; then
  as_ac_Lib=`echo "ac_cv_lib_sgemm_$sgemm" | $as_tr_sh`
echo "$as_me:2897: checking for $sgemm in -lsgemm" >&5
echo $ECHO_N "checking for $sgemm in -lsgemm... $ECHO_C" >&6
if eval "test \"\${$as_ac_Lib+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lsgemm -lblas $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 2905 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $sgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$sgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:2930: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:2933: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:2936: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:2939: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_Lib=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_Lib=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:2950: result: `eval echo '${'$as_ac_Lib'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_Lib'}'`" >&6
if test `eval echo '${'$as_ac_Lib'}'` = yes; then
  acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"
fi

fi

fi

fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
	as_ac_Lib=`echo "ac_cv_lib_cxml_$sgemm" | $as_tr_sh`
echo "$as_me:2965: checking for $sgemm in -lcxml" >&5
echo $ECHO_N "checking for $sgemm in -lcxml... $ECHO_C" >&6
if eval "test \"\${$as_ac_Lib+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lcxml  $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 2973 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $sgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$sgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:2998: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:3001: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:3004: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:3007: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_Lib=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_Lib=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:3018: result: `eval echo '${'$as_ac_Lib'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_Lib'}'`" >&6
if test `eval echo '${'$as_ac_Lib'}'` = yes; then
  acx_blas_ok=yes;BLAS_LIBS="-lcxml"
fi

fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
	as_ac_Lib=`echo "ac_cv_lib_dxml_$sgemm" | $as_tr_sh`
echo "$as_me:3029: checking for $sgemm in -ldxml" >&5
echo $ECHO_N "checking for $sgemm in -ldxml... $ECHO_C" >&6
if eval "test \"\${$as_ac_Lib+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-ldxml  $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 3037 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $sgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$sgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:3062: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:3065: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:3068: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:3071: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_Lib=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_Lib=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:3082: result: `eval echo '${'$as_ac_Lib'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_Lib'}'`" >&6
if test `eval echo '${'$as_ac_Lib'}'` = yes; then
  acx_blas_ok=yes;BLAS_LIBS="-ldxml"
fi

fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		echo "$as_me:3093: checking for acosp in -lsunmath" >&5
echo $ECHO_N "checking for acosp in -lsunmath... $ECHO_C" >&6
if test "${ac_cv_lib_sunmath_acosp+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lsunmath  $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 3101 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char acosp ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
acosp ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:3126: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:3129: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:3132: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:3135: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  ac_cv_lib_sunmath_acosp=yes
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
ac_cv_lib_sunmath_acosp=no
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:3146: result: $ac_cv_lib_sunmath_acosp" >&5
echo "${ECHO_T}$ac_cv_lib_sunmath_acosp" >&6
if test $ac_cv_lib_sunmath_acosp = yes; then
  as_ac_Lib=`echo "ac_cv_lib_sunperf_$sgemm" | $as_tr_sh`
echo "$as_me:3150: checking for $sgemm in -lsunperf" >&5
echo $ECHO_N "checking for $sgemm in -lsunperf... $ECHO_C" >&6
if eval "test \"\${$as_ac_Lib+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lsunperf -lsunmath $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 3158 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $sgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$sgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:3183: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:3186: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:3189: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:3192: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_Lib=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_Lib=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:3203: result: `eval echo '${'$as_ac_Lib'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_Lib'}'`" >&6
if test `eval echo '${'$as_ac_Lib'}'` = yes; then
  BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 acx_blas_ok=yes
fi

fi

	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
	as_ac_Lib=`echo "ac_cv_lib_scs_$sgemm" | $as_tr_sh`
echo "$as_me:3218: checking for $sgemm in -lscs" >&5
echo $ECHO_N "checking for $sgemm in -lscs... $ECHO_C" >&6
if eval "test \"\${$as_ac_Lib+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lscs  $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 3226 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $sgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$sgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:3251: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:3254: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:3257: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:3260: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_Lib=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_Lib=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:3271: result: `eval echo '${'$as_ac_Lib'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_Lib'}'`" >&6
if test `eval echo '${'$as_ac_Lib'}'` = yes; then
  acx_blas_ok=yes; BLAS_LIBS="-lscs"
fi

fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
	as_ac_Lib=`echo "ac_cv_lib_complib.sgimath_$sgemm" | $as_tr_sh`
echo "$as_me:3282: checking for $sgemm in -lcomplib.sgimath" >&5
echo $ECHO_N "checking for $sgemm in -lcomplib.sgimath... $ECHO_C" >&6
if eval "test \"\${$as_ac_Lib+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lcomplib.sgimath  $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 3290 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $sgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$sgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:3315: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:3318: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:3321: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:3324: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_Lib=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_Lib=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:3335: result: `eval echo '${'$as_ac_Lib'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_Lib'}'`" >&6
if test `eval echo '${'$as_ac_Lib'}'` = yes; then
  acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"
fi

fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	as_ac_Lib=`echo "ac_cv_lib_blas_$sgemm" | $as_tr_sh`
echo "$as_me:3346: checking for $sgemm in -lblas" >&5
echo $ECHO_N "checking for $sgemm in -lblas... $ECHO_C" >&6
if eval "test \"\${$as_ac_Lib+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lblas  $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 3354 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $sgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$sgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:3379: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:3382: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:3385: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:3388: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_Lib=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_Lib=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:3399: result: `eval echo '${'$as_ac_Lib'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_Lib'}'`" >&6
if test `eval echo '${'$as_ac_Lib'}'` = yes; then
  as_ac_Lib=`echo "ac_cv_lib_essl_$sgemm" | $as_tr_sh`
echo "$as_me:3403: checking for $sgemm in -lessl" >&5
echo $ECHO_N "checking for $sgemm in -lessl... $ECHO_C" >&6
if eval "test \"\${$as_ac_Lib+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lessl -lblas $FLIBS $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 3411 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $sgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$sgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:3436: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:3439: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:3442: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:3445: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_Lib=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_Lib=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:3456: result: `eval echo '${'$as_ac_Lib'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_Lib'}'`" >&6
if test `eval echo '${'$as_ac_Lib'}'` = yes; then
  acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"
fi

fi

fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
	as_ac_Lib=`echo "ac_cv_lib_blas_$sgemm" | $as_tr_sh`
echo "$as_me:3469: checking for $sgemm in -lblas" >&5
echo $ECHO_N "checking for $sgemm in -lblas... $ECHO_C" >&6
if eval "test \"\${$as_ac_Lib+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lblas  $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 3477 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $sgemm ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$sgemm ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:3502: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:3505: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:3508: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:3511: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_Lib=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_Lib=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:3522: result: `eval echo '${'$as_ac_Lib'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_Lib'}'`" >&6
if test `eval echo '${'$as_ac_Lib'}'` = yes; then
  acx_blas_ok=yes; BLAS_LIBS="-lblas"
fi

fi

LIBS="$acx_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blas_ok" = xyes; then

cat >>confdefs.h <<\EOF
#define HAVE_BLAS 1
EOF

        :
else
        acx_blas_ok=no
        { { echo "$as_me:3542: error: could not find your BLAS library" >&5
echo "$as_me: error: could not find your BLAS library" >&2;}
   { (exit 1); exit 1; }; }
fi

acx_lapack_ok=no

# Check whether --with-lapack or --without-lapack was given.
if test "${with_lapack+set}" = set; then
  withval="$with_lapack"

fi;
case $with_lapack in
	yes | "") ;;
	no) acx_lapack_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
	*) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
case $ac_cv_f77_mangling in
  upper*) ac_val="CHEEV" ;;
  lower*) ac_val="cheev" ;;
  *)      ac_val="unknown" ;;
esac
case $ac_cv_f77_mangling in *," underscore"*) ac_val="$ac_val"_ ;; esac

cheev="$ac_val"

# We cannot use LAPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
	acx_lapack_ok=noblas
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
	echo "$as_me:3579: checking for $cheev in $LAPACK_LIBS" >&5
echo $ECHO_N "checking for $cheev in $LAPACK_LIBS... $ECHO_C" >&6
	cat >conftest.$ac_ext <<_ACEOF
#line 3582 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $cheev ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$cheev ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:3607: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:3610: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:3613: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:3616: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  acx_lapack_ok=yes
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
LAPACK_LIBS=""
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
	echo "$as_me:3625: result: $acx_lapack_ok" >&5
echo "${ECHO_T}$acx_lapack_ok" >&6
	LIBS="$save_LIBS"
	if test acx_lapack_ok = no; then
		LAPACK_LIBS=""
	fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
	as_ac_var=`echo "ac_cv_func_$cheev" | $as_tr_sh`
echo "$as_me:3637: checking for $cheev" >&5
echo $ECHO_N "checking for $cheev... $ECHO_C" >&6
if eval "test \"\${$as_ac_var+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  cat >conftest.$ac_ext <<_ACEOF
#line 3643 "configure"
#include "confdefs.h"
/* System header to define __stub macros and hopefully few prototypes,
    which can conflict with char $cheev (); below.  */
#include <assert.h>
/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $cheev ();
char (*f) ();

#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
/* The GNU C library defines this for functions which it implements
    to always fail with ENOSYS.  Some functions are actually named
    something starting with __ and the normal name is an alias.  */
#if defined (__stub_$cheev) || defined (__stub___$cheev)
choke me
#else
f = $cheev;
#endif

  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:3680: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:3683: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:3686: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:3689: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_var=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_var=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
fi
echo "$as_me:3699: result: `eval echo '${'$as_ac_var'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_var'}'`" >&6
if test `eval echo '${'$as_ac_var'}'` = yes; then
  acx_lapack_ok=yes
fi

	LIBS="$save_LIBS"
fi

# Generic LAPACK library?
if test $acx_lapack_ok = no; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	as_ac_Lib=`echo "ac_cv_lib_lapack_$cheev" | $as_tr_sh`
echo "$as_me:3712: checking for $cheev in -llapack" >&5
echo $ECHO_N "checking for $cheev in -llapack... $ECHO_C" >&6
if eval "test \"\${$as_ac_Lib+set}\" = set"; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-llapack $FLIBS $LIBS"
cat >conftest.$ac_ext <<_ACEOF
#line 3720 "configure"
#include "confdefs.h"

/* Override any gcc2 internal prototype to avoid an error.  */
#ifdef __cplusplus
extern "C"
#endif
/* We use char because int might match the return type of a gcc2
   builtin and then its argument prototype would still apply.  */
char $cheev ();
#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif
int
main ()
{
$cheev ();
  ;
  return 0;
}
_ACEOF
rm -f conftest.$ac_objext conftest$ac_exeext
if { (eval echo "$as_me:3745: \"$ac_link\"") >&5
  (eval $ac_link) 2>&5
  ac_status=$?
  echo "$as_me:3748: \$? = $ac_status" >&5
  (exit $ac_status); } &&
         { ac_try='test -s conftest$ac_exeext'
  { (eval echo "$as_me:3751: \"$ac_try\"") >&5
  (eval $ac_try) 2>&5
  ac_status=$?
  echo "$as_me:3754: \$? = $ac_status" >&5
  (exit $ac_status); }; }; then
  eval "$as_ac_Lib=yes"
else
  echo "$as_me: failed program was:" >&5
cat conftest.$ac_ext >&5
eval "$as_ac_Lib=no"
fi
rm -f conftest.$ac_objext conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
echo "$as_me:3765: result: `eval echo '${'$as_ac_Lib'}'`" >&5
echo "${ECHO_T}`eval echo '${'$as_ac_Lib'}'`" >&6
if test `eval echo '${'$as_ac_Lib'}'` = yes; then
  acx_lapack_ok=yes; LAPACK_LIBS="-llapack"
fi

	LIBS="$save_LIBS"
fi

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then

cat >>confdefs.h <<\EOF
#define HAVE_LAPACK 1
EOF

        :
else
        acx_lapack_ok=no
        { { echo "$as_me:3784: error: could not find your LAPACK library" >&5
echo "$as_me: error: could not find your LAPACK library" >&2;}
   { (exit 1); exit 1; }; }
fi

ac_config_files="$ac_config_files src/Makevars"
cat >confcache <<\_ACEOF
# This file is a shell script that caches the results of configure
# tests run on this system so they can be shared between configure
# scripts and configure runs, see configure's option --config-cache.
# It is not useful on other systems.  If it contains results you don't
# want to keep, you may remove or edit it.
#
# config.status only pays attention to the cache file if you give it
# the --recheck option to rerun configure.
#
# `ac_cv_env_foo' variables (set or unset) will be overriden when
# loading this file, other *unset* `ac_cv_foo' will be assigned the
# following values.

_ACEOF

# The following way of writing the cache mishandles newlines in values,
# but we know of no workaround that is simple, portable, and efficient.
# So, don't put newlines in cache variables' values.
# Ultrix sh set writes to stderr and can't be redirected directly,
# and sets the high bit in the cache file unless we assign to the vars.
{
  (set) 2>&1 |
    case `(ac_space=' '; set | grep ac_space) 2>&1` in
    *ac_space=\ *)
      # `set' does not quote correctly, so add quotes (double-quote
      # substitution turns \\\\ into \\, and sed turns \\ into \).
      sed -n \
        "s/'/'\\\\''/g;
    	  s/^\\([_$as_cr_alnum]*_cv_[_$as_cr_alnum]*\\)=\\(.*\\)/\\1='\\2'/p"
      ;;
    *)
      # `set' quotes correctly as required by POSIX, so do not add quotes.
      sed -n \
        "s/^\\([_$as_cr_alnum]*_cv_[_$as_cr_alnum]*\\)=\\(.*\\)/\\1=\\2/p"
      ;;
    esac;
} |
  sed '
     t clear
     : clear
     s/^\([^=]*\)=\(.*[{}].*\)$/test "${\1+set}" = set || &/
     t end
     /^ac_cv_env/!s/^\([^=]*\)=\(.*\)$/\1=${\1=\2}/
     : end' >>confcache
if cmp -s $cache_file confcache; then :; else
  if test -w $cache_file; then
    test "x$cache_file" != "x/dev/null" && echo "updating cache $cache_file"
    cat confcache >$cache_file
  else
    echo "not updating unwritable cache $cache_file"
  fi
fi
rm -f confcache

test "x$prefix" = xNONE && prefix=$ac_default_prefix
# Let make expand exec_prefix.
test "x$exec_prefix" = xNONE && exec_prefix='${prefix}'

# VPATH may cause trouble with some makes, so we remove $(srcdir),
# ${srcdir} and @srcdir@ from VPATH if srcdir is ".", strip leading and
# trailing colons and then remove the whole line if VPATH becomes empty
# (actually we leave an empty line to preserve line numbers).
if test "x$srcdir" = x.; then
  ac_vpsub='/^[ 	]*VPATH[ 	]*=/{
s/:*\$(srcdir):*/:/;
s/:*\${srcdir}:*/:/;
s/:*@srcdir@:*/:/;
s/^\([^=]*=[ 	]*\):*/\1/;
s/:*$//;
s/^[^=]*=[ 	]*$//;
}'
fi

# Transform confdefs.h into DEFS.
# Protect against shell expansion while executing Makefile rules.
# Protect against Makefile macro expansion.
#
# If the first sed substitution is executed (which looks for macros that
# take arguments), then we branch to the quote section.  Otherwise,
# look for a macro that doesn't take arguments.
cat >confdef2opt.sed <<\EOF
t clear
: clear
s,^[ 	]*#[ 	]*define[ 	][ 	]*\([^ 	(][^ 	(]*([^)]*)\)[ 	]*\(.*\),-D\1=\2,g
t quote
s,^[ 	]*#[ 	]*define[ 	][ 	]*\([^ 	][^ 	]*\)[ 	]*\(.*\),-D\1=\2,g
t quote
d
: quote
s,[ 	`~#$^&*(){}\\|;'"<>?],\\&,g
s,\[,\\&,g
s,\],\\&,g
s,\$,$$,g
p
EOF
# We use echo to avoid assuming a particular line-breaking character.
# The extra dot is to prevent the shell from consuming trailing
# line-breaks from the sub-command output.  A line-break within
# single-quotes doesn't work because, if this script is created in a
# platform that uses two characters for line-breaks (e.g., DOS), tr
# would break.
ac_LF_and_DOT=`echo; echo .`
DEFS=`sed -n -f confdef2opt.sed confdefs.h | tr "$ac_LF_and_DOT" ' .'`
rm -f confdef2opt.sed

: ${CONFIG_STATUS=./config.status}
ac_clean_files_save=$ac_clean_files
ac_clean_files="$ac_clean_files $CONFIG_STATUS"
{ echo "$as_me:3899: creating $CONFIG_STATUS" >&5
echo "$as_me: creating $CONFIG_STATUS" >&6;}
cat >$CONFIG_STATUS <<_ACEOF
#! $SHELL
# Generated automatically by configure.
# Run this file to recreate the current configuration.
# Compiler output produced by configure, useful for debugging
# configure, is in config.log if it exists.

debug=false
SHELL=\${CONFIG_SHELL-$SHELL}
ac_cs_invocation="\$0 \$@"

_ACEOF

cat >>$CONFIG_STATUS <<\_ACEOF
# Be Bourne compatible
if test -n "${ZSH_VERSION+set}" && (emulate sh) >/dev/null 2>&1; then
  emulate sh
  NULLCMD=:
elif test -n "${BASH_VERSION+set}" && (set -o posix) >/dev/null 2>&1; then
  set -o posix
fi

# Name of the executable.
as_me=`echo "$0" |sed 's,.*[\\/],,'`

if expr a : '\(a\)' >/dev/null 2>&1; then
  as_expr=expr
else
  as_expr=false
fi

rm -f conf$$ conf$$.exe conf$$.file
echo >conf$$.file
if ln -s conf$$.file conf$$ 2>/dev/null; then
  # We could just check for DJGPP; but this test a) works b) is more generic
  # and c) will remain valid once DJGPP supports symlinks (DJGPP 2.04).
  if test -f conf$$.exe; then
    # Don't use ln at all; we don't have any links
    as_ln_s='cp -p'
  else
    as_ln_s='ln -s'
  fi
elif ln conf$$.file conf$$ 2>/dev/null; then
  as_ln_s=ln
else
  as_ln_s='cp -p'
fi
rm -f conf$$ conf$$.exe conf$$.file

as_executable_p="test -f"

# Support unset when possible.
if (FOO=FOO; unset FOO) >/dev/null 2>&1; then
  as_unset=unset
else
  as_unset=false
fi

# NLS nuisances.
$as_unset LANG || test "${LANG+set}" != set || { LANG=C; export LANG; }
$as_unset LC_ALL || test "${LC_ALL+set}" != set || { LC_ALL=C; export LC_ALL; }
$as_unset LC_TIME || test "${LC_TIME+set}" != set || { LC_TIME=C; export LC_TIME; }
$as_unset LC_CTYPE || test "${LC_CTYPE+set}" != set || { LC_CTYPE=C; export LC_CTYPE; }
$as_unset LANGUAGE || test "${LANGUAGE+set}" != set || { LANGUAGE=C; export LANGUAGE; }
$as_unset LC_COLLATE || test "${LC_COLLATE+set}" != set || { LC_COLLATE=C; export LC_COLLATE; }
$as_unset LC_NUMERIC || test "${LC_NUMERIC+set}" != set || { LC_NUMERIC=C; export LC_NUMERIC; }
$as_unset LC_MESSAGES || test "${LC_MESSAGES+set}" != set || { LC_MESSAGES=C; export LC_MESSAGES; }

# IFS
# We need space, tab and new line, in precisely that order.
as_nl='
'
IFS=" 	$as_nl"

# CDPATH.
$as_unset CDPATH || test "${CDPATH+set}" != set || { CDPATH=:; export CDPATH; }

exec 6>&1

_ACEOF

# Files that config.status was made for.
if test -n "$ac_config_files"; then
  echo "config_files=\"$ac_config_files\"" >>$CONFIG_STATUS
fi

if test -n "$ac_config_headers"; then
  echo "config_headers=\"$ac_config_headers\"" >>$CONFIG_STATUS
fi

if test -n "$ac_config_links"; then
  echo "config_links=\"$ac_config_links\"" >>$CONFIG_STATUS
fi

if test -n "$ac_config_commands"; then
  echo "config_commands=\"$ac_config_commands\"" >>$CONFIG_STATUS
fi

cat >>$CONFIG_STATUS <<\EOF

ac_cs_usage="\
\`$as_me' instantiates files from templates according to the
current configuration.

Usage: $0 [OPTIONS] [FILE]...

  -h, --help       print this help, then exit
  -V, --version    print version number, then exit
  -d, --debug      don't remove temporary files
      --recheck    update $as_me by reconfiguring in the same conditions
  --file=FILE[:TEMPLATE]
                   instantiate the configuration file FILE

Configuration files:
$config_files

Report bugs to <bug-autoconf@gnu.org>."
EOF

cat >>$CONFIG_STATUS <<EOF
ac_cs_version="\\
config.status
configured by $0, generated by GNU Autoconf 2.52,
  with options \\"`echo "$ac_configure_args" | sed 's/[\\""\`\$]/\\\\&/g'`\\"

Copyright 1992, 1993, 1994, 1995, 1996, 1998, 1999, 2000, 2001
Free Software Foundation, Inc.
This config.status script is free software; the Free Software Foundation
gives unlimited permission to copy, distribute and modify it."
srcdir=$srcdir
EOF

cat >>$CONFIG_STATUS <<\EOF
# If no file are specified by the user, then we need to provide default
# value.  By we need to know if files were specified by the user.
ac_need_defaults=:
while test $# != 0
do
  case $1 in
  --*=*)
    ac_option=`expr "x$1" : 'x\([^=]*\)='`
    ac_optarg=`expr "x$1" : 'x[^=]*=\(.*\)'`
    shift
    set dummy "$ac_option" "$ac_optarg" ${1+"$@"}
    shift
    ;;
  -*);;
  *) # This is not an option, so the user has probably given explicit
     # arguments.
     ac_need_defaults=false;;
  esac

  case $1 in
  # Handling of the options.
EOF
cat >>$CONFIG_STATUS <<EOF
  -recheck | --recheck | --rechec | --reche | --rech | --rec | --re | --r)
    echo "running $SHELL $0 " $ac_configure_args " --no-create --no-recursion"
    exec $SHELL $0 $ac_configure_args --no-create --no-recursion ;;
EOF
cat >>$CONFIG_STATUS <<\EOF
  --version | --vers* | -V )
    echo "$ac_cs_version"; exit 0 ;;
  --he | --h)
    # Conflict between --help and --header
    { { echo "$as_me:4066: error: ambiguous option: $1
Try \`$0 --help' for more information." >&5
echo "$as_me: error: ambiguous option: $1
Try \`$0 --help' for more information." >&2;}
   { (exit 1); exit 1; }; };;
  --help | --hel | -h )
    echo "$ac_cs_usage"; exit 0 ;;
  --debug | --d* | -d )
    debug=: ;;
  --file | --fil | --fi | --f )
    shift
    CONFIG_FILES="$CONFIG_FILES $1"
    ac_need_defaults=false;;
  --header | --heade | --head | --hea )
    shift
    CONFIG_HEADERS="$CONFIG_HEADERS $1"
    ac_need_defaults=false;;

  # This is an error.
  -*) { { echo "$as_me:4085: error: unrecognized option: $1
Try \`$0 --help' for more information." >&5
echo "$as_me: error: unrecognized option: $1
Try \`$0 --help' for more information." >&2;}
   { (exit 1); exit 1; }; } ;;

  *) ac_config_targets="$ac_config_targets $1" ;;

  esac
  shift
done

exec 5>>config.log
cat >&5 << _ACEOF

## ----------------------- ##
## Running config.status.  ##
## ----------------------- ##

This file was extended by $as_me 2.52, executed with
  CONFIG_FILES    = $CONFIG_FILES
  CONFIG_HEADERS  = $CONFIG_HEADERS
  CONFIG_LINKS    = $CONFIG_LINKS
  CONFIG_COMMANDS = $CONFIG_COMMANDS
  > $ac_cs_invocation
on `(hostname || uname -n) 2>/dev/null | sed 1q`

_ACEOF
EOF

cat >>$CONFIG_STATUS <<\EOF
for ac_config_target in $ac_config_targets
do
  case "$ac_config_target" in
  # Handling of arguments.
  "src/Makevars" ) CONFIG_FILES="$CONFIG_FILES src/Makevars" ;;
  *) { { echo "$as_me:4121: error: invalid argument: $ac_config_target" >&5
echo "$as_me: error: invalid argument: $ac_config_target" >&2;}
   { (exit 1); exit 1; }; };;
  esac
done

# If the user did not use the arguments to specify the items to instantiate,
# then the envvar interface is used.  Set only those that are not.
# We use the long form for the default assignment because of an extremely
# bizarre bug on SunOS 4.1.3.
if $ac_need_defaults; then
  test "${CONFIG_FILES+set}" = set || CONFIG_FILES=$config_files
fi

# Create a temporary directory, and hook for its removal unless debugging.
$debug ||
{
  trap 'exit_status=$?; rm -rf $tmp && exit $exit_status' 0
  trap '{ (exit 1); exit 1; }' 1 2 13 15
}

# Create a (secure) tmp directory for tmp files.
: ${TMPDIR=/tmp}
{
  tmp=`(umask 077 && mktemp -d -q "$TMPDIR/csXXXXXX") 2>/dev/null` &&
  test -n "$tmp" && test -d "$tmp"
}  ||
{
  tmp=$TMPDIR/cs$$-$RANDOM
  (umask 077 && mkdir $tmp)
} ||
{
   echo "$me: cannot create a temporary directory in $TMPDIR" >&2
   { (exit 1); exit 1; }
}

EOF

cat >>$CONFIG_STATUS <<EOF

#
# CONFIG_FILES section.
#

# No need to generate the scripts if there are no CONFIG_FILES.
# This happens for instance when ./config.status config.h
if test -n "\$CONFIG_FILES"; then
  # Protect against being on the right side of a sed subst in config.status.
  sed 's/,@/@@/; s/@,/@@/; s/,;t t\$/@;t t/; /@;t t\$/s/[\\\\&,]/\\\\&/g;
   s/@@/,@/; s/@@/@,/; s/@;t t\$/,;t t/' >\$tmp/subs.sed <<\\CEOF
s,@SHELL@,$SHELL,;t t
s,@exec_prefix@,$exec_prefix,;t t
s,@prefix@,$prefix,;t t
s,@program_transform_name@,$program_transform_name,;t t
s,@bindir@,$bindir,;t t
s,@sbindir@,$sbindir,;t t
s,@libexecdir@,$libexecdir,;t t
s,@datadir@,$datadir,;t t
s,@sysconfdir@,$sysconfdir,;t t
s,@sharedstatedir@,$sharedstatedir,;t t
s,@localstatedir@,$localstatedir,;t t
s,@libdir@,$libdir,;t t
s,@includedir@,$includedir,;t t
s,@oldincludedir@,$oldincludedir,;t t
s,@infodir@,$infodir,;t t
s,@mandir@,$mandir,;t t
s,@PACKAGE_NAME@,$PACKAGE_NAME,;t t
s,@PACKAGE_TARNAME@,$PACKAGE_TARNAME,;t t
s,@PACKAGE_VERSION@,$PACKAGE_VERSION,;t t
s,@PACKAGE_STRING@,$PACKAGE_STRING,;t t
s,@PACKAGE_BUGREPORT@,$PACKAGE_BUGREPORT,;t t
s,@build_alias@,$build_alias,;t t
s,@host_alias@,$host_alias,;t t
s,@target_alias@,$target_alias,;t t
s,@ECHO_C@,$ECHO_C,;t t
s,@ECHO_N@,$ECHO_N,;t t
s,@ECHO_T@,$ECHO_T,;t t
s,@PATH_SEPARATOR@,$PATH_SEPARATOR,;t t
s,@DEFS@,$DEFS,;t t
s,@LIBS@,$LIBS,;t t
s,@CC@,$CC,;t t
s,@CFLAGS@,$CFLAGS,;t t
s,@LDFLAGS@,$LDFLAGS,;t t
s,@CPPFLAGS@,$CPPFLAGS,;t t
s,@ac_ct_CC@,$ac_ct_CC,;t t
s,@EXEEXT@,$EXEEXT,;t t
s,@OBJEXT@,$OBJEXT,;t t
s,@F77@,$F77,;t t
s,@FFLAGS@,$FFLAGS,;t t
s,@ac_ct_F77@,$ac_ct_F77,;t t
s,@FLIBS@,$FLIBS,;t t
s,@BLAS_LIBS@,$BLAS_LIBS,;t t
s,@LAPACK_LIBS@,$LAPACK_LIBS,;t t
CEOF

EOF

  cat >>$CONFIG_STATUS <<\EOF
  # Split the substitutions into bite-sized pieces for seds with
  # small command number limits, like on Digital OSF/1 and HP-UX.
  ac_max_sed_lines=48
  ac_sed_frag=1 # Number of current file.
  ac_beg=1 # First line for current file.
  ac_end=$ac_max_sed_lines # Line after last line for current file.
  ac_more_lines=:
  ac_sed_cmds=
  while $ac_more_lines; do
    if test $ac_beg -gt 1; then
      sed "1,${ac_beg}d; ${ac_end}q" $tmp/subs.sed >$tmp/subs.frag
    else
      sed "${ac_end}q" $tmp/subs.sed >$tmp/subs.frag
    fi
    if test ! -s $tmp/subs.frag; then
      ac_more_lines=false
    else
      # The purpose of the label and of the branching condition is to
      # speed up the sed processing (if there are no `@' at all, there
      # is no need to browse any of the substitutions).
      # These are the two extra sed commands mentioned above.
      (echo ':t
  /@[a-zA-Z_][a-zA-Z_0-9]*@/!b' && cat $tmp/subs.frag) >$tmp/subs-$ac_sed_frag.sed
      if test -z "$ac_sed_cmds"; then
  	ac_sed_cmds="sed -f $tmp/subs-$ac_sed_frag.sed"
      else
  	ac_sed_cmds="$ac_sed_cmds | sed -f $tmp/subs-$ac_sed_frag.sed"
      fi
      ac_sed_frag=`expr $ac_sed_frag + 1`
      ac_beg=$ac_end
      ac_end=`expr $ac_end + $ac_max_sed_lines`
    fi
  done
  if test -z "$ac_sed_cmds"; then
    ac_sed_cmds=cat
  fi
fi # test -n "$CONFIG_FILES"

EOF
cat >>$CONFIG_STATUS <<\EOF
for ac_file in : $CONFIG_FILES; do test "x$ac_file" = x: && continue
  # Support "outfile[:infile[:infile...]]", defaulting infile="outfile.in".
  case $ac_file in
  - | *:- | *:-:* ) # input from stdin
        cat >$tmp/stdin
        ac_file_in=`echo "$ac_file" | sed 's,[^:]*:,,'`
        ac_file=`echo "$ac_file" | sed 's,:.*,,'` ;;
  *:* ) ac_file_in=`echo "$ac_file" | sed 's,[^:]*:,,'`
        ac_file=`echo "$ac_file" | sed 's,:.*,,'` ;;
  * )   ac_file_in=$ac_file.in ;;
  esac

  # Compute @srcdir@, @top_srcdir@, and @INSTALL@ for subdirectories.
  ac_dir=`$as_expr X"$ac_file" : 'X\(.*[^/]\)//*[^/][^/]*/*$' \| \
         X"$ac_file" : 'X\(//\)[^/]' \| \
         X"$ac_file" : 'X\(//\)$' \| \
         X"$ac_file" : 'X\(/\)' \| \
         .     : '\(.\)' 2>/dev/null ||
echo X"$ac_file" |
    sed '/^X\(.*[^/]\)\/\/*[^/][^/]*\/*$/{ s//\1/; q; }
  	  /^X\(\/\/\)[^/].*/{ s//\1/; q; }
  	  /^X\(\/\/\)$/{ s//\1/; q; }
  	  /^X\(\/\).*/{ s//\1/; q; }
  	  s/.*/./; q'`
  if test "$ac_dir" != "$ac_file" && test "$ac_dir" != .; then
    { case "$ac_dir" in
  [\\/]* | ?:[\\/]* ) as_incr_dir=;;
  *)                      as_incr_dir=.;;
esac
as_dummy="$ac_dir"
for as_mkdir_dir in `IFS='/\\'; set X $as_dummy; shift; echo "$@"`; do
  case $as_mkdir_dir in
    # Skip DOS drivespec
    ?:) as_incr_dir=$as_mkdir_dir ;;
    *)
      as_incr_dir=$as_incr_dir/$as_mkdir_dir
      test -d "$as_incr_dir" || mkdir "$as_incr_dir"
    ;;
  esac
done; }

    ac_dir_suffix="/`echo $ac_dir|sed 's,^\./,,'`"
    # A "../" for each directory in $ac_dir_suffix.
    ac_dots=`echo "$ac_dir_suffix" | sed 's,/[^/]*,../,g'`
  else
    ac_dir_suffix= ac_dots=
  fi

  case $srcdir in
  .)  ac_srcdir=.
      if test -z "$ac_dots"; then
         ac_top_srcdir=.
      else
         ac_top_srcdir=`echo $ac_dots | sed 's,/$,,'`
      fi ;;
  [\\/]* | ?:[\\/]* )
      ac_srcdir=$srcdir$ac_dir_suffix;
      ac_top_srcdir=$srcdir ;;
  *) # Relative path.
    ac_srcdir=$ac_dots$srcdir$ac_dir_suffix
    ac_top_srcdir=$ac_dots$srcdir ;;
  esac

  if test x"$ac_file" != x-; then
    { echo "$as_me:4323: creating $ac_file" >&5
echo "$as_me: creating $ac_file" >&6;}
    rm -f "$ac_file"
  fi
  # Let's still pretend it is `configure' which instantiates (i.e., don't
  # use $as_me), people would be surprised to read:
  #    /* config.h.  Generated automatically by config.status.  */
  configure_input="Generated automatically from `echo $ac_file_in |
                                                 sed 's,.*/,,'` by configure."

  # First look for the input files in the build tree, otherwise in the
  # src tree.
  ac_file_inputs=`IFS=:
    for f in $ac_file_in; do
      case $f in
      -) echo $tmp/stdin ;;
      [\\/$]*)
         # Absolute (can't be DOS-style, as IFS=:)
         test -f "$f" || { { echo "$as_me:4341: error: cannot find input file: $f" >&5
echo "$as_me: error: cannot find input file: $f" >&2;}
   { (exit 1); exit 1; }; }
         echo $f;;
      *) # Relative
         if test -f "$f"; then
           # Build tree
           echo $f
         elif test -f "$srcdir/$f"; then
           # Source tree
           echo $srcdir/$f
         else
           # /dev/null tree
           { { echo "$as_me:4354: error: cannot find input file: $f" >&5
echo "$as_me: error: cannot find input file: $f" >&2;}
   { (exit 1); exit 1; }; }
         fi;;
      esac
    done` || { (exit 1); exit 1; }
EOF
cat >>$CONFIG_STATUS <<EOF
  sed "$ac_vpsub
$extrasub
EOF
cat >>$CONFIG_STATUS <<\EOF
:t
/@[a-zA-Z_][a-zA-Z_0-9]*@/!b
s,@configure_input@,$configure_input,;t t
s,@srcdir@,$ac_srcdir,;t t
s,@top_srcdir@,$ac_top_srcdir,;t t
" $ac_file_inputs | (eval "$ac_sed_cmds") >$tmp/out
  rm -f $tmp/stdin
  if test x"$ac_file" != x-; then
    mv $tmp/out $ac_file
  else
    cat $tmp/out
    rm -f $tmp/out
  fi

done
EOF

cat >>$CONFIG_STATUS <<\EOF

{ (exit 0); exit 0; }
EOF
chmod +x $CONFIG_STATUS
ac_clean_files=$ac_clean_files_save

# configure is writing to config.log, and then calls config.status.
# config.status does its own redirection, appending to config.log.
# Unfortunately, on DOS this fails, as config.log is still kept open
# by configure, so config.status won't be able to write to it; its
# output is simply discarded.  So we exec the FD to /dev/null,
# effectively closing config.log, so it can be properly (re)opened and
# appended to by config.status.  When coming back to configure, we
# need to make the FD available again.
if test "$no_create" != yes; then
  ac_cs_success=:
  exec 5>/dev/null
  $SHELL $CONFIG_STATUS || ac_cs_success=false
  exec 5>>config.log
  # Use ||, not &&, to avoid exiting from the if with $? = 1, which
  # would make configure fail if this is the last instruction.
  $ac_cs_success || { (exit 1); exit 1; }
fi
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
  if test "x$GCC" != xyes; then
    # Check for BLAS in Sun Performance library:
    AC_CHECK_LIB(sunmath, acosp,
                 AC_CHECK_LIB(sunperf, $dgemm_func,
			      BLAS_LIBS="-xlic_lib=sunperf -lsunmath", ,
			      [-lsunmath $FLIBS]))
  fi
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
