#!/bin/sh
#
R=${R:-R-patched}
if [ x$R_HOME = x ] ; then R_HOME=`$R RHOME`; fi
RINC=${R_HOME}/include
# For consistency of *.c below: set locale to simple C
# MM: Setting locale here does *NOT* help me; setting in the calling shell does !?!?
# Do something like this in the shell
LC_CTYPE=C export LC_CTYPE
unset LANG
unset LANGUAGE
unset LC_ADDRESS
unset LC_COLLATE
unset LC_IDENTIFICATION
unset LC_MEASUREMENT
unset LC_MESSAGES
unset LC_NAME
unset LC_PAPER
unset LC_TELEPHONE
unset LC_NUMERIC
unset LC_MONETARY
unset LC_TIME
## for debugging:
# echo '------------------------'
# locale
# echo '------------------------'
#
MatrixDir=`dirname $0`/..; cd $MatrixDir; MatrixDir=`pwd`
if [ ! -d $MatrixDir ]
then echo "no directory '$MatrixDir' .. exiting"; exit 3
fi
cd $MatrixDir
## __begin fix__   CHOLMOD has one include for which gcc -MM fails below:
FIX=CHOLMOD/Include/cholmod.h
if [ -f $FIX ]
then
  sed '/^#include "SuiteSparse_config/s/\(.*\)/\/* \1 *\//' $FIX > ${FIX}_fixed
  mv $FIX ${FIX}_orig
  mv ${FIX}_fixed $FIX
fi
## __end fix__
out=scripts/DEPS.mkf_automade
(echo '#-*- Makefile -*-
#-------------'" produced by $0  (plus minimal emacs cleanup)
#"
 ls *.c | grep -v '^t_' | xargs gcc -I$RINC -MM |
     perl -pe "s{$RINC/[^.]*.h( \\\\\\n)?}{}g; s{^[ \\s]+Syms.h \\\\\\n[ \\s]*}{ Syms.h }"
) > $out
if [ -f ${FIX}_orig ] ; then mv ${FIX}_orig $FIX ; fi
echo '------------------------'
echo ''; echo "$0 done.  Resulting file is
 $MatrixDir/$out" ; echo

