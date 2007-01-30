#!/bin/sh
#
R=${R:-R-patched}
if [ x$R_HOME = x ] ; then R_HOME=`$R RHOME`; fi
#
MatrixDir=`dirname $0`; cd $MatrixDir; MatrixDir=`pwd`
if [ ! -d $MatrixDir ]
then echo "no directory '$MatrixDir' .. exiting"; exit 3
fi
cd $MatrixDir

RINC=${R_HOME}/include
#
out=DEPS.mkf_automade
gcc -I$RINC -MM *.c | perl -pe "s{$RINC/[^.]*.h( \\\\\\n)?}{}g" > $out
echo ''; echo "$0 done.  Resulting file is $MatrixDir/$out"
