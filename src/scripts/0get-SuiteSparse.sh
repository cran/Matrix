#!/bin/sh
## This *REPLACES* former ./UFsparse_download.sh
## Update Libraries from Tim Davis (formerly University of Florida (UF), now Texas A&M)'s  "SuiteSparse":
#
if [ ! -d ../src -o  ! -d ./scripts ]
then echo 'Must run in Matrix/src/ !' ; exit 1
fi
getSPQR=no
##     --- since late summer 2010, we no longer get SPQR
#
# Tim Davis moved to Texas A&M, on July 1, 2014
# ufl_URL=http://www.cise.ufl.edu/research/sparse/SuiteSparse/current/
# TGZ=SuiteSparse.tar.gz
#  wget -nc  $ufl_URL/$TGZ
## 2020-04-03: Pointed to from  http://faculty.cse.tamu.edu/davis/suitesparse.html
## ----------- SuiteSparse is on github currently:
GH_base=https://github.com/DrTimothyAldenDavis/SuiteSparse
GH_rel=${GH_base}/releases
GH_latest=${GH_rel}/latest
curl --dump-header hd_latest $GH_latest > curl_out 2>&1
# Careful: curl gives results "in MSDOS Format" with \cr\lf --> remove \r
VER=$(sed -n -e '/^location:/{p;q}' hd_latest | tee SS_location | sed 's#.*/v\([1-9]\.[0-9]\.[0-9]\)\r#\1#')
echo "SS_location:"; cat SS_location
echo "
  SuiteSparse version VER='$VER'
"
GH_tar_url=https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/v${VER}.tar.gz
SS=SuiteSparse-${VER}
TGZ=${SS}.tar.gz
if [ -f $TGZ ]
then
    echo 'Tarfile present; not downloading (remove it to change this!)'
    echo 'Maybe *continue* downloading by'; echo;
    echo "	  wget -c $GH_tar_url -O $TGZ" ; echo
## just for experimenting!
#    echo ' >> INTERRUPT (Ctrl C) within 7 sec !) if you want do that '
#    sleep 7
else
    echo '  ==> Trying to get it from '"$GH_tar_url :"
    wget -nc $GH_tar_url -O $TGZ
fi
ls -l $TGZ

SSdocDir=../inst/doc/SuiteSparse

## 1) SuiteSparse_config ---------------------------------------------
## NOTA BENE:  SuiteSparse_config/  is what  UFconfig/  used to be
Sdir=$SS/SuiteSparse_config
  ## install SuiteSparse_config.h file (now needed by some SuiteSparse libraries)
tar zxf $TGZ $Sdir/SuiteSparse_config.h $Sdir/SuiteSparse_config.c $Sdir/Makefile $Sdir/README.txt
  ## Move the SuiteSparse_config/README.txt file to docs:
mv $Sdir/README.txt $SSdocDir/SuiteSparse_config.txt
mv $Sdir/Makefile $Sdir/Makefile_orig
  ## touch the file $Sdir/SuiteSparse_config.mk.  We use other configuration
  ## environment variables but this name is embedded in some Makefiles
touch $Sdir/SuiteSparse_config.mk
  ## move directory *up*
dd=`basename $Sdir`; mv $Sdir/* $dd/
  ## Need to add the Matrix-specific changes to SuiteSparse_config/SuiteSparse_config.h :
## 2014-08: different patch:
patch -p0 < scripts/SuiteSparse_config.patch

## 2) COLAMD -----------------------------------------------
Sdir=$SS/COLAMD
   ## install COLAMD/Source and COLAMD/Include directories
tar zxf $TGZ $Sdir/Source/ $Sdir/Include/ $Sdir/Doc/ $Sdir/README.txt
## MM {2014-12}: following Makefile no longer exists
f=$Sdir/Source/Makefile
if [ -f $f ]
then Rscript --vanilla -e 'source("scripts/fixup-fn.R")' -e 'fixup("'$f'")'
fi
  ## install documentation for COLAMD
mv $Sdir/README.txt $SSdocDir/COLAMD.txt
mv $Sdir/Doc/ChangeLog $SSdocDir/COLAMD-ChangeLog.txt
rm -rf $Sdir/Doc
# 2014: no longer
# patch -p0 < scripts/COLAMD.patch
##          ---------------------
  ## move directory *up*
dd=`basename $Sdir`; rsync -auv $Sdir/ $dd/


## 3) AMD --------------------------------------------------
Sdir=$SS/AMD
  ## install AMD/Source, AMD/Include and AMD/Lib directories
tar zxf $TGZ $Sdir/Source $Sdir/Include $Sdir/Lib $Sdir/README.txt
  ## install AMD documentation
mv $Sdir/README.txt $SSdocDir/AMD.txt
  ## remove Fortran source files and GNUMakefile
rm $Sdir/Source/*.f $Sdir/Lib/GNUmakefile
#(for f in $Sdir/Include/amd_internal.h $Sdir/Source/amd_global.c; do diff -ubBw ${f}.~1~ $f ; done ) | tee scripts/AMD-noprint.patch
## 2014: no longer
# patch -p0 < scripts/AMD-noprint.patch
##          ---------------------
  ## move directory *up*
dd=`basename $Sdir`; rsync -auv $Sdir/ $dd/


## 4) CHOLMOD ----------------------------------------------
Sdir=$SS/CHOLMOD
  ## install CHOLMOD source files
for d in Check Cholesky Core Include Lib MatrixOps Modify Partition Supernodal README.txt
do
    tar zxf $TGZ $Sdir/$d
done
  ## install CHOLMOD documentation
mv $Sdir/README.txt $SSdocDir/CHOLMOD.txt

cp -p $Sdir/Lib/Makefile $Sdir/Lib/Makefile_CHOLMOD
Rscript --vanilla -e 'source("scripts/fixup-fn.R")' -e 'fixup("'$Sdir'/Lib/Makefile")'
## but typically, this is not good enough, so we need manual work:
mv $Sdir/Lib/Makefile $Sdir/Lib/Makefile_pre
  ## move directory *up*
dd=`basename $Sdir`; rsync -auv $Sdir/ $dd/

echo 'If there changes in the following you  ** MUST **  manually update
  <Matrix>/ inst/include/cholmod.h  --- to export what we have.

Also, RcppEigen headers may also need to be updated -- ask Doug Bates.
This can be VERY IMPORTANT,  not the least for lme4
'
svn diff --diff-cmd /usr/bin/diff -x "-bBw"  $dd/Include/cholmod_core.h
echo 'Did the above show any non trivial diffs? --> do update inst/include/cholmod.h !!
'

svn revert $dd/Lib/Makefile
ls -l $dd/Lib/Makefile_pre
echo "now   diff $dd/Lib/Makefile $dd/Lib/Makefile_pre  [I do it for you below]"
echo ' make changes as necessary, and then (later)'
echo " rm $dd"'/Lib/Makefile_*' ; echo
echo "Ok, now   diff $dd/Lib/Makefile $dd/Lib/Makefile_pre :"
diff $dd/Lib/Makefile $dd/Lib/Makefile_pre

## 5) CSparse -------------------------------------------------
Sdir=$SS/CSparse
  ## install CSparse/Source & CSparse/Include
tar zxf $TGZ $Sdir/Source $Sdir/Include $Sdir/README.txt

  ## Include:
echo -n "Moving from $Sdir/Include .. "
f=$Sdir/Include/cs.h
chmod a+r $f && mv $f .

  ## Source:
MatrixDir=`pwd`
cd $Sdir/Source
cat cs_*.c | sed -e '1 p' -e '/^#include/d' -e 's/\bprintf/Rprintf/g' > $MatrixDir/cs.c
cd $MatrixDir
patch -p0 < scripts/cs.patch
echo '[Ok]'
echo -n "removing $Sdir .."
rm -rf $Sdir
echo '[Ok]'



## 6) SPQR -------------------------------------------------
if [ $getSPQR = yes ]
then
  ## install SPQR source files
  for d in Source Include Lib
  do
      tar zxf ./SPQR.tar.gz SPQR/$d
  done
    ## install CHOLMOD documentation in ../inst/doc/UFsparse
  tar zxf ./SPQR.tar.gz SPQR/README.txt
  mv SPQR/README.txt $SSdocDir/SPQR.txt
    ## patch for Matrix:
  patch -p0 < scripts/SPQR.patch

  cp -p SPQR/Lib/Makefile SPQR/Lib/Makefile_SPQR
  Rscript --vanilla -e 'source("scripts/fixup-fn.R")' -e 'fixup("SPQR/Lib/Makefile")'
  mv    SPQR/Lib/Makefile SPQR/Lib/Makefile_pre
  svn revert SPQR/Lib/Makefile
    ##
  ls -l SPQR/Lib/Makefile_pre
  echo 'now diff SPQR/Lib/Makefile with SPQR/Lib/Makefile_pre'
  echo ' make changes as necessary, and then (later)'
  echo ' rm SPQR/Lib/Makefile_*' ; echo
fi

## ----- remove the downloaded tar file -------------------
echo 'You could (eventually) do

      rm '"$TGZ"'

but keeping it will not download it anew, if not changed.
Further, consider updating the doxygen docu on the R-forge web site via

 /u/maechler/R/Pkgs/Matrix-doxygen-update
'



