#!/bin/sh

ufl_URL=http://www.cise.ufl.edu/research/sparse
wget $ufl_URL/amd/current/AMD.tar.gz
wget $ufl_URL/cholmod/current/CHOLMOD.tar.gz
wget $ufl_URL/colamd/current/COLAMD.tar.gz
wget $ufl_URL/UFconfig/current/UFconfig.tar.gz
  ## install UFconfig.h file (now needed by some UFsparse libraries)
tar zxf UFconfig.tar.gz UFconfig/UFconfig.h UFconfig/README.txt
  ## Move the UFconfig/README.txt file to ../inst/doc/UFsparse/UFconfig.txt
mv UFconfig/README.txt ../inst/doc/UFsparse/UFconfig.txt
  ## touch the file UFconfig/UFconfig.mk.  We use other configuration
  ## environment variables but this name is embedded in some Makefiles
touch UFconfig/UFconfig.mk
  ## install COLAMD/Source and COLAMD/Include directories
tar zxf COLAMD.tar.gz COLAMD/Source/ COLAMD/Include/
  ## install documentation for COLAMD
tar zxf COLAMD.tar.gz COLAMD/README.txt COLAMD/ChangeLog
mv COLAMD/README.txt ../inst/doc/UFsparse/COLAMD.txt
  ## install AMD/Source and AMD/Include directories
tar zxf AMD.tar.gz AMD/Source AMD/Include AMD/README.txt
  ## restore the AMD/Source/Makefile
svn revert AMD/Source/Makefile
  ## install AMD documentation
mv AMD/README.txt ../inst/doc/UFsparse/AMD.txt
  ## remove Fortran source files and GNUMakefile
rm AMD/Source/*.f AMD/Source/GNUmakefile
  ## install CHOLMOD source files
for d in Check Cholesky Core Include Lib MatrixOps Modify Partition Supernodal
do
    tar zxf ./CHOLMOD.tar.gz CHOLMOD/$d
done
  ## install CHOLMOD documentation in ../inst/doc/UFsparse
tar zxf ./CHOLMOD.tar.gz CHOLMOD/README.txt
mv CHOLMOD/README.txt ../inst/doc/UFsparse/CHOLMOD.txt

mv CHOLMOD/Lib/Makefile CHOLMOD/Lib/Makefile_CHOLMOD
svn revert CHOLMOD/Lib/Makefile
  ##

ls -l CHOLMOD/Lib/Makefile_CHOLMOD
echo 'now diff CHOLMOD/Lib/Makefile with .... Makefile_CHOLMD'
echo ' make changes as necessary, and then'
echo ' rm CHOLMOD/Lib/Makefile_CHOLMOD'

  ## remove the downloaded tar files
rm CHOLMOD.tar.gz AMD.tar.gz COLAMD.tar.gz UFconfig.tar.gz
