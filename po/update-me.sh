#!/bin/sh
#
## Script for updating package-specific *.pot files
## written such that it should work for any package
#
thisdir=`dirname $0` ; cd $thisdir; thisdir=`pwd`
echo 'preliminary thisdir='$thisdir
pkgDIR=`dirname $thisdir`
pkg=`basename $pkgDIR`
echo '  -->        pkgDIR='$pkgDIR' ; pkg='$pkg
cd `R-devel RHOME`/po
make pkg-update PKG=$pkg PKGDIR=$pkgDIR