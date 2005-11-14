CHOLMOD: a sparse CHOLesky MODification package 
-----------------------------------------------

Version 0.6, Aug 30, 2005.  Copyright (c) 2005.

Some Modules of CHOLMOD are copyrighted by the University of Florida.  Most are
copyrighted by the authors: Timothy A. Davis (all of them), and William W. Hager
(the Modify Module).

CHOLMOD relies on several other packages:  AMD, COLAMD, CCOLAMD, METIS,
the BLAS, and LAPACK.

AMD is authored by T. Davis, Iain Duff, and Patrick Amestoy.
COLAMD is authored by T. Davis and Stefan Larimore, with algorithmic design
in collaboration with John Gilbert and Esmond Ng.
CCOLAMD is authored by T. Davis and Siva Rajamanickam.

LAPACK and the BLAS are authored by Jack Dongarra and many others.
LAPACK is available at http://www.netlib.org/lapack

METIS is authored by George Karypis, Univ. of Minnesota.  Its use in CHOLMOD
is optional.  See http://www-users.cs.umn.edu/~karypis/metis.
Place a copy of the metis-4.0 directory in the same directory that
contains the CHOLMOD, AMD, COLAMD, and CCOLAMD directories.  Then do:
cd metis-4.0 ; make
You may then compile CHOLMOD.

The CHOLMOD, AMD, COLAMD, CCOLAMD, and UFconfig directories must all reside
in a common parent directory.  To compile all these libraries,
edit UFconfig/UFconfig.mk to reflect your environment (C compiler, location
of the BLAS, and so on) and then type "make" in either the CHOLMOD directory
or in the parent directory of CHOLMOD.  See each package for more details on
how to compile them.

For use in MATLAB (on any system, including Windows):  start MATLAB,
cd to the CHOLMOD/MATLAB directory, and type cholmod_make in the MATLAB
Command Window.

On the Pentium, do NOT use the Intel MKL BLAS with CHOLMOD.  It has a bug in
dgemm when computing A*B'.  The bug generates a NaN result, when the inputs
are well-defined.  Use the Goto BLAS instead.  It's faster and more reliable.
See http://www.tacc.utexas.edu/~kgoto/ or
http://www.cs.utexas.edu/users/flame/goto/.
Sadly, the Intel MKL BLAS is the default for MATLAB 7.0.4.  See
http://www.mathworks.com/support/bugreports/details.html?rp=252103 for more
details.  To workaround this problem on Linux, set environment variable
BLAS_VERSION to libmkl_p3.so:libguide.so. On Windows, set environment variable
BLAS_VERSION to mkl_p3.dll.

Acknowledgements:  this work was supported in part by the National Science
Foundation (NFS CCR-0203270 and DMS-9803599), and a grant from Sandia National
Laboratories (Dept. of Energy) which supported the development of CHOLMOD's
Partition Module.
