//   R : A Computer Language for Statistical Data Analysis
//   Copyright (C) 2000  the R Development Core Team

//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.

//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "SVD.h"

SVD::SVD(LaGenMatDouble& a, VectorDouble s, LaGenMatDouble& u,
	 LaGenMatDouble& v) : s(s)
{
    this->a.ref(a);
    this->u.ref(u);
    this->v.ref(v);
    doDecomp();
}

SVD::SVD(LaGenMatDouble& a, VectorDouble s) : s(s)
{
    this->a.ref(a);
    doDecomp();
}

SVD::SVD(LaGenMatDouble& a) : s(min(a.size(0), a.size(1)))
{
    this->a.ref(a);
    doDecomp();
}

void SVD::doDecomp()
{
    int lwork = 5 * max(a.size(0), a.size(1));
    double *work = new double[lwork];
    char jobu = 'N', jobvt = 'N';

    if (u.size(0) >= a.size(0)) {
	if (u.size(1) >= a.size(0)) { jobu = 'A'; }
	else if(u.size(1) >= min(a.size(0), a.size(1))) { jobu = 'S'; }
    }
    if (v.size(0) >= a.size(1)) {
	if (v.size(1) >= a.size(0)) { jobvt = 'A'; }
	else if (v.size(1) >= min(a.size(0), a.size(1))) { jobvt = 'S'; }
    }
  
    F77_CALL(dgesvd)(jobu, jobvt, a.size(0), a.size(1), &a(0,0),
		     a.gdim(0), s.addr(), &u(0,0), max(1, u.gdim(0)),
		     &v(0,0), max(1, v.gdim(0)), work, lwork, info); 
    delete[] work;
}


