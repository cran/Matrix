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
#include "lapackd.h"
#include "laexcp.h"

SVD::SVD(LaGenMatDouble& a, int nu = 0, int nvt = 0) :
    s(min(a.size(0), a.size(1))), u(), vt()
{
    int m = a.size(0), n = a.size(1);
    char jobu = 'N', jobvt = 'N';
    
    if (nu != 0 && nu != m && nu != n)
	throw(LaException("SVD : nu must be 0, or nrow(a), or ncol(a)"));
    if (nvt != 0 && nvt != m && nvt != n)
	throw(LaException("SVD : nv must be 0, or nrow(a), or ncol(a)"));
    if (nu >= m) { jobu = 'A'; }
    else if(nu >= min(m, n)) { jobu = 'S'; }
    if (nvt >= n) { jobvt = 'A'; }
    else if (nvt >= min(m, n)) { jobvt = 'S'; }

    u.resize(m, nu);
    vt.resize(nvt, n);

    LaGenMatDouble acopy(a);
    int lwork = 5 * max(m, n), info;
    VectorDouble work(lwork);
    F77_CALL(dgesvd)(jobu, jobvt, m, n, &acopy(0,0),
		     acopy.gdim(0), &s(0), &u(0,0), max(1, u.gdim(0)),
		     &vt(0,0), max(1, vt.gdim(0)), &work(0), lwork, info);
    if (info != 0)
	throw(LaException("SVD : dgesvd returned a non-zero info"));
}


