// -*- c++ -*-
//
//   R : A Computer Language for Statistical Data Analysis
//   Copyright (C) 2000, 2002  the R Development Core Team

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

#ifndef _ORTHONORMAL_H_
#define _ORTHONORMAL_H_

#include "lafnames.h"
#include LA_GEN_MAT_DOUBLE_H

class LaColOrthogonalMatDouble : public LaGenMatDouble
{
 public:
    SEXP asSEXP() const;
    LaColOrthogonalMatDouble& operator=(double s)
    {
        LaGenMatDouble::operator=(s);
        return *this;
    }
};

class LaRowOrthogonalMatDouble : public LaGenMatDouble
{
 public:
    SEXP asSEXP() const;
    LaRowOrthogonalMatDouble& operator=(double s)
    {
        LaGenMatDouble::operator=(s);
        return *this;
    }
};

class LaColOrthonormalMatDouble : public LaColOrthogonalMatDouble
{
 public:
    SEXP asSEXP() const;
    LaColOrthonormalMatDouble& operator=(double s)
    {
        LaGenMatDouble::operator=(s);
        return *this;
    }
};

class LaRowOrthonormalMatDouble : public LaRowOrthogonalMatDouble
{
 public:
    SEXP asSEXP() const;
    LaRowOrthonormalMatDouble& operator=(double s)
    {
        LaGenMatDouble::operator=(s);
        return *this;
    }
};

class LaOrthogonalMatDouble : public LaColOrthonormalMatDouble
{
 public:
    SEXP asSEXP() const;
    LaOrthogonalMatDouble& operator=(double s)
    {
        LaGenMatDouble::operator=(s);
        return *this;
    }
};

#endif // _ORTHONORMAL_H_
