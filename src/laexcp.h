// -*- c++ -*-
//              LAPACK++ 1.1 Linear Algebra Package 1.1
//               University of Tennessee, Knoxvilee, TN.
//            Oak Ridge National Laboratory, Oak Ridge, TN.
//        Authors: J. J. Dongarra, E. Greaser, R. Pozo, D. Walker
//                 (C) 1992-1996 All Rights Reserved
//
//                             NOTICE
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Institutions (University of Tennessee, and Oak Ridge National
// Laboratory) nor the Authors make any representations about the suitability 
// of this software for any purpose.  This software is provided ``as is'' 
// without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.
//
// Modifications Copyright (C) 2000-2001 the R Development Core Team

#ifndef _LA_EXCEPTION_H_
#define _LA_EXCEPTION_H_

#ifdef length
#undef length   // override possible definition in Rinternals.h
#endif
#include <string>

class LaException {
    std::string s;
public:
    LaException(const std::string& s) : s(s) { };
    LaException(const std::string& where, const std::string& what) : s(where + what) { };
    LaException(const char* c) : s(c) { };
    LaException(const char* where, const char* what) : s(std::string(where)
	+ std::string(what)) { };

    const char* what() { return s.c_str(); }
};

#endif  

