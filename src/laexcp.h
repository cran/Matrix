//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.

//
//  C++ exception handling is not currently supported by most compilers.
//  The macros below allow one to "mimic" the throw expression of
//  character strings.  Transfer of control does not correspond
//  to a "try" block, but rather exits the program with the proper
//  banner.  This is a similar behavior to the real exception handling
//  if there was no explicit user-supplied try block.
//

#ifndef _LA_EXCEPTION_H_
#define _LA_EXCEPTION_H_

#ifdef length
#undef length   // override possible definition in Rinternals.h
#endif
#include <string>

//#include <iostream.h>
//#include <stdlib.h>

//#define LaException(where, what)    where , what
//#define throw throw_
//inline void throw(const char *where, const char *what)
//{
//    cerr << "Exception: " << where << "  " << what << endl;
//    exit(1);
//}

//#define LaException(where, what) where

class LaException {
    string s;
public:
    LaException(const string& s) : s(s) { };
    LaException(const string& where, const string& what) : s(where + what) { };
    LaException(const char* c) : s(c) { };
    LaException(const char* where, const char* what) : s(string(where)
	+ string(what)) { };

    const char* what() { return s.c_str(); }
};


#endif  

