//
// Copyright (C) 2000-2000 the R Development Core Team
//

#include "lamatrix.h"

int LaMatrix::debug_ = 0;	// turn off global deubg flag initially.
                                // use A.debug(1) to turn on/off,
                                // and A.debug() to check current status.

int* LaMatrix::info_= new int;  // turn off info print flag.

ostream& LaMatrix::Info(ostream& s)
{
    s << "Size: (" << size(0) << "x" << size(1) << ") " ;
    s << "Indices: " << index(0) << " " << index(1);
    return s;
}

ostream& LaMatDouble::Info(ostream& s)
{
    LaMatrix::Info(s);
    s << "addr: " << (unsigned) addr() << endl;
    return s;
}

ostream& LaMatInt::Info(ostream& s)
{
    LaMatrix::Info(s);
    s << "addr: " << (unsigned) addr() << endl;
    return s;
}
