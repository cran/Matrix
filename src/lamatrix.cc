//
// Copyright (C) 2000-2000 the R Development Core Team
//

#include "lamatrix.h"
#include "eigen.h"

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

LaEigenDouble* LaMatDouble::eigen(bool leftEV, bool rightEV,
				  char balanc, char rcond)
{
    LaGenMatDouble tmp(size(0), size(1));
    tmp.inject(*this);
    return new LaGenEigenDouble(tmp, leftEV, rightEV, balanc, rcond);
}
