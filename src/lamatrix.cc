//
// Copyright (C) 2000-2000, 2002 the R Development Core Team
//

#include "lamatrix.h"
#include "eigen.h"

int* LaMatrix::info_= new int;  // turn off info print flag.

std::ostream& LaMatrix::Info(std::ostream& s)
{
    s << "Size: (" << size(0) << "x" << size(1) << ") " ;
    s << "Indices: " << index(0) << " " << index(1);
    return s;
}

std::ostream& LaMatDouble::Info(std::ostream& s)
{
    LaMatrix::Info(s);
    s << "addr: " << addr() << std::endl;
    return s;
}

std::ostream& LaMatInt::Info(std::ostream& s)
{
    LaMatrix::Info(s);
    s << "addr: " << addr() << std::endl;
    return s;
}

LaEigenDouble* LaMatDouble::eigen(bool leftEV, bool rightEV,
				  char balanc, char rcond)
{
    LaGenMatDouble tmp(size(0), size(1));
    tmp.inject(*this);
    return new LaGenEigenDouble(tmp, leftEV, rightEV, balanc, rcond);
}
