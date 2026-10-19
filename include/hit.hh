#pragma once

#include<matrix3.hh>

class hit
{
protected:
    int ID{0};
    matrix<3,1> loc{0.0, 0.0, 0.0};
    matrix<3,1> x{0.0, 0.0, 0.0};
    long eID{0};

public:
    hit() {}
    hit(int iID, matrix<3,1> iloc, matrix<3,1> ix, long eID): ID(iID), loc(iloc), x(ix) {}

    operator bool () { return ID!=0; }
    const matrix<3,1>& get_loc() const { return loc; }
    const matrix<3,1>& get_x() const { return x; }
    const long& get_eID() const { return eID; }
    const int& get_ID() const { return ID; }
};