#include<beam_properties.hh>

#include<matrix3.hh>

#include "TRandom.h"

std::array<double, 2> beam_properties::spin_sample()
{
    matrix<2,1> p = {polarization[0],polarization[1]};

    auto r = matrix<2,1>{p.norm(), 1} / matrix<2,2>{1,-1,1,1};

    double u1 = gRandom->Rndm();
    
    return (u1 > r[1]) ? std::array<double, 2>{p[0],p[1]} : std::array<double, 2>{-p[0],-p[1]};
}