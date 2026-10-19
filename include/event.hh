#pragma once

#include<matrix3.hh>

#include<cross_sections.hh>
#include<beam_properties.hh>
// #include<target_properties.hh>

class event
{
protected:
    double q;

    matrix<3,1> vtx;
    matrix<4,1> p;

    long eID;

public:
    event() = delete;
    event(long ieID, double iq, matrix<3,1> vertex, matrix<4,1> four_momentum): eID(ieID), q(iq), vtx(vertex), p(four_momentum) {}

    const matrix<3,1>& vertex() const { return vtx; }
    const matrix<4,1>& four_momentum() const { return p; }
    const long& evID() const { return eID; } 

    double rest_mass() const { return pow(p.tp() * Lorentz::metric * p, 0.5); }
    double KE() const { return p[0] - rest_mass(); }
    double momentum() const { return p.sub<3,1>(1,0).norm(); }
    double charge() const { return q; }
};