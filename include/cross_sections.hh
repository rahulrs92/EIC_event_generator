#pragma once

#include<matrix3.hh>

#include<TMath.h>

const double m0_p = 0.938272; // GeV/c^2
const double rad = TMath::Pi() / 180.0;
const double deg = 1.0 / rad;
const double N_a = 6.02214076e23;

namespace Bonin
{
    const double m0_C = 11.17486; // GeV/c^2 
    
    auto dcs = [](double* x, double* par)
    {
        double T0 = par[0];             // Kinetic energy in GeV
        double E = T0 + m0_p;
        double p = sqrt(E*E - m0_p*m0_p);

        double th = x[0];               // theta in degrees

        double G = -2.98e-23    + 9.32e-23*p        - 2.39e-24*p*p;
        double B = -0.88*th     + 5.67e-02*th*th    - 1.39e-03*th*th*th;

        double result = (G * exp(B)) * 1e27;
        return result;                  // differential cross section in millibarns per steradian
    };

    auto f_theta = [](double* x, double* par)
    {
        double T0 = par[0];             // Kinetic energy in GeV
        double E = T0 + m0_p;
        double p = sqrt(E*E - m0_p*m0_p);

        double t = par[1];              // Thickness of target in m
        double a = par[2];              // Areal density of target in g/cm^2

        double th = x[0];               // theta in degrees

        double L = 35.3*p       - 1.58*p*p          - 2.62e-03*p*p*p;

        double result = 2*M_PI * sin(th * rad) * (a * N_a / 12) * exp(t/L * 100) * dcs(x,par) * 1e-27;
        return result;                  // differential count rate per degree (per particle incident on the target)
    };    

    auto Ay = [](double* x, double* par)
    {
        double T0 = par[0];             // Kinetic energy in GeV
        double E = T0 + m0_p;
        double p = sqrt(E*E - m0_p*m0_p);
        double pd = p - 1.4;

        double th = x[0];               // theta in degrees
        double r = p * sin(th * rad);

        double a = 1.626    - 1.786*pd      + 1.933*pd*pd   - 0.285*pd*pd*pd;
        double b = 14.909   - 0.594*pd      - 91.766*pd*pd  + 207.064*pd*pd*pd;
        double c = 574.719  - 560.410*pd    - 610.026*pd*pd + 1154.66*pd*pd*pd;
        double d = 0.1275   + 0.0550*pd     - 0.138*pd*pd   - 0.639*pd*pd*pd;

        double result = 1 + b*r*r + c*r*r*r*r;
        result = a*r / result;
        result = result + d*p*sin(5*th * rad);

        return result;
    };
    
    auto azim = [](double* x, double* par)
    {
        double th = x[0];       // theta in degrees
        double Px = par[1];     // polarization in x-direction
        double Py = par[2];     // polarization in y-direction
        double phi = x[1];      // phi in radians

        double result = dcs(x,par) * (1 + Ay(x, par) * (Py*cos(phi) - Px*sin(phi)));
        return result;
    };

    // Variables: theta, phi
    // Parameters: incident kinetic energy, target thickness, target areal density, Px, Py
    auto f_theta_phi = [](double* x, double* par)
    {
        double th = x[0];       // theta in degrees
        double Px = par[3];     // polarization in x-direction
        double Py = par[4];     // polarization in y-direction
        double phi = x[1];      // phi in radians

        double result = f_theta(x,par) * (1 + Ay(x, par) * (Py*cos(phi) - Px*sin(phi)));
        return result;
    };
}

namespace Distributions
{
    auto gaus2_standard = [](double* z, double* par)
    {
        double x = z[0], y = z[1];
        double r = pow(pow(x,2) + pow(y,2),0.5);

        double result = 1/(2*TMath::Pi()) * exp(-pow(r,2)/2);
        return result;
    };
}

namespace Lorentz
{
    const matrix<4,4> metric = {+1,0,0,0,
                                0,-1,0,0,
                                0,0,-1,0,
                                0,0,0,-1};
    
    
}