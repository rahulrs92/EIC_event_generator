#include<beam_properties.hh>
#include<cross_sections.hh>

using namespace std;

array<double,2> beam_properties::get_sigma()
{
    double T0 = kinetic_energy;
    double E = T0 + m0_p;
    double p = sqrt(E*E - m0_p*m0_p);
    
    double em_geo = emittance * m0_p / p;
    // return pow(em_geo * beta,0.5);

    sigma[0] = pow(em_geo * beta_function[0],0.5); //sig(BP.beta_function[0]);
    sigma[1] = pow(em_geo * beta_function[1],0.5); //sig(BP.beta_function[1]);

    return sigma;
}