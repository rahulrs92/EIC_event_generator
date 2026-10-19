#include<target_properties.hh>

#include<TF2.h>

double target_properties::get_efficiency(beam_properties& BP)
{
    // TF2* pdf = new TF2("pdf", Bonin::f_theta_phi, 0.0, 30.0, -M_PI, M_PI, 5);
    pdf->SetParameters(BP.kinetic_energy, d_eff, d_t, BP.polarization[0], BP.polarization[1]);
    return pdf->Integral(0.0, 30.0, -M_PI, M_PI);
}