#include<detector_properties.hh>

using namespace std;

matrix<N_R, N_R> detector_properties::cov_geom(const matrix<2, 1>& Ab0, const target_properties& TP, const beam_properties& BP)
{
    matrix<N_R, N_R> C;
    
    auto bs = BP.get_sigma();
    // double d_arc = atan2(sz_strip_r[0] / 2, F[-1].norm()); // Semi-azhimuthal coverage of the strip in the recoil detector

    //Variance of the beam along the strip
    double varpl = pow(bs[0] * bs[1], 2) / (pow(bs[0]*sin(TP.theta*rad), 2) + pow(bs[1]*cos(TP.theta*rad), 2));
    // Variance of the beam perpendicular to the strip
    double varpp = pow(TP.d_eff, 2) / 12; // Assuming uniform distribution across the strip width

    for(int i = 0; i < N_R; ++i)
    {
        // Target Angular Variance
        double phi = atan2(F[-(i+1)][1], F[-(i+1)][0]); // Azimuthal angle of the detector
        double angle = TP.theta*rad - phi;
        double vart = (varpl * pow(sin(angle), 2) + varpp * pow(cos(angle), 2)) / (F[-(i+1)].tp() * F[-(i+1)]);
        
        // Total Variance = Target Angular Variance + Detector Arc Variance
        // double var = vart + pow(d_arc,2) / 3; // Assuming uniform distribution across the detector arc
        
        // Assign initial asymmetry scaled by the total variance to the diagonal element of the covariance matrix
        C[i][i] = pow(Ab0[0]*cos(phi) - Ab0[1]*sin(phi), 2) * vart; 
    }

    return C;
}