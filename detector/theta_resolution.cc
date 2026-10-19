#include<detector_properties.hh>

matrix<2,1,valarray<double>> detector_properties::theta_resolution(const target_properties& TP, const beam_properties& BP)
{
    matrix<2,1,valarray<double>> res = {valarray<double>(0.0,n_fdet),valarray<double>(0.0,n_rdet)};
    
    auto bs = BP.get_sigma();
    
    // Variance of the beam along the strip
    double varpl = pow(bs[0] * bs[1], 2) / (pow(bs[0]*sin(TP.theta*rad), 2) + pow(bs[1]*cos(TP.theta*rad), 2));
    // Variance of the beam perpendicular to the strip
    double varpp = pow(TP.d_eff, 2) / 12; // Assuming uniform distribution across the strip width
    
    for (int i = -static_cast<int>(n_rdet); i < static_cast<int>(n_fdet); i++)
    {
        int j = i>=0 ? i+1 : i;         // Current detector index
        auto sz_pix = (i>=0 ? sz_strip_f : sz_strip_r);
        double d_arc = atan2(sz_pix[1] / 2, F[j].norm()); // Semi-polar angular coverage of a strip in the detector
        
        // Target Angular Variance
        double phi = atan2(F[j][1], F[j][0]); // Azimuthal angle of the detector
        double cth = F[j][2] / F[j].norm(); // Cosine of polar angle of the detector
        double angle = TP.theta*rad - phi;
        double vart = pow(cth, 2) * (varpp * pow(sin(angle), 2) + varpl * pow(cos(angle), 2)) / (F[j].tp() * F[j]);
        
        // Total Variance = Detector Arc Variance + Target Angular Variance
        double var = vart + pow(d_arc,2) / 3; // Assuming uniform distribution across the strip arc

        // Assign deviations to resolution matrix
        res[i>=0 ? 0 : 1][i>=0 ? i : -i-1] = sqrt(var);
    }
    return res;
}