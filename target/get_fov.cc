#include<target_properties.hh>

#include<matrix3.hh>

std::array<double,5> target_properties::get_fov(const beam_properties& BP)
{
    if(!f_fov)
    {
        // Calculate the transformed target properties based on the loaded values
        auto sig = BP.get_sigma();
        double sig_x = sig[0];
        double sig_y = sig[1];

        matrix<2,1> t = {position[0], position[1]}, r = {cos(theta * rad), sin(theta * rad)};

        // Transformation matrix to go from symmetric to asymmetric frame
        S = {sig_x, 0.0, 0.0, sig_y};                               // matrix<2,2>
        
        auto rbS = r/S;

        // free parameter value corresponding to the foot of the perpendicular from the origin to the strip
        auto ku = - ((t/S).tp() * rbS) / (rbS.tp() * rbS);
        double nku = ku / rbS.norm();
        // cout << "normalized ku: " << nku << endl;                   // in symmetric frame units 

        // location of the foot of the perpendicular
        auto u = (t/S) + ku * rbS;
        // cout << "u.r': " << u.tp() * rbS << endl;                   // Should be zero if u is perpendicular to rbS

        // perpendicular distance of the strip from the center of the beam
        double d2 = u.norm();
        
        // unit vector perpendicular to strip in asymmetric frame
        matrix<2,1> rt = {-r[1], r[0]};

        // width of the strip in symmetric frame
        double w2 = d_eff * (rt/S).norm();

        // length of the strip in symmetric frame
        double l2 = l_t * rbS.norm();
        
        // antizlockwise angle (from +ve x-axis) of the strip in symmetric frame 
        double thS = atan2(rbS[1],rbS[0]);
        // cout << "rbS = " << rbS << endl;
        // cout << "Theta: " << thS * 180 / TMath::Pi() << "°" << endl;

        // Rotation matrix of strip in symmetric frame
        R = {cos(thS), -sin(thS), sin(thS), cos(thS)};              // matrix<2,2>

        auto fov_calc = [](double d, double w)
        {
            // 1. If w is relatively large, direct subtraction is cheap and perfectly safe.
            if (abs(w) >= 1e-5 * abs(d) || abs(d) < 1e-30) return 0.5 * (erf((d + w/2)/sqrt(2.0)) - erf((d - w/2)/sqrt(2.0)));

            // 2. If w is extremely small, use a truncated 3-term Taylor series to avoid catastrophic cancellation.
            double lead = pow(2 * TMath::Pi(), -0.5) * std::exp(-pow(d, 2) / 2.0);
            
            double term1 = w;
            double term2 = ((pow(d, 2) - 1) / 24) * pow(w, 3);
            double term3 = ((pow(d, 4) - 6*pow(d, 2) + 3) / 1920) * pow(w, 5);
            // At w < 1e-5, the 4th term (w^7) is always smaller than 1e-20 and can be safely ignored.

            return lead * (term1 + term2 + term3);
        };
        
        // return the coordinate range of the strip in the symmetric frame...
        fov[0] = -l2/2 + ku;    // x_min
        fov[1] = l2/2 + ku;     // x_max
        fov[2] = -d2 - w2/2;    // y_min
        fov[3] = -d2 + w2/2;    // y_max

        // ...and the calculated overlap fraction
        fov[4] = fov_calc(d2,w2) * fov_calc(ku,l2);

        f_fov = true;
    }
    return fov;
}