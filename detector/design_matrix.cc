#include<detector_properties.hh>

using namespace std;

matrix<N_R, 2> detector_properties::design_matrix(const target_properties& TP, const beam_properties& BP)
{
    matrix<N_R, 2> R;
    
    auto func = [&](double* x, double* p)
    {
        matrix<2, 1> pos = {x[0], x[1]};        // position vector of point from beam centre in symmetric frame

        matrix<2, 1> ff = {p[0], p[1]};         // position vector of detector face centre in asymmetric frame
        // matrix<2, 1> pos = {x[0], x[1]};
        
        // Transform the point from the symmetric frame to the asymmetric frame
        matrix<2, 1> pos_asym = TP.S * TP.R * pos;  // position vector of point from beam centre in asymmetric frame
        // Calculate the position vector of the detector face centre from the point
        ff = (ff - pos_asym); ff /= ff.norm();

        bool p2 = p[2] > 0; // Determine if we are calculating the cosine or sine component based on the third parameter
        return TP.gaus2(pos[0], pos[1]) * (p2 ? ff[0] : ff[1]); // Return the weighted value based on the cosine or sine component
    };

    auto fov = TP.get_fov();
    // Extract strip properties {l2, ku, d2, w2} in symmetric frame from limits of integration
    auto stp = matrix<4, 1>{fov[0], fov[1], fov[2], fov[3]} / matrix<4,4>{
       -0.5,    1,  0,  0, 
        0.5,    1,  0,  0, 
        0,      0, -1, -0.5, 
        0,      0, -1,  0.5  };

    // TP.gaus2.set_range(fov[0], fov[1], fov[2], fov[3]);
    TF2* wavg_cos = new TF2("wavg_cos", func, fov[0], fov[1], fov[2], fov[3], 3);
    wavg_cos->FixParameter(2, 1.0); // Set the third parameter to 1.0 for the cosine calculation
    TF2* wavg_sin = new TF2("wavg_sin", func, fov[0], fov[1], fov[2], fov[3], 3);
    wavg_sin->FixParameter(2, -1.0); // Set the third parameter to -1.0 for the sine calculation

    double D = fov[4]; // Strip integral in the symmetric frame
    double d_arc = atan2(sz_r[0] / 2, F[-1].norm()); // Semi-azhimuthal coverage of the strip in the recoil detector

    for(int i = 0; i < N_R; ++i)
    {
        wavg_cos->SetParameter(0, F[-(i+1)][0]); // Set the x-coordinate of the recoil detector face centre
        wavg_cos->SetParameter(1, F[-(i+1)][1]); // Set the y-coordinate of the recoil detector face centre

        wavg_sin->SetParameter(0, F[-(i+1)][0]); // Set the x-coordinate of the recoil detector face centre
        wavg_sin->SetParameter(1, F[-(i+1)][1]); // Set the y-coordinate of the recoil detector face centre

        // 1. Calculate integrals of the function in the symmetric frame
        // 2. Scale by sines/cosines in the asymmetric frame
        // 3. Normalize by the symmetric frame strip integral
        R[i][1] = wavg_cos->Integral(fov[0], fov[1], fov[2], fov[3]) / D; // Calculate the cosine component of the design matrix
        R[i][0] = wavg_sin->Integral(fov[0], fov[1], fov[2], fov[3]) / D; // Calculate the sine component of the design matrix

        R[i] *= sin(d_arc) / d_arc; // Scale the row by sinc of the strip's angular coverage
    }

    delete wavg_cos;
    delete wavg_sin;

    return R;
}