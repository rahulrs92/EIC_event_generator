#pragma once

#include<map>

#include<cross_sections.hh>

#include<beam_properties.hh>
#include<target_properties.hh>
#include<event.hh>
#include<hit.hh>

#include<matrix3.hh>

constexpr std::size_t N_F = N_FDET; // from input/detector.in
constexpr std::size_t N_R = N_RDET; // from input/detector.in

// class hit;

// /**
//  * @brief Struct to hold extracted flux, individual asymmetries, and full covariance matrix.
//  * 
//  * @tparam M Number of azimuthal detectors in the polarimeter array.
//  */
// template<size_t M>
// struct PolarimeterMetrics
// {
//     double L;               // Relative beam flux ratio: L = lambda_1 / lambda_2
//     double L_e;             // 1-sigma uncertainty on relative flux ratio (sigma_L)
//     matrix<M, 1> A;         // Column vector of exact detector asymmetries A_i (M x 1)
//     matrix<M, M> cov_A;     // Full M x M covariance matrix cov(A_i, A_j) including flux correlations
//     matrix<2, 2> cov_Ab;    // 2 x 2 covariance matrix cov(A_bs, A_bc) of physical beam parameters
//     int status = 0;         // Status code: 0 = success, -1 = non-positive counts/errors
// };

class detector_properties : public config_properties
{
public:
    // // Default initializations
    array<double,3> sz_f{0.0,0.0,0.0};                      // (length (azhimuthal), height (polar), thickness) in m
    array<double,3> pos_f{0.0,0.0,0.0};                     // (x, y, z) of face centre in m
    double tilt_f{0.0};                                     // in degrees

    size_t n_fdet{0};                                       // number of detectors
    array<size_t,2> n_fstrip{0,0};                          // number of strips

    array<double,3> sz_r{0.0,0.0,0.0};                      // (length (azhimuthal), height (polar), thickness) in m
    array<double,3> pos_r{0.0,0.0,0.0};                     // (x, y, z) of face centre in m
    double tilt_r{0.0};                                     // in degrees

    size_t n_rdet{0};                                       // number of detectors
    array<size_t,2> n_rstrip{0,0};                          // number of strips (theta-aligned, phi-aligned)

    // // Calculated properties

    // Size of strip in each detector (length along phi, length along theta) in m
    array<double,2> sz_strip_f{0.0,0.0}, sz_strip_r{0.0,0.0};

    // "centre of tilt" of the detectors on the z-axis
    // it will not coincide with the origin in case of non-standard detector tilt
    matrix<3,1> Oz_f, Oz_r;

    // rotation matrix to shift from one detector to the next consecutive anticlockwise one.
    // it follows the rule R^n_det = I
    matrix<3,3> R1_f, R1_r;

    // position vector of first detector face centre from centre of tilt
    matrix<3,1> p1_f, p1_r;

    // position vectors of all detector face centres
    // indices are recoil:[-n_rdet,-1] and forward:[1,n_fdet]
    map<int, matrix<3,1>> F;

    // "detector coordinate" unit vectors (r, theta, phi) of all detectors, arranged column-wise
    // indices are recoil:[-n_rdet,-1] and forward:[1,n_fdet]
    map<int, matrix<3,3>> U;
    
    explicit detector_properties(const std::string& path = "")
    {
        if (!path.empty()) read(path);

        sz_strip_f = { sz_f[0] / n_fstrip[0], sz_f[1] / n_fstrip[1] };
        sz_strip_r = { sz_r[0] / n_rstrip[0], sz_r[1] / n_rstrip[1] };

        // cout << "Forward detector strip size: " << sz_strip_f[0] << " m x " << sz_strip_f[1] << " m" << endl;
        // cout << "Recoil detector strip size: " << sz_strip_r[0] << " m x " << sz_strip_r[1] << " m" << endl;

        auto coft = [](array<double,3> pos, double tilt)
        {
            double rt = pow(pow(pos[0],2) + pow(pos[1],2),0.5);
            double Oz = pos[2] - rt / tan(tilt * rad);
            return Oz;
        };

        Oz_f[2] = coft(pos_f,tilt_f);
        Oz_r[2] = coft(pos_r,tilt_r);

        p1_f = {pos_f[0],pos_f[1],pos_f[2]};
        p1_r = {pos_r[0],pos_r[1],pos_r[2]};
        
        double phi_f = 2*M_PI / n_fdet;
        double phi_r = 2*M_PI / n_rdet;
        
        R1_f.sub<2,2>(0,0) = matrix<2,2>{cos(phi_f), -sin(phi_f), sin(phi_f), cos(phi_f)};
        R1_f[2][2] = 1;
        
        R1_r.sub<2,2>(0,0) = matrix<2,2>{cos(phi_r), -sin(phi_r), sin(phi_r), cos(phi_r)};
        R1_r[2][2] = 1;
        
        // minus sign multiplied to r coordinate keep the coordinate system right-handed (r, theta, phi)
        U[1](0) = Oz_f - p1_f; U[1](0) /= U[1](0).norm();
        U[1](2) = {-U[1](1,0), U[1](0,0), 0}; U[1](2) /= U[1](2).norm();
        U[1](1) = U[1](2) * U[1](0);

        // // minus sign multiplied to r coordinate keep the coordinate system right-handed (r, theta, phi)
        // U[-1][0] = (Oz_r - p1_r).tp(); U[-1][0] /= U[-1][0].norm();
        // U[-1][2] = {-U[-1](0,1), U[-1](0,0), 0}; U[-1][2] /= U[-1][2].norm();
        // U[-1][1] = (U[-1][2].tp() * U[-1][0].tp()).tp();

        // minus sign multiplied to r coordinate keep the coordinate system right-handed (r, theta, phi)
        U[-1](0) = Oz_r - p1_r; U[-1](0) /= U[-1](0).norm();
        U[-1](2) = {-U[-1](1,0), U[-1](0,0), 0}; U[-1](2) /= U[-1](2).norm();
        U[-1](1) = U[-1](2) * U[-1](0);

        for (int i = -static_cast<int>(n_rdet); i < static_cast<int>(n_fdet); i++)
        {
            int j = i<0 ? i : i+1;
            auto RR = (i<0 ? R1_r : R1_f) ^ (i<0 ? -i-1 : i);
            auto oz = (i<0 ? Oz_r : Oz_f);
        
            F[j] = RR * (i<0 ? p1_r : p1_f);
            // auto nn = F[j] - oz;

            // minus sign multiplied to keep the coordinate system right-handed (r, theta, phi)
            // nn = - nn / nn.norm();

            U[j] = RR * U[i<0 ? -1 : 1];// * RR.tp();
            // U[j] = U[j].tp();
            
            // double rxy = sqrt(pow(F[j][0],2) + pow(F[j][1],2));
            // matrix<3,1> uphi = {-F[j][1]/rxy, F[j][0]/rxy, 0};
            // U[j](0) = nn; U[j](1) = uphi*nn; U[j](2) = uphi;            // ^y = ^z cross ^x, ^z = ^phi, ^x = -^r, ^y = ^theta
            // cout<<"U["<<j<<"] = \n"<<U[j]<<endl;
            // cout<<"U[1/-1] = \n"<<U[j].divide_left(RR)<<endl;
        }
        
        // cout<<"p1_f = "<<(R1_f^n_fdet) * p1_f<<"p1_r = "<<(R1_r^n_rdet) * p1_r<<endl;
    }

    int regime(double, const matrix<3,1>&, const matrix<3,1>&);
    array<int,2> pixel(int, const matrix<3,1>&);
    matrix<3,1> pixel_to_loc(int, const array<int,2>&);
    hit read_hit(const event&);

    // template<size_t M>
    // static PolarimeterMetrics<M> scale_matrix(const matrix<M, 2>&, const matrix<M, 2>&, const matrix<M, 2>&, matrix<2, 1>&, matrix<2, 1>&);
    // template<size_t M>
    // static matrix<M,M> C2W_singular(const matrix<M,M>&);

    matrix<2,1,valarray<double>> theta_resolution(const target_properties&, const beam_properties&);
    matrix<N_R, 2> design_matrix(const target_properties&, const beam_properties&);
    matrix<N_R, N_R> cov_geom(const matrix<2, 1>&, const target_properties&, const beam_properties&);
protected:
    // List of keys required for a valid target_properties file
    std::vector<std::string> getRequiredKeys() const override
    {
        return {"sz_fdet", "pos_fdet", "tilt_fdet", "n_fdet", "n_fstrip",
                "sz_rdet", "pos_rdet", "tilt_rdet", "n_rdet", "n_rstrip"};
    }

    bool parseKey(const std::string& key, std::stringstream& ss) override
    {
        if (key == "sz_fdet")           ss >> sz_f[0] >> sz_f[1] >> sz_f[2];
        else if (key == "pos_fdet")     ss >> pos_f[0] >> pos_f[1] >> pos_f[2];
        else if (key == "tilt_fdet")    ss >> tilt_f;
        else if (key == "n_fdet")       ss >> n_fdet;
        else if (key == "n_fstrip")     ss >> n_fstrip[0] >> n_fstrip[1];
        else if (key == "sz_rdet")      ss >> sz_r[0] >> sz_r[1] >> sz_r[2];
        else if (key == "pos_rdet")     ss >> pos_r[0] >> pos_r[1] >> pos_r[2];
        else if (key == "tilt_rdet")    ss >> tilt_r;
        else if (key == "n_rdet")       ss >> n_rdet;
        else if (key == "n_rstrip")     ss >> n_rstrip[0] >> n_rstrip[1];
        else return false;

        return true;
    }

    void init();
};