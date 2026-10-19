#pragma once

#include<sampler.hh>

#include<cross_sections.hh>
#include<beam_properties.hh>
#include<event.hh>

class target_properties : public config_properties
{
public:
    // Configuration properties
    double d_eff{0.0};                                          // in m
    double d_t{0.0};                                            // in g/cm^2
    std::array<double, 2> position{0.0, 0.0};                   // (tpx, tpy)
    double l_t{0.0};                                            // in m
    double theta{0.0};                                          // in degrees

    // Default initializations
    sampler<2,5> pdf;
    sampler<2,0> gaus2;

    // Calculated properties
    std::array<double, 5> fov{0.0, 0.0, 0.0, 0.0, 0.0};         // (x_min, x_max, y_min, y_max, fov) in m
    matrix<2,2> S, R;
    
    // flags
    bool f_fov = false;                                         // Flag to indicate if fov has been calculated

    explicit target_properties(const std::string& path = ""): 
    pdf("target_CS",Bonin::f_theta_phi), 
    gaus2("target_Gaus2",Distributions::gaus2_standard) 
    {
        if (!path.empty()) read(path);
        // double x[2] = {1,0};
        // double p[5] = {1.2, d_eff, d_t, 1, 0};
        // std::cout<<"function value at min = "<<pdf(x,p)<<std::endl;
    }

    std::array<double, 5> get_fov(const beam_properties&);
    const std::array<double, 5> get_fov() const { return fov; }

    std::array<event,2> generate_events(const beam_properties&, long);

    double get_efficiency(beam_properties&);

protected:
    // List of keys required for a valid target_properties file
    std::vector<std::string> getRequiredKeys() const override
    {
        return {"d_eff", "d_t", "position", "l_t", "theta"};
    }

    bool parseKey(const std::string& key, std::stringstream& ss) override
    {
        if (key == "d_eff")             ss >> d_eff;
        else if (key == "d_t")          ss >> d_t;
        else if (key == "position")     ss >> position[0] >> position[1];
        else if (key == "l_t")          ss >> l_t;
        else if (key == "theta")        ss >> theta;
        else return false;

        return true;
    }
};