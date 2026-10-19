#pragma once

#include<config_properties.hh>

class beam_properties : public config_properties
{
public:
    // Default initializations
    std::string particle_id;
    double kinetic_energy{0.0};                     // in GeV
    std::array<double, 2> polarization{0.0, 0.0};   // (Px, Py)
    double emittance{0.0};                         // in m.rad
    std::array<double, 2> beta_function{0.0, 0.0}; // (beta_x, beta_y) in m
    uint64_t n_cyc{0};                              // Particles per bunch

    explicit beam_properties(const std::string& path = "") { if (!path.empty()) read(path); get_sigma(); }

    // Calculated Properties
    std::array<double,2> sigma{0.0,0.0};

    std::array<double, 2> spin_sample();

    std::array<double,2> get_sigma();
    const std::array<double,2> get_sigma() const { return sigma; };

    // Spin Polarization Flipping
    beam_properties& operator- () { for(auto& a : polarization) a = -a; return *this; }
    
protected:
    
    // List of keys required for a valid beam_properties file
    std::vector<std::string> getRequiredKeys() const override
    {
        return {"pID", "kinetic_energy", "polarization", "emittance", "beta_function", "N_cyc"};
    }

    bool parseKey(const std::string& key, std::stringstream& ss) override
    {
        if (key == "pID")                  ss >> particle_id;
        else if (key == "kinetic_energy")  ss >> kinetic_energy;
        else if (key == "polarization")    ss >> polarization[0] >> polarization[1];
        else if (key == "emittance")       ss >> emittance;
        else if (key == "beta_function")   ss >> beta_function[0] >> beta_function[1];
        else if (key == "N_cyc")           ss >> n_cyc;
        else return false;

        return true;
    }
};