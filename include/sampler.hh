#pragma once

#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <optional>
#include <functional>

#include "TF1.h"
#include "TF2.h"
#include "TRandom3.h"

template <size_t Nx, size_t Np>
class sampler
{
    using func_t = std::function<double(double*, double*)>;
    using result_t = std::conditional_t<Nx == 1, double, std::array<double, Nx>>;
    using callable_t = std::conditional_t<Nx == 1, TF1*, std::conditional_t<Nx == 2, TF2*, func_t>>;

private:
    const func_t f;
    std::optional<std::array<double, Np>> par; // Store parameters for the function
    std::optional<std::array<double, Nx>> x_min, x_max; // Store minimum and maximum range for each dimension

    // TRandom3* rand_gen; // Random number generator
    TF1* pdf1D; // For 1D random sampling
    TF2* pdf2D; // For 2D random sampling

    bool f_r1 = true; // Flag to indicate if it is the first random sample

    void set_parameters()
    {
        if constexpr (Nx == 1)
        {
            pdf1D->SetParameters(par->data());
        }
        else if constexpr (Nx == 2)
        {
            pdf2D->SetParameters(par->data());
        }
    }

    void set_range()
    {
        if constexpr (Nx == 1)
        {
            pdf1D->SetRange((*x_min)[0], (*x_max)[0]);
            pdf1D->SetNpx(1000); // Set number of points for graphical representation
        }
        else if constexpr (Nx == 2)
        {
            pdf2D->SetRange((*x_min)[0], (*x_min)[1], (*x_max)[0], (*x_max)[1]);
            pdf2D->SetNpx(500); // Set number of points along x for graphical representation
            pdf2D->SetNpy(500); // Set number of points along y for graphical representation
        }
    }

    
    public:
    // State checking
    bool hasRange() const { return x_min.has_value() && x_max.has_value(); }
    bool hasParameters() const { return Np==0 ? true : par.has_value(); }
    bool isReady() const { return hasRange() && hasParameters(); }

    // Construct from an existing callable. We store it as a const reference,
    // so the caller must ensure the lifetime of 'f' outlives this object.
    explicit sampler(std::string iname, const func_t &iff) : f(iff)
    {
        // rand_gen = new TRandom3(12345); // Seed the random number generator for reproducibility
        if constexpr (Nx == 1) pdf1D = new TF1(("pdf_" + iname).data(), f, 0, 1, Np);
        if constexpr (Nx == 2) pdf2D = new TF2(("pdf_" + iname).data(), f, 0, 1, 0, 1, Np);
    }
    explicit sampler(std::string iname, const func_t &iff, const double *params) : sampler(iname, iff)
    {
        set_parameters(params);
    }
    explicit sampler(std::string iname, const func_t &iff, double* ix_min, double* ix_max) : sampler(iname, iff)
    {
        set_range(ix_min,ix_max);
    }
    explicit sampler(std::string iname, const func_t &iff, const double *params, double *ix_min, double *ix_max) : sampler(iname, iff)
    {
        set_parameters(params);
        set_range(ix_min,ix_max);
    }
    // Copy constructor, defaulted (reference member will be copied as a reference)
    sampler(const sampler &) = default;
    // Deleted default constructor to enforce providing a callable
    sampler() = delete;


    callable_t operator -> () const
    {
        if constexpr (Nx == 1) return pdf1D;
        else if constexpr (Nx == 2) return pdf2D;
        else return f;
    }
    

    void set_parameters(double* params)
    {
        // Compare elements if already set
        if (par.has_value() && std::equal(par->begin(), par->end(), params)) return;
        // {
        //     std::cout << "[Params] New parameters match existing parameters. Skipping assignment.\n";
        //     return;
        // }
        par.emplace();
        for(int i=0; i<Np; i++) (*par)[i] = params[i];
        set_parameters();
    }

    template <typename... Args, typename = std::enable_if_t
                  <(sizeof...(Args) == Np) && (std::is_convertible_v<Args, double> && ...)>>
    void set_parameters(Args... args)
    {
        std::array<double, Np> values{static_cast<double>(args)...};
        set_parameters(values.data());
        // for (size_t i = 0; i < Np; ++i) par[i] = values[i];
        // f_lpar = true;
        // if(f_rng) set_parameters();
    }

    const double* get_parameters() const { return (Nx==1 ? pdf1D : pdf2D)->GetParameters(); }



    void set_range(double *ix_min, double *ix_max)
    {
        bool minEqual = x_min.has_value() && std::equal(x_min->begin(), x_min->end(), ix_min);
        bool maxEqual = x_max.has_value() && std::equal(x_max->begin(), x_max->end(), ix_max);

        // Skip assignment if both min and max match existing values
        if (minEqual && maxEqual) return;
        // {
        //     std::cout << "[Range] New x_min and x_max match existing range. Skipping assignment.\n";
        //     return;
        // }
        x_min.emplace();
        for(int i=0; i<Nx; i++) (*x_min)[i] = ix_min[i];
        x_max.emplace();
        for(int i=0; i<Nx; i++) (*x_max)[i] = ix_max[i];
        set_range();
    }

    template <typename... Args, typename = std::enable_if_t
                <(sizeof...(Args) == 2*Nx) && (std::is_convertible_v<Args, double> && ...)>>
    void set_range(Args... args)
    {
        std::array<double, 2*Nx> values{static_cast<double>(args)...};
        double ix_min[Nx], ix_max[Nx];
        for (size_t i = 0; i < Nx; ++i)
        {
            ix_min[i] = values[2*i];
            ix_max[i] = values[2*i + 1];
        }
        set_range(ix_min,ix_max);
    }

    const double* get_range_mins() const { return x_min.has_value() ? x_min->data() : nullptr; }
    const double* get_range_maxs() const { return x_max.has_value() ? x_max->data() : nullptr; }
    



    // Invoke the stored callable
    double operator()(double *x, double *par) const { return f(x, par); }

    template <typename... Args, 
              typename = std::enable_if_t<
                  (sizeof...(Args) == Nx || sizeof...(Args) == Nx + Np) && 
                  (std::is_convertible_v<Args, double> && ...)
              >>
    double operator()(Args... args) const
    {
        if(sizeof...(Args) == Nx)
        {
            std::array<double, Nx> values{static_cast<double>(args)...};
            double x[Nx], p_t[Np];
            for (size_t i = 0; i < Nx; ++i) x[i] = values[i];
            if(hasParameters() || Np == 0) // If parameters are set or there are no parameters
            {
                std::copy(par->begin(), par->end(), p_t);
                return f(x, p_t);
            }
            else
            {
                std::cerr << "warning: Parameters not set. Assuming zeroes.\n";
                std::vector<double> zero_par(Np, 0);
                return f(x, zero_par.data());
            }
        }
        else if(sizeof...(Args) == Nx + Np)
        {
            std::array<double, Nx + Np> values{static_cast<double>(args)...};
            double x[Nx], p[Np];
            for (size_t i = 0; i < Nx; ++i) x[i] = values[i];
            for (size_t i = 0; i < Np; ++i) p[i] = values[Nx + i];
            return f(x, p);
        }
    }

    
    result_t random_sample()
    {
        result_t x;
        if (!isReady())
        {
            throw std::runtime_error("Range and parameters must be set before calling random_sample().");
        }
        if constexpr (Nx == 1)
        {
            x = pdf1D->GetRandom();
        }
        else if constexpr (Nx == 2)
        {
            pdf2D->GetRandom2(x[0], x[1]);
            // if(f_r1) // If it's the first random sample, print the sample for debugging
            // {
            //     std::cout << "First random sample: " << x[0] << ", " << x[1] << "\n";
            //     std::cout << "Parameters: ";
            //     for (size_t i = 0; i < Np; ++i) std::cout << pdf2D->GetParameter(i) << " ";
            //     std::cout << "\n";
            //     std::cout << "Range: x = [" << pdf2D->GetXmin() << ", " << pdf2D->GetXmax() << "], "
            //                 << "y = [" << pdf2D->GetYmin() << ", " << pdf2D->GetYmax() << "]\n";
            // }
            // std::cout << "Random sample: " << x[0] << ", " << x[1] << "\n"; // Debugging output
        }
        else
        {
            throw std::runtime_error("random_sample() is for now only implemented for Nx=1 or Nx=2.");
        }
        if(f_r1) f_r1 = false; // Mark that the first random sample has been generated
        return x;
    }

    result_t random_sample(double *ix_min, double *ix_max)
    {
        set_range(ix_min,ix_max);
        return random_sample();
    }


};
