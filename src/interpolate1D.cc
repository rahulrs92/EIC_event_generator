#include <vector>
#include <algorithm>
#include <stdexcept>

/**
 * Fast 1D Linear Interpolator
 * @param x_grid Vector of sorted x-coordinates (e.g., your CDF values from 0 to 1)
 * @param f_grid Vector of corresponding function values (e.g., your theta angles)
 * @param query_x The value you want to look up
 * @return The interpolated value
 */
double interpolate1D(const std::vector<double>& x_grid, 
                      const std::vector<double>& f_grid, 
                      double query_x) {
    
    // 1. Handle boundary conditions safely
    if (query_x <= x_grid.front()) return f_grid.front();
    if (query_x >= x_grid.back()) return f_grid.back();

    // 2. Binary search to find the upper bound of the interval
    // std::lower_bound finds the first element that is >= query_x
    auto it = std::lower_bound(x_grid.begin(), x_grid.end(), query_x);
    
    // 3. Get the indices of the enclosing interval [idx - 1, idx]
    size_t idx = std::distance(x_grid.begin(), it);
    size_t idx0 = idx - 1;
    size_t idx1 = idx;

    // 4. Perform linear interpolation math
    double x0 = x_grid[idx0];
    double x1 = x_grid[idx1];
    double f0 = f_grid[idx0];
    double f1 = f_grid[idx1];

    // Weight factor t (how far query_x is between x0 and x1, from 0 to 1)
    double t = (query_x - x0) / (x1 - x0);

    // Return the weighted average
    return f0 + t * (f1 - f0);
}