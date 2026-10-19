#include<detector_properties.hh>

int detector_properties::regime(double q, const matrix<3,1> &v, const matrix<3,1> &r)
{
    auto &f = (q==1 ? p1_f : p1_r);
    auto &nd = (q==1 ? n_fdet : n_rdet);

    matrix<3,1> o = {0, 0, f[2]};
    double k = (f[2] - v[2]) / r[2];
    matrix<3,1> a = v + k*r - o;
    
    double phi = atan2(a[1],a[0]) - atan2(f[1],f[0]);
    phi += (phi<0 ? 2*M_PI : 0);
    
    double arc = phi*nd/(2*M_PI) + 0.5;
    int reg = int(arc)%nd + 1;
    
    return (q==1 ? 1 : -1) * reg;
}