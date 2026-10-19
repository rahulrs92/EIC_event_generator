#include<detector_properties.hh>

hit detector_properties::read_hit(const event& E)
{
    auto v = E.vertex();
    // cout<<"v = "<<v<<endl;
    auto p = E.four_momentum();
    matrix<3,1> r = {p[1],p[2],p[3]};
    r = r / r.norm();
    // cout<<"r = "<<r<<endl;
    // auto m0 = pow(E.four_momentum().tp() * Lorentz::metric * E.four_momentum(),0.5);
    // cout<<"m0 = "<<m0<<endl;
    
    matrix<3,1> loc, x;
    
    auto it_hit = [&v, &r, &loc, &x](matrix<3,1> f, matrix<3,1> n, array<double,3> sz)
    {
        double k = (f.tp()*n - v.tp()*n) / (r.tp()*n);
        // cout<<"k = "<<k<<endl;
        
        double rxy = sqrt(pow(f[0],2)+pow(f[1],2));
        matrix<3,1> uphi = {-f[1]/rxy, f[0]/rxy};
        
        auto hpt = v + k*r - f;
        bool a = abs(hpt.tp() * uphi) < sz[0]/2;
        bool b = abs(hpt.tp() * (n*uphi)) < sz[1]/2;
        
        if(a && b)
        {
            // cout<<"v = "<<v<<"r = "<<r<<"f = "<<f<<endl;
            // cout<<"uphi = "<<uphi<<"(n*uphi) = "<<(n*uphi)<<endl;
            loc = v + k*r;
            x = sz[2] * r / (r.tp()*n);
            return true;
        }
        else return false;
    };
    
    // cout<<"n_rdet = "<<n_rdet<<", n_fdet = "<<n_fdet<<endl;
    for (int i = -static_cast<int>(n_rdet); i < static_cast<int>(n_fdet); i++)
    {
        auto RR = (i<0 ? R1_r : R1_f) ^ (i<0 ? -i-1 : i);
        
        auto ff = RR * (i<0 ? p1_r : p1_f);
        auto nn = ff - (i<0 ? Oz_r : Oz_f);
        nn = nn / nn.norm();

        auto &ssz = i<0 ? sz_r : sz_f;

        if(it_hit(ff, nn, ssz))
        {
            // cout<<"detector id: "<<(i<0 ? i : i+1)<<endl;
            return hit((i<0 ? i : i+1), loc, x, E.evID());
        }
    }
    return hit();
}