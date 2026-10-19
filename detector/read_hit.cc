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
    
    auto it_hit = [&v, &r, this, &loc, &x](int reg)
    {
        auto n = U[reg](0), th = U[reg](1), ph = U[reg](2);
        auto &f = F[reg];
        auto sz = reg<0 ? sz_r : sz_f;

        if(r.tp()*n == 0) return false;
        
        double k = (f.tp()*n - v.tp()*n) / (r.tp()*n);
        // cout<<"k = "<<k<<endl;
        
        auto hpt = v + k*r - F[reg];

        bool a = abs(hpt.tp() * ph) < sz[0]/2;
        bool b = abs(hpt.tp() * th) < sz[1]/2;
        
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
    
    int reg = regime(E.charge(), v, r);
    
    // cout<<"detector id: "<<(i<0 ? i : i+1)<<endl;
    return it_hit(reg) ? hit(reg, loc, x, E.evID()) : hit();
}