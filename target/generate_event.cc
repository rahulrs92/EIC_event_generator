#include<target_properties.hh>

#include<sampler.hh>

array<event,2> target_properties::generate_events(const beam_properties &BP, long eID)
{
    double x_min[2] = {0.0, -M_PI};
    double x_max[2] = {30.0, M_PI};

    pdf.set_parameters(BP.kinetic_energy, d_eff, d_t, BP.polarization[0], BP.polarization[1]);
    // cout<<"Kinetic Energy set to "<<pdf.get_parameters()[0]<<" GeV"<<endl;
    // cout<<"Has parameters: "<<pdf.hasParameters()<<endl;
    pdf.set_range(x_min[0], x_max[0], x_min[1], x_max[1]);
    // cout<<"Has range: "<<pdf.hasRange()<<endl;
    // cout<<"yay!!"<<endl;

    double E1 = m0_p + BP.kinetic_energy;
    double E2 = Bonin::m0_C;
    matrix<4,1> p1 = {E1, 0, 0, sqrt(pow(E1,2) - pow(m0_p,2))};
    matrix<4,1> p2 = {E2, 0, 0, 0};

    auto s1 = pdf.random_sample();
    // cout<<"s1[0]: "<<s1[0]<<", s1[1]: "<<s1[1]<<endl;
    double a = pow(E1 + E2, 2) - pow(p1[3], 2)*pow(cos(s1[0] * rad), 2);
    double b = -2 * (E1 + E2) * (pow(m0_p,2) + E1*E2);
    double c = pow(pow(m0_p,2) + E1*E2,2) + pow(m0_p,2)*pow(p1[3],2)*pow(cos(s1[0] * rad),2);

    if(b*b - 4*a*c < 0)
    {
        cout<<"Discriminant is not positive!!! ("<<b*b - 4*a*c<<")"<<endl;
        cout<<"Resampling..."<<endl;
        return generate_events(BP, eID);
    }
    double E3 = (-b + sqrt(b*b - 4*a*c)) / (2*a);
    double p3_m = sqrt(pow(E3,2) - pow(m0_p,2));
    matrix<4,1> p3 = {E3, p3_m*sin(s1[0]*rad)*cos(s1[1]), p3_m*sin(s1[0]*rad)*sin(s1[1]), p3_m*cos(s1[0]*rad)};
    matrix<4,1> p4 = p1 + p2 - p3;

    auto ov = get_fov();
    // cout << "Theoretical overlap fraction: " << a[4] << endl;
    gaus2.set_range(ov[0], ov[1], ov[2], ov[3]);
    auto s2 = gaus2.random_sample();
    auto ts2 = S * R * matrix<2,1>{s2[0],s2[1]};
    matrix<3,1> v = {ts2[0], ts2[1], 0};

    return {event(eID,1,v,p3), event(eID,12,v,p4)};
}