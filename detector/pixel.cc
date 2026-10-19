#include<detector_properties.hh>
#include<hit.hh>

array<int,2> detector_properties::pixel(int ID, const matrix<3,1>& loc)
{
    auto &sz = ID>0 ? sz_strip_f : sz_strip_r;
    if(ID==0) return {0,0};

    int pix_theta, pix_phi;
    auto pix = U[ID].tp() * (loc - F[ID]);

    if((ID>0 ? n_fstrip[1] : n_rstrip[1]) % 2 == 0) pix_theta = pix[1] / sz[1] + (pix[1]<0 ? -1 : 1);
    else pix_theta = pix[1] / sz[1] + (pix[1]<0 ? -0.5 : 0.5);

    if((ID>0 ? n_fstrip[0] : n_rstrip[0]) % 2 == 0) pix_phi = pix[2] / sz[0] + (pix[2]<0 ? -1 : 1);
    else pix_phi = pix[2] / sz[0] + (pix[2]<0 ? -0.5 : 0.5);
    // if(ID>0) cout << "pix_phi (forward): " << n_fstrip[0] << endl;

    return {pix_theta, pix_phi};
}

matrix<3,1> detector_properties::pixel_to_loc(int ID, const array<int,2>& pixel)
{
    auto &sz = ID>0 ? sz_strip_f : sz_strip_r;
    if(ID==0) return {0,0};

    matrix<3,1> loc = {0,0,0};

    if((ID>0 ? n_fstrip[1] : n_rstrip[1]) % 2 == 0) loc[1] = (pixel[0]<0 ? pixel[0] + 0.5 : pixel[0] - 0.5) * sz[1];
    else loc[1] = pixel[0] * sz[1];

    if((ID>0 ? n_fstrip[0] : n_rstrip[0]) % 2 == 0) loc[2] = (pixel[1]<0 ? pixel[1] + 0.5 : pixel[1] - 0.5) * sz[0];
    else loc[2] = pixel[1] * sz[0];

    loc = F[ID] + U[ID] * loc;
    return loc;
}