#include <iostream>
#include <chrono>

#include "cross_sections.hh"
#include "asymmetry.hh"
#include "sampler.hh"

#include "beam_properties.hh"
#include "target_properties.hh"
#include "detector_properties.hh"

#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLine.h"

// template PolarimeterMetrics<N_R> detector_properties::scale_matrix<N_R>(
//     const matrix<N_R, 2>&, 
//     const matrix<N_R, 2>&, 
//     const matrix<N_R, 2>&, 
//     matrix<2, 1>&, 
//     matrix<2, 1>&
// );

int main(int argc, char** argv)
{
    TApplication app("app", &argc, argv);

    beam_properties BP("input/beam.in");
    target_properties TP("input/target.in");
    detector_properties DP("input/detector.in");

    auto a = TP.get_fov(BP);
    cout << "Theoretical overlap fraction: " << a[4] << endl;
    double ep = TP.get_efficiency(BP);
    cout << "Theoretical polarimeter efficiency: " << ep << endl;
    auto thres = DP.theta_resolution(TP,BP);
    cout << "Theoretical polar angular resolution by detector (°):\n" << thres * deg << endl;

    double rthc = atan2(DP.sz_r[1]/2, DP.F[-1].norm());
    double rth = acos(DP.F[-1][2] / DP.F[-1].norm());
    size_t max_bins = 2 * rthc / thres[1].max();
    // cout << "strip size forward: " << DP.sz_strip_f[0] << ", " << DP.sz_strip_f[1] << endl;
    cout << "Maximum number of bins: " << max_bins << endl;
    cout << "Recoil theta range: " << (rth - rthc)*deg << ", " << (rth + rthc)*deg << endl;

    size_t nbins = max_bins/5;
    cout << "Number of bins set to: " << nbins << endl;

    TH3F *H3 = new TH3F("H3", "H3", 200, -0.2, 0.2, 200, -0.2, 0.2, 500, -0.1, 1);
    TH1F *H1p = new TH1F("H1p", "H1p", 513, -256, 256);
    TH1F *H1 = new TH1F("H1", "H1", nbins, (rth - rthc)*deg, (rth + rthc)*deg);
    TH2F *H2 = new TH2F("H2", "H2", 500, 6, 16, 500, 79, 87);
    TGraph *G1 = new TGraph();

    matrix<N_R,2,valarray<double>> rec_hits = valarray<double>(0.0,nbins);
    // map<array<int,3>, int> pixel_map;
    auto event_loop = [&](int i, int j)
    {
        auto evs = TP.generate_events(BP, i);
        auto h1 = DP.read_hit(evs[0]);
        auto h2 = DP.read_hit(evs[1]);
        if(h1)
        {
            auto l = h1.get_loc();
            int id = h1.get_ID();
            auto pixel = DP.pixel(id, l);
            matrix<3,1> vpix = DP.pixel_to_loc(id, pixel);
            H3->Fill(vpix[0], vpix[1], vpix[2]);
        }
        if(h2)
        {
            auto l = h2.get_loc();
            int id = h2.get_ID();
            auto pixel = DP.pixel(id, l);
            matrix<3,1> vpix = DP.pixel_to_loc(id, pixel);
            H3->Fill(vpix[0], vpix[1], vpix[2]);
            double th = acos(vpix[2]/vpix.norm())*180/M_PI;
            H1->Fill(th);
            int k = H1->FindBin(th);
            if(k>0 && k<=nbins) rec_hits(-id-1,j)[k-1] += 1.0;
            H1p->Fill(pixel[0]);
        }
        if(h1 && h2)
        {
            auto l1 = h1.get_loc(), l2 = h2.get_loc();
            int id1 = h1.get_ID(), id2 = h2.get_ID();
            auto pixel1 = DP.pixel(id1, l1), pixel2 = DP.pixel(id2, l2);
            matrix<3,1> vpix1 = DP.pixel_to_loc(id1, pixel1), vpix2 = DP.pixel_to_loc(id2, pixel2);
            H2->Fill(acos(vpix1[2]/vpix1.norm())*180/M_PI, acos(vpix2[2]/vpix2.norm())*180/M_PI);
        }
    };

    std::chrono::duration<double> ttk;
    auto t1 = std::chrono::high_resolution_clock::now();
    for(int i=0; i<100000000; i++) event_loop(i,0);
    BP = -BP;
    for(int i=0; i<100000000; i++) event_loop(i,1);
    BP = -BP;
    auto t2 = std::chrono::high_resolution_clock::now();
    ttk = t2 - t1;
    cout << "\nEvent Simulation Time: " << ttk.count() << " seconds\n" << endl;
    
    cout << "Total hit count per recoil detector:" << endl;
    cout << sum(rec_hits) << endl;

    auto get_theta_f = [&](double theta_r)
    {
        double m0_C = Bonin::m0_C;
        const double E1 = BP.kinetic_energy + m0_p;
        const double s = m0_p * m0_p + m0_C * m0_C + 2.0 * m0_C * E1;
        const double gamma = (E1 + m0_C) / std::sqrt(s);
        const double g = (m0_C * E1 + m0_p * m0_p) / (m0_C * (E1 + m0_C));
        
        const double A = 0.5 * (gamma * gamma * (1.0 + g) - (1.0 - g));
        const double B = 0.5 * (gamma * gamma * (1.0 + g) + (1.0 - g));

        // Returned lambda evaluates theta_f (in radians) for a given theta_r (in radians)
        const double sin_2tr = std::sin(2.0 * theta_r);
        const double cos_2tr = std::cos(2.0 * theta_r);
        return std::atan2(sin_2tr, A - B * cos_2tr);
    };
    
    TMultiGraph *mg = new TMultiGraph("Polarization by bin","Polarization by bin;bin number;P_{x,y}"); 
    TGraphErrors *g1u = new TGraphErrors(), *g1d = new TGraphErrors();
    TF1 *AP = new TF1("AP",Bonin::Ay,0,30,1);
    AP->SetParameter(0,BP.kinetic_energy);
    valarray<double> Px(nbins), Pxe(nbins), Py(nbins), Pye(nbins);
    for(int i=0; i<nbins; i++)
    {
        matrix<N_R,2> yld = pull(rec_hits,i), yld_e;
        for(int j=0; j<yld.M(); j++) for(int k=0; k<yld.N(); k++) yld_e(j,k) = sqrt(yld(j,k));
        double thr = H1->GetBinCenter(i+1) * rad;
        double thf = get_theta_f(thr) * deg;
        double Ay = AP->Eval(thf);
        auto R = Ay * DP.design_matrix(TP,BP);
        // cout<<"Design matrix:\n"<<R<<endl;
        
        matrix<2, 1> Ab;   // Output physical asymmetries: Column 0 = A_bs, Column 1 = A_bc
        matrix<2, 1> Ab_e; // Output physical asymmetry uncertainties: Column 0 = sigma(A_bs), Column 1 = sigma(A_bc)
        auto metrics = Polarimeter::scale_matrix(yld, yld_e, R, Ab, Ab_e);
        auto cg = DP.cov_geom(Ab,TP,BP);
        // cout<<"Covariances:\n"<<cat(metrics.cov_A,cg)<<endl;
        auto C = metrics.cov_A + cg;    // Complete Covariance matrix for asymmetries
        
        double sf = C.norm() / C.M(); // Scale factor for numerical stability
        auto W = (C / sf).pseudo_inverse(1e-10); // Use pseudo-inverse with threshold for numerical stability
        W /= sf; // Rescale back to original units
        
        auto Cb = (R.tp() * W * R) ^ -1;
        Ab = Cb * R.tp() * W * metrics.A;
        Px[i] = Ab[0]; Py[i] = Ab[1];
        Ab_e = {sqrt(Cb(0,0)),sqrt(Cb(1,1))};
        Pxe[i] = Ab_e[0]; Pye[i] = Ab_e[1];
        
        g1u->AddPointError(i,Ab[0],0,Ab_e[0]); g1d->AddPointError(i,Ab[1],0,Ab_e[1]);
    }
    double PX = (Px * pow(Pxe,-2)).sum() / pow(Pxe,-2).sum();
    double PXE = 1.0 / sqrt(pow(Pxe,-2).sum());
    double PY = (Py * pow(Pye,-2)).sum() / pow(Pye,-2).sum();
    double PYE = 1.0 / sqrt(pow(Pye,-2).sum());

    cout<<"Polarization: ("<<PX<<", "<<PY<<") +/- ("<<PXE<<", "<<PYE<<")"<<endl;
    
    gStyle->SetNumberContours(99);  // More color levels = smoother gradient
    
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    g1u->SetTitle("P_{x}"); g1d->SetTitle("P_{y}");
    g1u->SetMarkerColor(4); g1u->SetMarkerStyle(21);
    g1d->SetMarkerColor(2); g1d->SetMarkerStyle(22);
    // g1u->SetDrawOption("P");g1d->SetDrawOption("P");
    mg->Add(g1u); mg->Add(g1d);
    mg->Draw("AP");
    c2->SetGrid();
    auto leg = c2->BuildLegend(0.7,0.9,0.9,0.8);
    leg->SetFillColorAlpha(10,0.2);
    leg->Draw();

    TLine* pxl = new TLine(0,BP.polarization[0],nbins-1,BP.polarization[0]);
    pxl->SetLineColor(kBlack);
    pxl->SetLineStyle(2);
    pxl->Draw("f");

    TLine* pyl = new TLine(0,BP.polarization[1],nbins-1,BP.polarization[1]);
    pxl->SetLineColor(kBlack);
    pxl->SetLineStyle(2);
    pxl->Draw("f");

    TBox* pxb = new TBox(0, PX-PXE, nbins-1, PX+PXE);
    pxb->SetFillColorAlpha(kBlue, 0.1);
    pxb->Draw();

    TBox* pyb = new TBox(0, PY-PYE, nbins-1, PY+PYE);
    pyb->SetFillColorAlpha(kRed, 0.1);
    pyb->Draw();
    
    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    H2->Draw("colz");
    H2->SetTitle("Distribution of coincident elastic p-C scattered polar angles (MC);#theta_{f} [deg];#theta_{r} [deg]; frequency");
    c3->SetGrid();
    c3->SaveAs("coincidence.pdf", "pdf");

    TCanvas *c4 = new TCanvas("c4", "c4", 400, 1200);
    // c4->SetLogz();
    H3->Draw("scat");
    H3->SetStats(0);
    c4->SetTheta(5.0);
    c4->SetPhi(20.0);
    // c4->SetGrid();
    H3->SetTitle("Distribution for elastic p-C scattered hits (MC);x [m];y [m];z [m]");
    c4->SaveAs("hits.pdf", "pdf");

    TCanvas *c5 = new TCanvas("c5", "c5", 800, 600);
    // c5->SetLogz();
    H1->Draw();
    c5->SetGrid();
    H1->SetTitle("Distribution of Recoil Angles from Pixel Readout;theta [deg];counts");
    c5->SaveAs("recoil_angles.pdf", "pdf");

    TCanvas *c6 = new TCanvas("c6", "c6", 800, 600);
    // c6->SetLogz();
    H1p->Draw();
    c6->SetGrid();
    H1p->SetTitle("Distribution of Raw Recoil Pixel hits;pixel;counts");
    c6->SaveAs("recoil_pixels.pdf", "pdf");

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
    
    return 0;
}

int main_asymmetry(int argc, char** argv)
{
    size_t constexpr M = 6;     // Number of detectors

    // // Yield counts matrix: Column 0 = y_{r+}, Column 1 = y_{r-}
    // matrix<M, 2> yields = {
    // 96059.0, 102215.0,  // Chip 1
    // 94532.0,  99893.0,  // Chip 2
    // 163075.0, 146204.0,  // Chip 3
    // 153642.0, 132488.0,  // Chip 4
    // 161036.0, 107847.0,  // Chip 5
    // 165097.0, 107054.0,  // Chip 6
    // 137717.0,  81521.0,  // Chip 7
    // 147319.0,  87814.0,  // Chip 8
    // 158090.0, 111069.0,  // Chip 9
    // 152154.0, 111142.0,  // Chip 10
    // 119777.0, 112654.0,  // Chip 11
    // 120060.0, 116626.0   // Chip 12
    // }, yields_err;
    
    // // Geometric weights matrix: Column 0 = c_r (<Ay cos phi>_r), Column 1 = s_r (<Ay sin phi>_r)
    // matrix<M, 2> cs = {
    // +0.239782, -0.017039,  // Chip 1
    // +0.239782, +0.017039,  // Chip 2
    // +0.134647, +0.199138,  // Chip 3
    // +0.105135, +0.216177,  // Chip 4
    // -0.105135, +0.216177,  // Chip 5
    // -0.134647, +0.199138,  // Chip 6
    // -0.239782, +0.017039,  // Chip 7
    // -0.239782, -0.017039,  // Chip 8
    // -0.134647, -0.199138,  // Chip 9
    // -0.105135, -0.216177,  // Chip 10
    // +0.105135, -0.216177,  // Chip 11
    // +0.134647, -0.199138   // Chip 12
    // };

    // matrix<M, 2> yields = {
    //     4., 1.,
    //     3., 3.,
    //     6., 1.,
    //     4., 5.,
    //     2., 7.,
    //     7., 3.
    // }, yields_err;
    
    // matrix<M, 2> cs;
    // for(size_t i = 0; i < M; ++i)
    // {
    //     double phi = 2 * M_PI * i / M;
    //     cs[i][0] = 0.1070 * cos(phi);
    //     cs[i][1] = 0.1070 * sin(phi);
    // }
    
    // Raw Yield Matrix: Col 0 = y_{r+}, Col 1 = y_{r-}
    matrix<6, 2> yields = {
        85680.0,  92000.0,  // Det S1 (0 deg)
        143134.0, 126402.0,  // Det S2 (60 deg)
        151996.0,  99368.0,  // Det S3 (120 deg)
        130410.0,  76500.0,  // Det S4 (180 deg)
        139358.0,  99398.0,  // Det S5 (240 deg)
        108131.0, 104182.0   // Det S6 (300 deg)
    }, yields_err;

    // Geometric Weights Matrix: Col 0 = c_r (<Ay cos phi>), Col 1 = s_r (<Ay sin phi>)
    matrix<6, 2> cs = {
        +0.250000, +0.000000,  // Det S1 (0 deg)
        +0.125000, +0.216506,  // Det S2 (60 deg)
        -0.125000, +0.216506,  // Det S3 (120 deg)
        -0.250000, +0.000000,  // Det S4 (180 deg)
        -0.125000, -0.216506,  // Det S5 (240 deg)
        +0.125000, -0.216506   // Det S6 (300 deg)
    };
    cs = cs * matrix<2,2>{0,1,1,0}; // Swap columns to match the expected order: Col 0 = s_r, Col 1 = c_r

    for(size_t i = 0; i < M; ++i) for(size_t j = 0; j < 2; ++j) yields_err[i][j] = sqrt(yields[i][j]);

    matrix<2, 1> Ab;   // Output physical asymmetries: Column 0 = A_bs, Column 1 = A_bc
    matrix<2, 1> Ab_e; // Output physical asymmetry uncertainties: Column 0 = sigma(A_bs), Column 1 = sigma(A_bc)

    // detector_properties DP("input/detector.in");
    // Geometric Weight Matrix R: +ve for forward detectors, -ve for recoil detectors 
    auto metrics = Polarimeter::scale_matrix(yields, yields_err, cs, Ab, Ab_e);

    cout << "Scaled physical asymmetries:" << endl;
    cout << "P_x = " << Ab[0] << " ± " << Ab_e[0] << endl;
    cout << "P_y = " << -Ab[1] << " ± " << Ab_e[1] << endl;
    cout << "Relative beam flux ratio: L = " << metrics.L << " ± " << metrics.L_e << endl;
    cout << "Detector asymmetries (A_i):\n" << metrics.A << endl;
    // cout << "Covariance matrix of physical asymmetries:" << endl;
    // cout << metrics.cov_Ab << endl;

    matrix<M,2> X = cat(1.0 + metrics.A, 1.0 - metrics.A), tr = 1.0/M, ep;
    double s[2] = { 0.0, 0.0 };
    for(size_t i = 0; i < M; ++i) for(size_t j = 0; j < 2; ++j) {ep[i][j] = yields[i][j] / X[i][j]; s[j] += ep[i][j];}
    ep(0) /= s[0] / M; ep(1) /= s[1] / M;
    cout << "Normalized detector efficiencies (ep_i):\n" << ep << endl;
    return 0;
}

