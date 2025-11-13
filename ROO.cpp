#include <iostream>
#include <fstream>
#include <cmath>
#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TString.h"
#include "/alf/data/mleibelt/LogPlotLib/Plotting.h"




double pT_bins[13] = {1.8, 2.2, 2.6, 3.2, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 16.0, 20.0};
std::vector<double> ptCenters;
std::vector<double> rooValues;

void ROO(const TString& datarec_OO_file, 
        const TString& datarec_pp_file,
        const TString& eff_file_OO,
        const TString& eff_file_pp, 
        int Nevt_OO, int Nevt_pp,
        const std::vector<int>& goodbins,
        double linethickness){
    
    
    std::ofstream logfile("logfiles/logfile_ROO.txt");
    std::streambuf* originalCoutBuffer = std::cout.rdbuf();
    std::cout.rdbuf(logfile.rdbuf());

    std::cout << "Opening input files..." << std::endl;

    //input parameter
    double deltay = 0.8; 
    double eff_TVX_pp = 0.72;
    double eff_TVX_OO = 0.935;           
    double sig_TVX_pp = 0.0503; //b
    double sig_TVX_OO = 1.13;   //b

    //get histograms
    TFile* frec_OO = TFile::Open(datarec_OO_file, "READ");
    TFile* frec_pp = TFile::Open(datarec_pp_file, "READ");
    TFile *input_eff_OO = TFile::Open(eff_file_OO, "READ");
    TFile *input_eff_pp = TFile::Open(eff_file_pp, "READ");
    
    if (!frec_OO || !frec_pp) {
        std::cerr << "Fehler beim Öffnen der Dateien!" << std::endl;
        std::cout.rdbuf(originalCoutBuffer);
        return;
    }
    
    // Create histogram with variable binning
    TH1D *hROO = new TH1D("hROO", "R_{OO} vs p_{T};p_{T} (GeV/c);R_{OO}", 
                      12, pT_bins);  // 12 bins, edges from pT_bins

    for (int bin : goodbins) {

        TH1* hDataRec_OO = dynamic_cast<TH1*>(frec_OO->Get("h_omegacounts"));
        TH1* hDataRec_pp = dynamic_cast<TH1*>(frec_pp->Get("h_omegacounts"));
        TH1* heff_OO = dynamic_cast<TH1*>(input_eff_OO->Get("hEfficiency"));
        TH1* heff_pp = dynamic_cast<TH1*>(input_eff_pp->Get("hEfficiency"));
        

        double ptMin = pT_bins[bin];
        std::cout<<"\npT_bins[bin]: "<<pT_bins[bin]<<std::endl;
        double ptMax = pT_bins[bin+1];
        double ptCenter = 0.5 * (ptMin + ptMax);
        double deltapT = ptMax-ptMin;

        int binIndex_OO = hDataRec_OO->GetXaxis()->FindBin(ptCenter);
        int binIndex_pp = hDataRec_pp->GetXaxis()->FindBin(ptCenter);
        double nomega_OO = hDataRec_OO->GetBinContent(binIndex_OO);
        double nomega_pp = hDataRec_pp->GetBinContent(binIndex_pp);

        int binIndex_effOO = heff_OO->GetXaxis()->FindBin(ptCenter);
        int binIndex_effPP = heff_pp->GetXaxis()->FindBin(ptCenter);
        double eff_OO = heff_OO->GetBinContent(binIndex_effOO);
        double eff_pp = heff_pp->GetBinContent(binIndex_effPP);

        
        double err_nomega_OO = hDataRec_OO->GetBinError(binIndex_OO);
        double err_nomega_pp = hDataRec_pp->GetBinError(binIndex_pp);
        double err_eff_OO    = heff_OO->GetBinError(binIndex_effOO);
        double err_eff_pp    = heff_pp->GetBinError(binIndex_effPP);

        std::cout<<"binIndex_OO: "<<binIndex_OO<<", binIndex_pp: "<<binIndex_pp<<std::endl;
        std::cout<<"binIndex_effOO: "<<binIndex_effOO<<", binIndex_effPP: "<<binIndex_effPP<<std::endl;
        

        double yield_OO = nomega_OO/Nevt_OO * 1/(eff_OO*2*M_PI*deltay) * 1/(ptCenter*deltapT) * eff_TVX_OO;
        double yield_pp = nomega_pp/Nevt_pp * 1/(eff_pp*2*M_PI*deltay) * 1/(ptCenter*deltapT) * eff_TVX_pp; 
        
        double crosssection_OO = nomega_OO/Nevt_OO * 1/(eff_OO*2*M_PI*deltay) * 1/(ptCenter*deltapT) * sig_TVX_OO;
        double crosssection_pp = nomega_pp/Nevt_pp * 1/(eff_pp*2*M_PI*deltay) * 1/(ptCenter*deltapT) * sig_TVX_pp;
        int nuclear_thickness_OO = 188;    //b^-1

        double ROO = yield_OO / (nuclear_thickness_OO*crosssection_pp);
        std::cout<<"ptCenter: "<<ptCenter<<", deltapT: "<<deltapT<<std::endl;
        std::cout<<"nomega_OO: "<<nomega_OO<<", eff_OO: "<<eff_OO<<", yield_OO: "<<yield_OO<<", crosssection_OO: "<<crosssection_OO<<std::endl;
        std::cout<<"nomega_pp: "<<nomega_pp<<", eff_pp: "<<eff_pp<<", yield_pp: "<<yield_pp<<", crosssection_pp: "<<crosssection_pp<<std::endl;
        std::cout<<"ROO: "<<ROO<<std::endl;

        ptCenters.push_back(ptCenter);
        rooValues.push_back(ROO);

        double relErr_OO = (nomega_OO > 0 && eff_OO > 0) ?
            sqrt( pow(err_nomega_OO/nomega_OO,2) + pow(err_eff_OO/eff_OO,2) ) : 0;
        double relErr_pp = (nomega_pp > 0 && eff_pp > 0) ?
            sqrt( pow(err_nomega_pp/nomega_pp,2) + pow(err_eff_pp/eff_pp,2) ) : 0;
        double errROO = ROO * sqrt( pow(relErr_OO,2) + pow(relErr_pp,2) );
        std::cout<<"errROO: "<<errROO<<std::endl;
        int histBin = hROO->FindBin(ptCenter);
        hROO->SetBinContent(histBin, ROO);
        hROO->SetBinError(histBin, errROO);
    }

    Plotting1D Plot_ROO;

    // --- RpPb 5 TeV omega ---
    TFile *RpPb_nico = TFile::Open("Nico_results/Final_Omega_5TeV_Results.root","read");
    if (!RpPb_nico || RpPb_nico->IsZombie()) {
        std::cerr << "ERROR: could not open Final_Omega_5TeV_Results.root" << std::endl;
    } else {
        TDirectory* dir_RpPb = (TDirectory*)RpPb_nico->Get("RpPb");
        if (!dir_RpPb) {
            std::cerr << "ERROR: directory 'RpPb' not found!" << std::endl;
        } else {
            TGraphAsymmErrors* g_RpPb_stat = dynamic_cast<TGraphAsymmErrors*>(dir_RpPb->Get("RpAStat"));
            TGraphAsymmErrors* g_RpPb_sys  = dynamic_cast<TGraphAsymmErrors*>(dir_RpPb->Get("RpASys"));
            if (!g_RpPb_stat || !g_RpPb_sys) {
                std::cerr << "ERROR: graphs 'RpAStat' or 'RpASys' not found!" << std::endl;
            } else {
                Plot_ROO.New(g_RpPb_stat,"#omega p-Pb, #sqrt{s_{NN}} = 5.02 TeV - ALICE Preliminary", 20, linethickness-2, kGreen-8, "ep");
                Plot_ROO.New(g_RpPb_sys,"", 20, linethickness-2, kGreen-8, "Box");
            }
        }
    }

    // --- ROO 5.36 TeV pi0 ---
    TFile *ROO_nico = TFile::Open("Nico_results/LightMesonsInLightIons.root" ,"read");
    if (!ROO_nico || ROO_nico->IsZombie()) {
        std::cerr << "ERROR: could not open LightMesonsInLightIons.root" << std::endl;
    } else {
        TDirectory* dir_OO = (TDirectory*)ROO_nico->Get("5.36TeV_OO/Pi0/ROO");
        if (!dir_OO) {
            std::cerr << "ERROR: directory '5.36TeV_OO/Pi0/ROO' not found!" << std::endl;
        } else {
            TGraphAsymmErrors* g_ROO_stat = dynamic_cast<TGraphAsymmErrors*>(dir_OO->Get("NuclearModification_Stat"));
            TGraphAsymmErrors* g_ROO_sys  = dynamic_cast<TGraphAsymmErrors*>(dir_OO->Get("NuclearModification_Sys"));
            if (!g_ROO_stat || !g_ROO_sys) {
                std::cerr << "ERROR: graphs 'NuclearModification_Stat' or 'NuclearModification_Sys' not found!" << std::endl;
            } else {
                Plot_ROO.New(g_ROO_stat,"#pi^{0} OO, #sqrt{s_{NN}} = 5.36 TeV - ALICE Preliminary", 20, linethickness-2, kPink+1, "ep");
                Plot_ROO.New(g_ROO_sys,"", 20, linethickness-2, kPink+1, "Box");
            }
        }
    }

    double relNormUnc_OO = 0.1355765; // ~13.6%
    double relNormUnc_pPb = 0.0556057551; // ~5.6%;
    double yCenter = 1.0;

    // Nico OO
    double xMin_OO_Nico = 0.7;
    double xMax_OO_Nico = 0.9;
    double yLow_OO  = yCenter - relNormUnc_OO;
    double yHigh_OO = yCenter + relNormUnc_OO;

    double x_Nico[2] = {xMin_OO_Nico, xMax_OO_Nico};
    double y_Nico[2] = {yCenter, yCenter};
    double ex_Nico[2] = {0., 0.};
    double ey_Nico[2] = {relNormUnc_OO, relNormUnc_OO};

    TGraphAsymmErrors *gNorm_Nico = new TGraphAsymmErrors(2, x_Nico, y_Nico, ex_Nico, ex_Nico, ey_Nico, ey_Nico);
    gNorm_Nico->SetFillColor(kPink+1);  
    gNorm_Nico->SetLineColor(kPink+1);
    gNorm_Nico->SetLineWidth(1);

    Plot_ROO.New(gNorm_Nico, "", -1, linethickness-2, kPink+1, "E3");

    // Nico pPb
    double xMin_pPb = 0.5;
    double xMax_pPb = 0.7;

    double x_pPb[2] = {xMin_pPb, xMax_pPb};
    double y_pPb[2] = {yCenter, yCenter};
    double ex_pPb[2] = {0., 0.};
    double ey_pPb[2] = {relNormUnc_pPb, relNormUnc_pPb};

    TGraphAsymmErrors *gNorm_pPb = new TGraphAsymmErrors(2, x_pPb, y_pPb, ex_pPb, ex_pPb, ey_pPb, ey_pPb);
    gNorm_pPb->SetFillColor(kGreen-8); 
    gNorm_pPb->SetLineColor(kGreen-8);
    gNorm_pPb->SetLineWidth(1);

    Plot_ROO.New(gNorm_pPb, "", -1, linethickness-2, kGreen-8, "E3");

    // Mary
    double xMin_OO_Mary = 0.9;
    double xMax_OO_Mary = 1.1;

    double x_Mary[2] = {xMin_OO_Mary, xMax_OO_Mary};
    double y_Mary[2] = {yCenter, yCenter};
    double ex_Mary[2] = {0., 0.};
    double ey_Mary[2] = {relNormUnc_OO, relNormUnc_OO};

    TGraphAsymmErrors *gNorm_Mary = new TGraphAsymmErrors(2, x_Mary, y_Mary, ex_Mary, ex_Mary, ey_Mary, ey_Mary);
    gNorm_Mary->SetFillColor(kBlack); 
    gNorm_Mary->SetLineColor(kBlack);
    gNorm_Mary->SetLineWidth(1);

    Plot_ROO.New(gNorm_Mary, "", -1, linethickness, kBlack, "E3");





    Plot_ROO.New(hROO,"#omega OO, #sqrt{s_{NN}} = 5.36 TeV - this work", 20, linethickness, kBlack, "ep");

    Plot_ROO.SetMargins(0.12, 0.1, 0.05, 0.025);
    Plot_ROO.SetLegend(0.12, 0.45, 0.15, 0.4);
    Plot_ROO.SetAxisLabel("#bf{#it{p}_{T} (GeV/c)}", "#bf{R_{OO}}");
    Plot_ROO.SetAxisRange(0.,20.,0.,1.6);
    Plot_ROO.NewLine(0.,1.,20.,1.,2,linethickness);
    // Plot speichern
    Plot_ROO.Plot("Figures/ROO.png");
    Plot_ROO.Plot("Figures/ROO.pdf");
    
    std::cout<<"Done Plotting"<<std::endl;

    

    TGraph *gRooVsPt = new TGraph(ptCenters.size(), ptCenters.data(), rooValues.data());
    gRooVsPt->SetName("roo_vs_pt");
    gRooVsPt->SetTitle("Roo value vs p_{T};p_{T} (GeV/c);Roo value");

    TFile *outFile = new TFile("roo_vs_pt.root", "RECREATE");
    gRooVsPt->Write();
    hROO->Write();
    outFile->Close();


    frec_OO->Close();
    frec_pp->Close();
    input_eff_OO->Close();
    input_eff_pp->Close();

    std::cout.rdbuf(originalCoutBuffer);
}
