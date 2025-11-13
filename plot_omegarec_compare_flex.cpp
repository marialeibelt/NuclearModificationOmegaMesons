#include <iostream>
#include <fstream>
#include <cmath>
#include "TCanvas.h"
#include "TAxis.h"  
#include "TLegend.h" 
#include "TF1.h"
#include "TGraph.h"
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include "/alf/data/mleibelt/LogPlotLib/Plotting.h"



Color_t colorHist_list[3] = {kRed, kAzure+6, kMagenta};
Color_t colorBG_list[3] = {kRed+2, kAzure+2, kMagenta+2};
Color_t colorSignal_list[3] = {kPink-4, kGreen-7, kBlue-9};
Color_t colorSignalFit_list[3] = {kPink+10, kGreen+2, kViolet+1};

TString dataset[5] = {"small dataset","medium dataset","Test","new","large"};
TString collisionsystem[3] = {"pp","Pb-Pb","OO"};

double hline_grid[13] = {2.*1e3, 2.*1e3, 2.*1e3,2.*1e3, 2.*1e3, 2.*1e3, 2.*1e3, 2.*1e3, 2.*1e3, 1.*1e3, 1.*1e3, 0.4*1e3};
double hline_SignalOnlyGrid[13] = {100., 200., 200.,100., 200., 200., 200., 200., 200., 100., 100., 40.};

double yindividual1[13] = {-2*1e3, -2*1e3, -2*1e3,  -2*1e3,  -2*1e3, -1.5*1e3, -1.0*1e3, -0.7*1e3, -0.6*1e3, -0.5*1e3, -0.3*1e3, -200};
double yindividual2[13] = {25*1e3, 25*1e3, 25*1e3, 30.*1e3, 25.*1e3, 14.0*1e3,  8.0*1e3,  5.0*1e3,  5.0*1e3,  2.5*1e3,  1.8*1e3,  700};
double y_sig_individual1[13] = {-2*1e3, -2*1e3, -2*1e3,  -2*1e3,  -2*1e3, -1.5*1e3, -1.0*1e3, -0.7*1e3, -0.6*1e3, -0.5*1e3, -0.3*1e3, -200};
double y_sig_individual2[13] = {25*1e3, 25*1e3, 25*1e3, 30.*1e3, 25.*1e3, 14.0*1e3,  8.0*1e3,  5.0*1e3,  5.0*1e3,  2.5*1e3,  1.8*1e3,  700};

double Momega_EMCal = 0.7822; //in GeV/c2
double sigmaomega_EMCal = 0.014; //in GeV/c2
double low_cut  = (Momega_EMCal - 3 * sigmaomega_EMCal);  // 0.7402 in GeV/c
double high_cut = (Momega_EMCal + 3 * sigmaomega_EMCal);  // 0.8242 in GeV/c
double pT_bins[13] = {1.8, 2.2, 2.6, 3.2, 4.0,5.0,6.0,7.0,8.0,10.0,12.0,16.0,20.0};
int npT = sizeof(pT_bins) / sizeof(pT_bins[0]); // will be 13
double Momega_PDG = 0.78265; //in GeV/c2 /pm 0.16
double sigmaomega_PDG = 0.00849; ////in GeV/c2 /pm 0.00008
double xlow = 0.6;
double xup = 1.;

//labels
//Mgg,pTgg,M3pi,pT3pi,Ncoll,Ntracks,zvtx,NV0s,Nclust,Ngg,Nggsel,Nhnm,Nhnmsel,pTPC
//0  ,1   ,2   ,3    ,4    ,5      ,6   ,7   ,8     ,9  ,10    ,11  ,12     ,13  

//p ,pT,eta,phi,dE/dx,ppi+TPC,ppi-TPC,ppi+,Npi+,pTpi+,ppi+TPC-ppi+,ppi+TPC-ppi+/ppi+
//14,15,16 ,17 ,18   ,19     ,20     ,21  ,22  ,23   ,24          ,25

//nsigmapi+TPC,nsigmapi+TOF,nsigmapi+comb,nsigmapi+ITS,DCAxy,DCAz,TPCsharedclu,
//26          ,27          ,28           ,29          ,30   ,31  ,32

//TPCcrossedrows,TPCfind/crossedrows,TPCclu,ppi-,Npi-,pTpi-,ppi-TPC-ppi-,ppi-TPC-ppi-/ppi-
//33            ,34                 ,35    ,36  ,37  ,38   ,39

//ppi-TPC-ppi-/ppi-,nsigmapi-TPC,nsigmapi-TOF,nsigmapi-comb,nsigmapi-ITS,N3pi,Counts,Momega
//40               ,41          ,42          ,43           ,44          ,45  ,46    ,47    

//Chi2/NDf,dN/dpT,
//48      ,49    , 
TString axislabels[50] = {"#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})}", "#bf{#it{p}^{#gamma#gamma}_{T} (GeV/#it{c})}",
                        "#bf{#it{M}^{#pi^{+}#pi^{-}#gamma#gamma} (GeV/#it{c}^{2})}", "#bf{#it{p}^{#pi^{+}#pi^{-}#gamma#gamma}_{T} (GeV/#it{c})}",
                        "#bf{#it{N}_{collisions}}","#bf{#it{N}_{tracks}}","#bf{#it{z}_{vtx}} (cm)","#bf{#it{N}_{V0s}}","#bf{#it{N}_{clusters}}",
                        "#bf{#it{N}_{#gamma#gamma}}","#bf{#it{N}_{#gamma#gamma}^{selected}}","#bf{#it{N}_{HNM}}","#bf{#it{N}_{HNM}^{selected}}",
                        "#bf{#it{p}_{TPC} (GeV/#it{c})}","#bf{#it{p} (GeV/#it{c})}","#bf{#it{p}_{T} (GeV/#it{c})}","#bf{#eta}","#bf{#phi}",
                        "#bf{dE/dx}", "#bf{#it{p}^{#pi^{+}}_{TPC} (GeV/#it{c})}", "#bf{#it{p}^{#pi^{-}}_{TPC} (GeV/#it{c})}",
                        "#bf{#it{p}^{#pi^{+}} (GeV/#it{c})}","#bf{#it{N}^{#pi^{+}}}","#bf{#it{p}^{#pi^{+}}_{T} (GeV/#it{c})}",
                        "#bf{#it{p}^{#pi^{+}}_{TPC} - #it{p}^{#pi^{+}} (GeV/#it{c})}","#bf{#it{p}^{#pi^{+}}_{TPC} - #it{p}^{#pi^{+}} / p^{#pi^{+}}}",
                        "#bf{n#sigma^{#pi^{+}}_{TPC}}","#bf{n#sigma^{#pi^{+}}_{TOF}}","#bf{n#sigma^{#pi^{+}}_{comb}}","#bf{n#sigma^{#pi^{+}}_{ITS}}",
                        "#bf{DCA_{xy}}","#bf{DCA_{z}}","#bf{TPC Shared Clusters}","#bf{TPC Crossed Rows}","#bf{TPC Findable/CrossedRows}","#bf{TPC Clusters}",
                        "#bf{#it{p}^{#pi^{-}} (GeV/#it{c})}","#bf{#it{N}^{#pi^{-}}}","#bf{#it{p}^{#pi^{+}}_{T} (GeV/#it{c})}",
                        "#bf{#it{p}^{#pi^{-}}_{TPC} - #it{p}^{#pi^{-}} (GeV/#it{c})}","#bf{#it{p}^{#pi^{-}}_{TPC} - #it{p}^{#pi^{-}} / p^{#pi^{-}}}",
                        "#bf{n#sigma^{#pi^{-}}_{TPC}}","#bf{n#sigma^{#pi^{-}}_{TOF}}","#bf{n#sigma^{#pi^{-}}_{comb}}","#bf{n#sigma^{#pi^{-}}_{ITS}}",
                        "#bf{N^{#pi^{+}#pi^{-}#gamma#gamma}}","#bf{Counts}","#bf{#it{M}_{#omega} (GeV/#it{c}^{2})}","#bf{Chi2/NDf}",
                        "#bf{#frac{d#it{N}}{d#it{p}_{T}} (#it{c}/GeV) }"};


// 1D Comparison plots
std::vector<TString> hist_names = {
    "h_chi2ndf", "h_chi2ndf_signalfit", "h_mean_signalfit",
    "h_sigma_signalfit", "h_FWHM_signalfit", "h_omegacounts"
};
std::vector<TString> ylabels = {
    axislabels[48], axislabels[48], axislabels[47],
    "#bf{Width (GeV/#it{c}^{2})}", "#bf{FWHM (GeV/#it{c}^{2})}", axislabels[49]
};
std::vector<TString> filenames = {
    "chi2ndf_comparison", "chi2ndf_signalfit_comparison",
    "mean_signalfit_comparison", "sigma_signalfit_comparison",
    "width_FWHM_comparison", "countsperbinwidth_comparison"
};




Double_t FWHM_Func(Double_t sigma){
  return 2.0 * sqrt(2.0 * log(2.0)) * sigma;
}



void AddHistogramWithFits(TH1D* h, PlottingGrid& grid, const char* labelIn,
                             double signalScaleFactor = 1.0, int floop = 0, TString sigorproj = "huhu", double linethickness = 1.) {
    //h->Print("all");
    TString label(labelIn);
    TH1D* histToPlot = nullptr;

    if(sigorproj == "sig") {
        histToPlot = (TH1D*)h->Clone(Form("%s_Signal_scaled", h->GetName()));
        histToPlot->Scale(signalScaleFactor);
    }
    else if(sigorproj == "proj") {
        histToPlot = (TH1D*)h->Clone(Form("%s_proj", h->GetName()));
    }
    histToPlot->GetListOfFunctions()->Clear();

    grid.New(histToPlot,
                       (sigorproj == "sig") ? Form("%s Signal x %.1f", label.Data(), signalScaleFactor) : label,
                       1, linethickness, (sigorproj == "sig") ? colorSignal_list[floop] : colorHist_list[floop], "ep");

    TList* funcs = h->GetListOfFunctions();
    if (!funcs) {std::cout << "No functions found in histogram: " << h->GetName() << std::endl;return;}

    for (int i = 0; i < funcs->GetSize(); ++i) {
        TObject* obj = funcs->At(i);
        if (!obj) {std::cout<<"Skipping because obj in AddHistogramWithFits not found"; continue;}

        if (auto* fcast = dynamic_cast<TF1*>(obj)) {
            TF1* fclone = (TF1*)fcast->Clone(Form("%s_%s_clone", fcast->GetName(), label.Data()));
            TString fname = fclone->GetName();
            fclone->SetNpx(1000);  // 1000 points instead of 100
            Color_t funccolor = kBlack;
            TString funcname = "Function";

            if (fname.Contains("gausFit")) {
                funccolor = colorSignalFit_list[floop];
                if (fclone->GetNpar() > 0) {
                    fclone->SetParameter(0, fclone->GetParameter(0) * signalScaleFactor);
                }
                funcname = Form("%s Signal fit x %.1f", label.Data(), signalScaleFactor);
            }
            else if (fname.Contains("BG")) {
                funccolor = colorBG_list[floop];
                funcname = Form("%s Fitted BG", label.Data());
            }

            grid.NewFunc(fclone, funcname, 1, linethickness, funccolor, "l");
        }
    }
}

void AddSignalOnly(TH1D* h, PlottingGrid& grid, const char* labelIn, int floop = 0, double linethickness = 1.) {

    TString label(labelIn);

    TH1D* histToPlot = (TH1D*)h->Clone(Form("%s_Signal", h->GetName()));
    histToPlot->GetListOfFunctions()->Clear();
    grid.New(histToPlot, Form("%s Signal", label.Data()), 1, linethickness, colorSignal_list[floop], "ep");

    TList* funcs = h->GetListOfFunctions();
    if (!funcs || funcs->GetSize() == 0) { std::cout << "No TF1 functions found in histogram: " << h->GetName() << std::endl;return;}

    for (int j = 0; j < funcs->GetSize(); ++j) {
        TF1* fcast = dynamic_cast<TF1*>(funcs->At(j));
        if (!fcast) continue;

        TString fname(fcast->GetName());
        if (!fname.Contains("gausFit")) continue;

        TF1* fclone = (TF1*)fcast->Clone(Form("%s_%s_sigclone", fname.Data(), label.Data()));
        fclone->SetNpx(1000);

        grid.NewFunc(fclone, Form("%s Signal fit", label.Data()), 1, linethickness, colorSignalFit_list[floop], "l");
    }
    std::cout<<""<<std::endl;
}


void AddHistogramWithFits_1D(TH1D* h, Plotting1D& plotindividual, const char* labelIn,
                             double signalScaleFactor = 1.0, int floop = 0, TString sigorproj = "huhu", double linethickness = 1.) {
    //h->Print("all");
    TString label(labelIn);
    TH1D* histToPlot = nullptr;

    if(sigorproj == "sig") {
        histToPlot = (TH1D*)h->Clone(Form("%s_Signal_scaled", h->GetName()));
        histToPlot->Scale(signalScaleFactor);
    }
    else if(sigorproj == "proj") {
        histToPlot = (TH1D*)h->Clone(Form("%s_proj", h->GetName()));
    }
    histToPlot->GetListOfFunctions()->Clear();

    plotindividual.New(histToPlot,
                       (sigorproj == "sig") ? Form("%s Signal x %.1f", label.Data(), signalScaleFactor) : label,
                       1, linethickness, (sigorproj == "sig") ? colorSignal_list[floop] : colorHist_list[floop], "ep");

    TList* funcs = h->GetListOfFunctions();
    if (!funcs) {std::cout << "No functions found in histogram: " << h->GetName() << std::endl;return;}

    for (int i = 0; i < funcs->GetSize(); ++i) {
        TObject* obj = funcs->At(i);
        if (!obj) {std::cout<<"Skipping because obj in AddHistogramWithFits_1D not found";continue;}

        if (auto* fcast = dynamic_cast<TF1*>(obj)) {
            TF1* fclone = (TF1*)fcast->Clone(Form("%s_%s_clone", fcast->GetName(), label.Data()));
            fclone->SetNpx(1000);  // 1000 points instead of 100
            TString fname = fclone->GetName();
            Color_t funccolor = kBlack;
            TString funcname = "Function";

            if (fname.Contains("gausFit")) {
                funccolor = colorSignalFit_list[floop];
                if (fclone->GetNpar() > 0) {
                    fclone->SetParameter(0, fclone->GetParameter(0) * signalScaleFactor);
                }
                funcname = Form("%s Signal fit x %.1f", label.Data(), signalScaleFactor);
            }
            else if (fname.Contains("BG")) {
                funccolor = colorBG_list[floop];
                funcname = "Fitted BG";
            }

            plotindividual.New(fclone, funcname, 1, linethickness, funccolor, "l");
        }
    }
}

void AddSignalOnly_1D(TH1D* h, Plotting1D& plotindividual, const char* labelIn, int floop = 0, double linethickness = 1.) {

    TString label(labelIn);

    TH1D* histToPlot = (TH1D*)h->Clone(Form("%s_Signal", h->GetName()));
    histToPlot->GetListOfFunctions()->Clear();
    plotindividual.New(histToPlot, Form("%s Signal", label.Data()), 1, linethickness, colorSignal_list[floop], "ep");

    TList* funcs = h->GetListOfFunctions();
    if (!funcs || funcs->GetSize() == 0) { std::cout << "No TF1 functions found in histogram: " << h->GetName() << std::endl;return;}

    for (int j = 0; j < funcs->GetSize(); ++j) {
        TF1* fcast = dynamic_cast<TF1*>(funcs->At(j));
        if (!fcast) continue;

        TString fname(fcast->GetName());
        if (!fname.Contains("gausFit")) continue;

        TF1* fclone = (TF1*)fcast->Clone(Form("%s_%s_sigclone", fname.Data(), label.Data()));
        fclone->SetNpx(1000);

        plotindividual.New(fclone, Form("%s Signal fit", label.Data()), 1, linethickness, colorSignalFit_list[floop], "l");
    }
    std::cout<<""<<std::endl;
}




// Function to show results from different files in one figure
void plot_omegarec_compare_flex(
    const std::vector<TString>& file_vec,
    const std::vector<int>& datasetnumbers,
    const std::vector<int>& collisionsystemnumbers,
    const TString& output_name,
    const std::vector<int>& goodbins,
    double linethickness = 1.){

    std::ofstream logfile("logfiles/logfile_plot_compare_flex.txt");
    std::streambuf* originalCoutBuffer = std::cout.rdbuf();
    std::cout.rdbuf(logfile.rdbuf());

    std::vector<TFile*> input_files;
    std::array<double,13> yindividual1;
    std::array<double,13> yindividual2;
    std::array<double,13> y_sig_individual1;
    std::array<double,13> y_sig_individual2;
    int n_bins = npT - 1;
    int n_pTstart = 0;
    int n_end = n_bins-2;

    for (const auto& fname : file_vec) {
        std::cout << "fname: " << fname << std::endl;

        if (fname.Contains("OO")) {
            yindividual1      = {-2000, -2000, -2000, -2000, -2000, -1500, -1000, -700, -600, -500, -300, -200, -200};
            yindividual2      = {25000, 25000, 25000, 4.8*1e6, 1.3*1e6, 4.5*1e5, 1.5*1e5, 5.3*1e4, 4.2*1e4, 1.7*1e4, 1.2*1e4, 1*1e4, 1*1e4};
            y_sig_individual1 = {-1, -1, -1, -7*1e3,   -5*1e3, -1.6*1e3, -0.6*1e3, -0.45*1e3, -0.5*1e3, -200, -190, -120, -120};
            y_sig_individual2 = { 1,  1,  1,  2.2*1e4,  1.4*1e4,    5*1e3,  2.5*1e3,   1.5*1e3, 1.18*1e3,  460,  490,  270,  270};
            n_pTstart = 3;
            n_end = n_bins-2;
            // n_pTstart = 7;
            // n_end = n_bins-3;
        }
        else if (fname.Contains("pp")) {
            yindividual1      = {-2000, -2000, -2000, -2000, -2000, -1500, -1000, -700, -600, -500, -300, -200, -200};
            yindividual2      = {25000, 25000, 25000, 30000, 25000, 14000, 8000, 5000, 5000, 2500, 1800, 700, 700};
            y_sig_individual1 = {-1, -1, -1, -0.3*1e3,   -0.3*1e3, -0.22*1e3, -150, -170, -120, -95, -70, -39, -39};
            y_sig_individual2 = { 1,  1,  1,  1.3*1e3,   1.25*1e3,  1.05*1e3,  550,  400,  320, 220, 150,  80,  80};
            n_pTstart = 3;
            n_end = n_bins-2;
        }
        input_files.push_back(new TFile(fname, "READ"));
    }

    
    TString legendtext = TString::Format(
        "#bf{ALICE} - this work;%s, #sqrt{#it{s}} = 5.36 TeV;#omega#rightarrow#pi^{+}#pi^{-}#pi^{0}, #pi^{0} rec w/ EMC; Pol3 + Exclusion",
        collisionsystem[collisionsystemnumbers[0]].Data()
    );

    PlottingGrid grid;
    PlottingGrid SignalOnlyGrid;
    grid.SetAxisLabel(axislabels[2], axislabels[46]);
    SignalOnlyGrid.SetAxisLabel(axislabels[2], axislabels[46]);

    
    for (int i = n_pTstart; i < n_end; ++i) {
        double pTlow = pT_bins[i];
        double pThigh = pT_bins[i+1];

        TString legendtext2 = TString::Format(
            "#bf{ALICE} - this work;%s, #sqrt{#it{s}} = 5.36 TeV;%.1f #leq p_{T} (GeV/#it{c}) < %.1f;#omega#rightarrow#pi^{+}#pi^{-}#pi^{0}, #pi^{0} rec w/ EMC",
            collisionsystem[collisionsystemnumbers[0]].Data(),
            pTlow,
            pThigh
        );

        Plotting1D IndividualPlot;
        Plotting1D SignalOnlyIndividualPlot;
        IndividualPlot.SetAxisLabel(axislabels[2], axislabels[46]);
        SignalOnlyIndividualPlot.SetAxisLabel(axislabels[2], axislabels[46]);

        for (size_t f = 0; f < input_files.size(); ++f) {
            TH1D* proj_bg  = static_cast<TH1D*>(input_files[f]->Get(Form("proj_pT_%d", i)));
            TH1D* sig_gaus = static_cast<TH1D*>(input_files[f]->Get(Form("signal_clone_%d", i)));

            
            if (!sig_gaus) {std::cout << "Histogram signal_clone_" << i << " not found for f = " << f << std::endl;}
            if (!proj_bg) {std::cout << "Histogram proj_bg" << i << " not found for f = " << f << std::endl;}

            TString infilename = input_files[f]->GetName();
            std::cout<<"infilename: "<<infilename<<std::endl;
            TString label;
            if (infilename.Contains("MC")) {label = "MC";}
            else if (infilename.Contains("Data")) {label = "Data";}

            AddHistogramWithFits(proj_bg, grid, label, 5.0, f,"proj",linethickness);
            AddHistogramWithFits(sig_gaus, grid, label, 5.0, f,"sig",linethickness);
            AddSignalOnly(sig_gaus, SignalOnlyGrid, label,f,linethickness);
            AddHistogramWithFits_1D(proj_bg, IndividualPlot, label, 5.0, f,"proj",linethickness);
            AddHistogramWithFits_1D(sig_gaus, IndividualPlot, label, 5.0, f,"sig",linethickness);
            AddSignalOnly_1D(sig_gaus, SignalOnlyIndividualPlot, label, f,linethickness);
        }
        grid.NewLine(low_cut, 0., low_cut, hline_grid[i], 2, linethickness, 1., "Integration range");
        grid.NewLine(high_cut, 0., high_cut, hline_grid[i], 2, linethickness, 1., "");
        grid.NewLine(xlow, 0., low_cut, 0., 1, linethickness, 1., "");
        grid.NewLine(high_cut, 0., xup, 0., 1, linethickness, 1., "");
        grid.SetMargins(0.12, 0.1, 0.04, 0.13);
        grid.SetLegend(0.15, 0.5, 0., 0.5);
        grid.NewLatex(0.17, 0.8, legendtext, StdTextSize * 0.3, 0.08);
        grid.NextPad(Form("%.1f < #it{p}_{T} < %.1f", pT_bins[i], pT_bins[i + 1]));

        SignalOnlyGrid.NewLine(low_cut, 0., low_cut, hline_SignalOnlyGrid[i], 2, linethickness, 1., "Integration range");
        SignalOnlyGrid.NewLine(high_cut, 0., high_cut, hline_SignalOnlyGrid[i], 2, linethickness, 1., "");
        SignalOnlyGrid.NewLine(xlow, 0., low_cut, 0., 1, linethickness, 1., "");
        SignalOnlyGrid.NewLine(high_cut, 0., xup, 0., 1, linethickness, 1., "");
        SignalOnlyGrid.SetMargins(0.12, 0.1, 0.04, 0.13);
        SignalOnlyGrid.SetLegend(0.15, 0.5, 0., 0.5);
        SignalOnlyGrid.NewLatex(0.17, 0.8, legendtext, StdTextSize * 0.3, 0.08);
        SignalOnlyGrid.NextPad(Form("%.1f < #it{p}_{T} < %.1f", pT_bins[i], pT_bins[i + 1]));


        IndividualPlot.NewLine(low_cut, 0., low_cut, hline_grid[i], 2, linethickness, 1., "");
        IndividualPlot.NewLine(high_cut, 0., high_cut, hline_grid[i], 2, linethickness, 1., "");
        IndividualPlot.NewLine(xlow, 0., low_cut, 0., 1, linethickness, 1., "");
        IndividualPlot.NewLine(high_cut, 0., xup, 0., 1, linethickness, 1., "");
        IndividualPlot.SetMargins(0.12, 0.12, 0.05, 0.025);
        IndividualPlot.SetLegend(0.17, 0.45, 0.6, 0.93, false, 0.8);
        IndividualPlot.NewLatex(0.65, 0.87, legendtext2, StdTextSize * 0.8, 0.04);
        IndividualPlot.SetAxisRange(xlow,xup);
        IndividualPlot.SetAxisRange(xlow,xup,yindividual1[i],yindividual2[i]);

        IndividualPlot.Plot(Form("Figures/%s/individuals/png/individual_pT_%d_%s.png", output_name.Data(), i, output_name.Data()));
        IndividualPlot.Plot(Form("Figures/%s/individuals/pdf/individual_pT_%d_%s.pdf", output_name.Data(), i, output_name.Data()));


        SignalOnlyIndividualPlot.SetMargins(0.12, 0.12, 0.05, 0.025);
        SignalOnlyIndividualPlot.SetLegend(0.17, 0.45, 0.7, 0.93, false, 0.8);
        SignalOnlyIndividualPlot.NewLatex(0.65, 0.85, legendtext2, StdTextSize * 0.8, 0.04);
        SignalOnlyIndividualPlot.SetAxisRange(xlow,xup,y_sig_individual1[i],y_sig_individual2[i]);

        SignalOnlyIndividualPlot.Plot(Form("Figures/%s/individuals/png/signalonly_individual_pT_%d_%s.png", output_name.Data(), i, output_name.Data()));
        SignalOnlyIndividualPlot.Plot(Form("Figures/%s/individuals/pdf/signalonly_individual_pT_%d_%s.pdf", output_name.Data(), i, output_name.Data()));
    }
    //grid.Plot(Form("Figures/%s/png/all_pT_projections_%s_grid.png", output_name.Data(), output_name.Data()));
    grid.Plot(Form("Figures/%s/pdf/all_pT_projections_%s_grid.pdf", output_name.Data(), output_name.Data()));
    //SignalOnlyGrid.Plot(Form("Figures/%s/png/all_pT_projections_%s_signalonly_grid.png", output_name.Data(), output_name.Data()));
    SignalOnlyGrid.Plot(Form("Figures/%s/pdf/all_pT_projections_%s_signalonly_grid.pdf", output_name.Data(), output_name.Data()));





    for (size_t h = 0; h < hist_names.size(); ++h) {
        Plotting1D Plot1d;

        TString legendtext2 = TString::Format(
            "#bf{ALICE} - this work;%s, #sqrt{#it{s}} = 5.36 TeV;#omega#rightarrow#pi^{+}#pi^{-}#pi^{0}, #pi^{0} rec w/ EMC",
            collisionsystem[collisionsystemnumbers[0]].Data()
        );

        for (size_t f = 0; f < input_files.size(); ++f) {
            TH1* h1d_orig = dynamic_cast<TH1*>(input_files[f]->Get(hist_names[h]));
            if (!h1d_orig) continue;
            std::cout << "Histogram " << h1d_orig->GetName()
          << " has " << h1d_orig->GetNbinsX() << " bins" << std::endl;
            TString infilename = input_files[f]->GetName();
            std::cout<<"infilename: "<<infilename<<std::endl;

            TString label;

            //-----------------Alignment check-----------------------
            // if(f==0){label = "Run 2 Alignment";}
            // else if(f==2){label = "Run 3 Alignment";} //watch out for right label!!! different for 2vsno und 2vs3vsno
            // else if(f==1){label = "No Alignment";}
            //-------------------------------------------------------

            if (infilename.Contains("MC")) {label = "MC";}
            else if (infilename.Contains("Data")) {label = "Data";}

            TH1* h1d = h1d_orig;

            if ((hist_names[h].Contains("mean") ||
                hist_names[h].Contains("sigma") ||
                hist_names[h].Contains("FWHM") ||
                hist_names[h].Contains("counts")) &&
                !goodbins.empty()) 
            {
                TString newname = TString::Format("%s_goodbins_%zu", h1d_orig->GetName(), f);
                TH1* h_reduced = (TH1*) h1d_orig->Clone(newname);
                h_reduced->Reset("ICES");

                for (int gb : goodbins) {
                    if (gb < 0 || gb >= npT-1) {std::cout << "Skipping invalid goodbin: " << gb << std::endl;continue;}
                    double xmid = 0.5 * (pT_bins[gb] + pT_bins[gb+1]);
                    int bin = h1d_orig->GetXaxis()->FindBin(xmid);
                    h_reduced->SetBinContent(bin, h1d_orig->GetBinContent(bin));
                    h_reduced->SetBinError(bin, h1d_orig->GetBinError(bin));
                }
                h1d = h_reduced;
            }

            Plot1d.New(h1d, label, 20 + f, linethickness, colorSignalFit_list[f], "ep");
            Plot1d.SetMargins(0.12, 0.15, 0.05, 0.025);
            Plot1d.SetLegend(0.45, 0.65, 0.7, 0.9);
        }

        Plot1d.SetAxisLabel(axislabels[15], ylabels[h]);
        Plot1d.SetAxisRange(0., pT_bins[npT - 1]);

        if (hist_names[h].Contains("mean")) {
            Plot1d.NewLine(0., Momega_EMCal, pT_bins[npT - 1], Momega_EMCal, 1, linethickness, 1., "EMCal Run 2 ref");
            Plot1d.NewLine(0., Momega_PDG, pT_bins[npT - 1], Momega_PDG, 2, linethickness, 1., "PDG");
            //Plot1d.NewLine(0., low_cut, pT_bins[npT - 1], low_cut, 1, linethickness, kRed + 2, "Excluded region");
            //Plot1d.NewLine(0., high_cut, pT_bins[npT - 1], high_cut, 1, linethickness, kRed + 2, "");
            Plot1d.SetMargins(0.12, 0.15, 0.05, 0.025);
            Plot1d.SetLegend(0.6, 0.85, 0.15, 0.35);
            Plot1d.NewLatex(0.6, 0.45, legendtext2, StdTextSize, 0.04);
            Plot1d.SetAxisRange(0., pT_bins[npT - 1],0.763,0.786);
        }
        if (hist_names[h].Contains("sigma")) {
            Plot1d.NewLine(0., sigmaomega_EMCal, pT_bins[npT - 1], sigmaomega_EMCal, 1, linethickness, 1., "EMCal Run 2 ref");
            Plot1d.NewLine(0., sigmaomega_PDG, pT_bins[npT - 1], sigmaomega_PDG, 2, linethickness, 1., "PDG");
            Plot1d.SetMargins(0.12, 0.15, 0.05, 0.025);
            Plot1d.SetLegend(0.6, 0.85, 0.5, 0.75);
            Plot1d.NewLatex(0.6, 0.87, legendtext2, StdTextSize, 0.04);
        }
        if (hist_names[h].Contains("FWHM")) {
            Plot1d.NewLine(0., FWHM_Func(sigmaomega_EMCal), pT_bins[npT - 1], FWHM_Func(sigmaomega_EMCal), 1, linethickness, 1., "EMCal Run 2 ref");
            Plot1d.NewLine(0., FWHM_Func(sigmaomega_PDG), pT_bins[npT - 1], FWHM_Func(sigmaomega_PDG), 2, linethickness, 1., "PDG");
            Plot1d.SetMargins(0.12, 0.15, 0.05, 0.025);
            Plot1d.SetLegend(0.6, 0.85, 0.5, 0.75);
            // Plot1d.SetLegend(0.6, 0.85, 0.15, 0.35); //Alignment check
            Plot1d.NewLatex(0.6, 0.87, legendtext2, StdTextSize, 0.04);
            Plot1d.SetAxisRange(0., pT_bins[npT - 1],18*1e-3,80*1e-3);
        }
        if (hist_names[h].Contains("counts")) {
            Plot1d.SetMargins(0.12, 0.15, 0.05, 0.025);
            Plot1d.SetLegend(0.6, 0.85, 0.55, 0.75);
            Plot1d.NewLatex(0.6, 0.87, legendtext2, StdTextSize, 0.04);
        }
        Plot1d.Plot("Figures/" + output_name + "/" + filenames[h]+".png");
        Plot1d.Plot("Figures/" + output_name + "/" + filenames[h]+".pdf");
    }

    for (auto* file : input_files)
        if (file) file->Close();

    std::cout << "byebye" << std::endl;
    std::cout.rdbuf(originalCoutBuffer);
}
