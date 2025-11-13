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

void FillCorrPlot(TH1* h, Plotting1D& plot1d, double xMin, double xMax, const char* name, Color_t color, double linethickness) {
    int binMin = h->FindBin(xMin);
    int binMax = h->FindBin(xMax);

    for (int i = 1; i <= h->GetNbinsX(); i++) {
        if (i < binMin || i > binMax) {
            h->SetBinContent(i, 0);
            h->SetBinError(i, 0);
        }
    }

    double actualLow = h->GetBinLowEdge(binMin);
    double actualUp  = h->GetBinLowEdge(binMax) + h->GetBinWidth(binMax);

    std::cout << "Histogram: " << name << std::endl;
    std::cout << "Requested range: " << xMin << " – " << xMax << std::endl;
    std::cout << "Effective bin range: " << actualLow << " – " << actualUp << std::endl;
    std::cout << "Bin range (indices): " << binMin << " – " << binMax << std::endl;

    for (int i = 1; i <= h->GetNbinsX(); i++) {
        double low  = h->GetBinLowEdge(i);
        double up   = low + h->GetBinWidth(i);
        double cent = h->GetBinCenter(i);
        double val  = h->GetBinContent(i);

        std::cout << " Bin " << i
                  << " : [" << low << ", " << up << "]"
                  << " center=" << cent
                  << " content=" << val;

        if (i >= binMin && i <= binMax) std::cout << "   <-- kept";
        else                            std::cout << "   <-- zeroed";

        std::cout << std::endl;
    }
    
    plot1d.New(h, name, 20, linethickness, color, "ep"); 
}




void plot_omegacorr(const TString& datacorr_OO_file, 
                    const TString& datacorr_pp_file, 
                    int Nevt_OO, int Nevt_ppm,
                    double xMin_OO, double xMax_OO,
                    double xMin_pp, double xMax_pp,
                    const TString& output_name,
                    double linethickness) {

    std::ofstream logfile("logfiles/logfile_plot_corr.txt");
    std::streambuf* originalCoutBuffer = std::cout.rdbuf();
    std::cout.rdbuf(logfile.rdbuf());

    TFile* fcorr_OO  = TFile::Open(datacorr_OO_file, "READ");
    TFile* fcorr_pp = TFile::Open(datacorr_pp_file, "READ");

    if (!fcorr_OO || !fcorr_pp) { std::cerr << "Fehler beim Öffnen der Dateien!" << std::endl;}

    TH1D* hDataCorr_OO = (TH1D*) fcorr_OO->Get("hCrosssec");
    TH1D* hDataCorr_pp = (TH1D*) fcorr_pp->Get("hCrosssec");

    if (!hDataCorr_OO || !hDataCorr_pp) {std::cerr << "Fehler beim Lesen der Histos!" << std::endl;}

    // Scale to get nb
    hDataCorr_OO->Scale(1e9);
    hDataCorr_pp->Scale(1e9);

    TString legendtext = TString::Format(
        "#bf{ALICE} - this work;#sqrt{#it{s}} = 5.36 TeV;"
        "#omega#rightarrow#pi^{+}#pi^{-}#pi^{0}, #pi^{0} rec w/ EMC"
    );

    Plotting1D Plot1d;

    FillCorrPlot(hDataCorr_OO, Plot1d, xMin_OO, xMax_OO, "OO", kPink+10,linethickness);
    FillCorrPlot(hDataCorr_pp, Plot1d, xMin_pp, xMax_pp, "pp", kGreen+2,linethickness);

    std::cout<<"pp: BinContent(5): "<<hDataCorr_pp->GetBinContent(5)<<", pp: BinCenter(5): "<<hDataCorr_pp->GetBinCenter(5)<<std::endl;
    std::cout<<"pp: BinContent(6): "<<hDataCorr_pp->GetBinContent(6)<<", pp: BinCenter(6): "<<hDataCorr_pp->GetBinCenter(6)<<std::endl;
    std::cout<<"pp: BinContent(6): "<<hDataCorr_pp->GetBinContent(7)<<", pp: BinCenter(6): "<<hDataCorr_pp->GetBinCenter(7)<<std::endl;
    std::cout<<"pp: BinContent(6): "<<hDataCorr_pp->GetBinContent(8)<<", pp: BinCenter(6): "<<hDataCorr_pp->GetBinCenter(8)<<std::endl;

    Plot1d.SetAxisLabel("#bf{#it{p}_{T} (GeV/#it{c})}", "#bf{#it{E}#frac{d^{3} #sigma}{d#it{p}^{3}} (nbGeV^{-2}#it{c}^{3})}");
    Plot1d.SetMargins(0.12, 0.15, 0.05, 0.025);
    Plot1d.SetLegend(0.6, 0.85, 0.38, 0.5);
    Plot1d.NewLatex(0.6, 0.6, legendtext, StdTextSize, 0.04);   
    Plot1d.SetAxisRange(0, 20., 3, 1.*1e5);
        
    Plot1d.Plot("Figures/corrected/" + output_name + ".png", false, true);
    Plot1d.Plot("Figures/corrected/" + output_name + ".pdf", false, true);

    std::cout << "Done Plotting" << std::endl;

    fcorr_OO->Close();
    fcorr_pp->Close();

    std::cout.rdbuf(originalCoutBuffer);
    logfile.close();    
}
