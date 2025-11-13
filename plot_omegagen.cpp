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

TObject* obj = nullptr;
TString xaxislabel = "";
TString yaxislabel = "";
TString collisionsystem[3] = {"pp","Pb-Pb","OO"};
TString dataset[5] = {"small dataset","medium dataset","Test","new","large"};
double hline_grid[13] = {2.*1e3, 2.*1e3, 2.*1e3,2.*1e3, 2.*1e3, 2.*1e3, 2.*1e3, 2.*1e3, 2.*1e3, 1.*1e3, 1.*1e3, 0.4*1e3};
double hline_SignalOnlyGrid[13] = {100., 200., 200.,100., 200., 200., 200., 200., 200., 100., 100., 40.};
double Momega_EMCal = 0.7822; //in GeV/c2
double sigmaomega_EMCal = 0.014; //in GeV/c2
double low_cut  = (Momega_EMCal - 3 * sigmaomega_EMCal);  // 0.7402 in GeV/c
double high_cut = (Momega_EMCal + 3 * sigmaomega_EMCal);  // 0.8242 in GeV/c
double pT_bins[13] = {1.8, 2.2, 2.6, 3.2, 4.0,5.0,6.0,7.0,8.0,10.0,12.0,16.0,20.0};
int npT = sizeof(pT_bins) / sizeof(pT_bins[0]); // will be 13
//double pT_bins[] = {7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7,8.8, 8.9};
double Momega_PDG = 0.78265; //in GeV/c2 /pm 0.16
double sigmaomega_PDG = 0.00849; ////in GeV/c2 /pm 0.00008
int n_pTstart = -999;
double xlow = 0.6;
double xup = 1.;
int nevents = 1;

TString file[13] = {"OldOriginalTask_SmallppRefDataSet",
                "OldTask_SmallppRefDataSet", 
                "NewTask_SmallppRefDataSet",
                "NewTask_MediumppRefDataSet_Run2Alignment",
                "NewTask_MediumppRefDataSet_Run3Alignment",
                "NewTask_SmallOODataSet_Run3Alignment",
                "ppRefNotAligned",
                "ppMC_fixed","ppMC_fixed_gen","OOMC_fixed","OOMC_fixed_gen",
                "ppMC_large","ppMC_large_gen"};

TString legend_names[13] = {"Old task",
                "Old task - loose pi0 cut", 
                "New task",
                "Run 2 alignment",
                "Run 3 alignment",
                "New task-small OO dataset - run 3 alignment",
                "No Alignment",
                "MC reconstructed",
                "MC generated",
                "MC reconstructed",
                "MC generated",
                "MC reconstructed",
                "MC generated"};

Color_t Colors[2] = {kGreen+1,kPink+7};

void plot_omegagen(const std::vector<int>& file_indices,
                   const std::vector<int>& datasetnumbers,
                   const std::vector<int>& collisionsystemnumbers,
                   const TString& output_name)
{
    // --- gewünschtes pT-Binning ---
    double pT_bins[10] = {3.2, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 16.0, 20.0};
    int nbins_rec = sizeof(pT_bins) / sizeof(pT_bins[0]) - 1;

    std::ofstream logfile("logfiles/logfile_plot_gen.txt");
    std::streambuf* originalCoutBuffer = std::cout.rdbuf();
    std::cout.rdbuf(logfile.rdbuf());

    std::vector<TFile*> input_files;
    std::vector<TString> hist_names = {"Yield_TZSGUE", "h_omegacounts"}; // 0 = gen, 1 = rec
    std::vector<TString> filenames = {"NvspT"};

    // Dateien öffnen

    for (int idx : file_indices) {
        TString filename= file[idx];
        if (filename.Contains("gen")){
            TString fname = "outputs/omegagen/" + file[idx] + "_hists.root";
            std::cout << "fname: " << fname << std::endl;
            input_files.push_back(TFile::Open(fname, "read"));
        }
        else{
            TString fname = "outputs/omegarec/" + file[idx] + "_hists.root";
            std::cout << "fname: " << fname << std::endl;
            input_files.push_back(TFile::Open(fname, "read"));
        }
    }


    // Effizienz berechnen mit Rebinning
    TH1* hGen_fine = dynamic_cast<TH1*>(input_files[0]->Get(hist_names[0]));
    TH1* hRec = dynamic_cast<TH1*>(input_files[1]->Get(hist_names[1]));

    if (!hGen_fine || !hRec) {
        std::cout << "Histos nicht gefunden!" << std::endl;
        return;
    }

    hGen_fine->Sumw2();
    hRec->Sumw2();

    TH1* hGen = hGen_fine->Rebin(nbins_rec, "hGen_rebinned", pT_bins);
    
    if (nevents > 0) {
        hGen->Scale(1.0 / nevents);
        hRec->Scale(1.0 / nevents);
    }

    // MC gen vs rec plotten
    Plotting1D Plot1d;

    Plot1d.New(hGen, legend_names[file_indices[0]], 20, 2.0, Colors[0], "ep");
    Plot1d.New(hRec, legend_names[file_indices[1]], 21, 2.0, Colors[1], "ep");

    Plot1d.SetAxisLabel("#bf{#it{p}_{T} (GeV/#it{c})}", "#bf{N^{#omega}}");
    Plot1d.SetMargins(0.12, 0.15, 0.05, 0.025);
    Plot1d.SetLegend(0.5, 0.75, 0.7, 0.9);
    Plot1d.Plot("Figures/omegagen/" + output_name + "/" + filenames[0] + ".png", false, true);


    // --- Vor Division Werte ausgeben ---
    std::cout << "\n--- Vor Division ---\n";
    for (int bin = 1; bin <= nbins_rec; ++bin) {
        double pt_low  = hGen->GetXaxis()->GetBinLowEdge(bin);
        double pt_high = hGen->GetXaxis()->GetBinUpEdge(bin);
        double bin_width = pt_high - pt_low;

        double gen_counts = hGen->GetBinContent(bin);
        double rec_counts = hRec->GetBinContent(bin);

        std::cout << Form("pT: %.2f - %.2f | Gen: %.3g | Rec: %.3g | Width: %.2f",
                          pt_low, pt_high, gen_counts, rec_counts, bin_width)
                  << std::endl;
    }

    // --- Effizienz berechnen ---
    TH1* hEff = dynamic_cast<TH1*>(hRec->Clone("hEfficiency"));
    hEff->Divide(hRec, hGen, 1.0, 1.0, "B");

    // --- Effizienz pro Bin ausgeben ---
    std::cout << "\n--- Effizienz pro Bin ---\n";
    for (int bin = 1; bin <= nbins_rec; ++bin) {
        double pt_low  = hEff->GetXaxis()->GetBinLowEdge(bin);
        double pt_high = hEff->GetXaxis()->GetBinUpEdge(bin);
        double eff     = hEff->GetBinContent(bin);
        double err     = hEff->GetBinError(bin);

        std::cout << Form("pT: %.2f - %.2f | Effizienz: %.4f ± %.4f",
                          pt_low, pt_high, eff, err)
                  << std::endl;
    }

    // Effizienz-Plot
    Plotting1D PlotEff;
    PlotEff.New(hEff, "Efficiency", 20, 1.5, kBlack, "ep");
    PlotEff.SetAxisLabel("#bf{#it{p}_{T} (GeV/#it{c})}", "#bf{Efficiency}");
    PlotEff.SetMargins(0.12, 0.15, 0.05, 0.025);
    PlotEff.Plot("Figures/omegagen/" + output_name + "/efficiency.png");

    // In ROOT-Datei speichern
    TFile outFile("outputs/efficiencies/" + output_name + "_efficiency.root", "RECREATE");
    hEff->Write();
    outFile.Close();

    for (auto* file : input_files)
        if (file) { file->Close(); delete file; }

    std::cout << "byebye" << std::endl;
    std::cout.rdbuf(originalCoutBuffer);
}
