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


// Function to get invariant cross sections
void Inv_Crosssec(TH1* hist, TFile* outputfile, TH1* eff_hist, int Nevt, double deltay, double sig_TVX, double eff_TVX){
    int nbins = hist->GetNbinsX();

    TH1D* hCrosssec = (TH1D*) hist->Clone("hCrosssec");
    hCrosssec->Reset(); // nur die Achsen behalten

    for (int i = 1; i <= nbins; i++) {
        double pT = hist->GetBinCenter(i);
        double deltapT = hist->GetBinWidth(i);
        double nomega = hist->GetBinContent(i);
        double eff = eff_hist->GetBinContent(i);

        double sigmaN   = hist->GetBinError(i);
        double sigmaEff = eff_hist->GetBinError(i);

        double yield = nomega / (Nevt * eff * 2 * M_PI * deltay * pT * deltapT) * eff_TVX; 
        double crossec = nomega / (Nevt * eff * 2 * M_PI * deltay * pT * deltapT) * sig_TVX;
        
        // Fehlerfortpflanzung
        double error = crossec * sqrt( pow(sigmaN/nomega, 2) + pow(sigmaEff/eff, 2) );

        std::cout<<"pT: "<<pT<<", nomega: "<<nomega<<", eff: "<<eff<<", crossec: "<<crossec<<", yield: "<<yield<<std::endl;

        hCrosssec->SetBinContent(i, crossec);
        hCrosssec->SetBinError(i, error);
        if (eff <= 0 || nomega <= 0) {
            std::cout << "Skipping bin " << i << " due to zero content or efficiency\n";
            continue;
        }


    }
    outputfile->cd();
    hCrosssec->Write();
}




double deltay_OO = 0.8;   
double deltay_pp = 0.8;  
double eff_TVX_pp = 0.72;
double eff_TVX_OO = 0.935;           
double sig_TVX_pp = 0.0503; //b
double sig_TVX_OO = 1.13;   //b


// Function to fully correct omega spectra
void omegacorr(const TString& datarec_file, 
               const TString& eff_file, 
               const TString& output_name,
                const int Nevt_pp,
                const int Nevt_OO){

    std::ofstream logfile("logfiles/logfile_corr.txt");
    std::streambuf* originalCoutBuffer = std::cout.rdbuf();
    std::cout.rdbuf(logfile.rdbuf());

    std::cout << "Opening input files..." << std::endl;

    // open input files
    TFile *input = TFile::Open(datarec_file, "READ"); 
    TFile *input_eff = TFile::Open(eff_file, "READ");
    if(!input || input->IsZombie() || !input_eff || input_eff->IsZombie()){
        std::cerr << "Error opening input files!" << std::endl;
        return;
    }

    TFile* output = new TFile(output_name, "RECREATE");

    // get histograms (reconstructed data + MC efficiencies)
    TH1* hdatarec = dynamic_cast<TH1*>(input->Get("h_omegacounts"));
    TH1* heff = dynamic_cast<TH1*>(input_eff->Get("hEfficiency"));
    if(!hdatarec || !heff){
        std::cerr << "Error: could not retrieve histograms!" << std::endl;
        return;
    }

    // clone so original isn’t modified
    TH1* h_omegacounts_uncorr = dynamic_cast<TH1*>(hdatarec->Clone("h_omegacounts_uncorr"));

    if (datarec_file.Contains("pp")) {
        std::cout<<datarec_file<<std::endl;
        std::cout<<"pp"<<std::endl;
        Inv_Crosssec(h_omegacounts_uncorr, output, heff, Nevt_pp, deltay_pp, sig_TVX_pp, eff_TVX_pp);
    }
    else if (datarec_file.Contains("OO")) {
        std::cout<<datarec_file<<std::endl;
        std::cout<<"OO"<<std::endl;
        Inv_Crosssec(h_omegacounts_uncorr, output, heff, Nevt_OO, deltay_OO, sig_TVX_OO, eff_TVX_OO);
    } 
    else{std::cout<<"filename does not contain pp or OO"<<std::endl;}
    

    output->Close();
    input->Close();
    input_eff->Close();

    std::cout.rdbuf(originalCoutBuffer);
}
