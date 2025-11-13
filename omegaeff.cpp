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


//Function to correct reconstructed data with MC efficiencies
void omegaeff(const TString& datarec_file, 
               const TString& eff_file, 
               const TString& output_name){

    std::ofstream logfile("logfiles/logfile_eff.txt");
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

    // get histograms
    TH1* hdatarec = dynamic_cast<TH1*>(input->Get("h_omegacounts"));
    TH1* heff = dynamic_cast<TH1*>(input_eff->Get("hEfficiency"));
    if(!hdatarec || !heff){
        std::cerr << "Error: could not retrieve histograms!" << std::endl;
        return;
    }

    // clone so original isn’t modified
    TH1* hcorr = dynamic_cast<TH1*>(hdatarec->Clone("h_omegacounts_corrected"));

    // multiply with efficiency
    hcorr->Multiply(heff);

    // save to output
    TFile* output = new TFile(output_name, "RECREATE");
    hcorr->Write();
    output->Close();

    std::cout << "Corrected histogram written to " << output_name << std::endl;

    // cleanup
    input->Close();
    input_eff->Close();
    std::cout.rdbuf(originalCoutBuffer);
}
