#include <iostream>
#include <fstream>
#include <cmath>
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>


// Function to fit background
Double_t FunctionBGExclusionPol3(Double_t *x, Double_t *par){   
  double Momega_EMCal = 0.7822;
  double sigmaomega_EMCal = 0.014;                         

  double low_cut  = Momega_EMCal - 3 * sigmaomega_EMCal;
  double high_cut = Momega_EMCal + 3 * sigmaomega_EMCal;

  if (x[0] > low_cut && x[0] < high_cut){
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}

Double_t FunctionBGPol3(Double_t *x, Double_t *par){   
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}

Double_t LinFunc(Double_t *x, Double_t *par){
  return par[0]*x[0]+par[1];
}

Double_t FWHM_Func(Double_t sigma){
  return 2.0 * sqrt(2.0 * log(2.0)) * sigma;
}


//-----------------------------------------------------------------------------------------------------------------
//definitions
const int nGG = 3;
const int nEvent = 8;
int nTCTracks = 0;
TObject* obj = nullptr;
TString file[4] = {"MC/AnalysisResults_MC.root",
                    "MC/AnalysisResults_fixed_ppMC.root",
                    "MC/AnalysisResults_fixed_OOMC.root",
                    "MC/AnalysisResults_pp_MC_LargeStatistics.root"};
double Momega_EMCal = 0.7822; //in GeV/c2
double sigmaomega_EMCal = 0.014; //in GeV/c2
double low_cut  = (Momega_EMCal - 3 * sigmaomega_EMCal);  // 0.7402 in GeV/c
double high_cut = (Momega_EMCal + 3 * sigmaomega_EMCal);  // 0.8242 in GeV/c
double pT_bins[13] = {1.8, 2.2, 2.6, 3.2, 4.0,5.0,6.0,7.0,8.0,10.0,12.0,16.0,20.0};
int npT = sizeof(pT_bins) / sizeof(pT_bins[0]); // will be 13
double Momega_PDG = 0.78265; //in GeV/c2 /pm 0.16
double sigmaomega_PDG = 0.00849; //in GeV/c2 /pm 0.00008
int n_pTstart = -999;
double xlow = 0.6;
double xup = 1.;
double sigma_BG = 0.;
std::vector<double> index_list;

double hline_y2[12]={0.5*1e3,0.8*1e3,1.5*1e3,1.8*1e3,1.8*1e3,1.*1e3,0.6*1e3,0.5*1e3,0.5*1e3,0.3*1e3,0.2*1e3,50};

//Root Folder structure in detail
TString dir_names_ppMC[8] = {"track-propagation","lumi-task","event-selection-task","bc-selection-task",
  "mc-generator-studies","emcal-correction-task","skimmer-gamma-calo","omega-meson-emc"};

TString dir_names_ppMC_fixed[11] = {"track-propagation","skimmer-gamma-calo",
  "mc-generator-studies_MinBias","mc-generator-studies","emcal-correction-task",
  "lumi-task","event-selection-task","bc-selection-task","omega-meson-emc",
  "omega-meson-emc_L0Triggered","omega-meson-emc_MinBias"};

TString dir_names_OOMC_fixed[11] = {"track-propagation","lumi-task",
  "event-selection-task","bc-selection-task","mc-generator-studies_MinBias",
  "mc-generator-studies","emcal-correction-task","skimmer-gamma-calo",
  "omega-meson-emc","omega-meson-emc_L0Triggered","omega-meson-emc_MinBias"};

TString dir_names_ppMC_large[11] = {"track-propagation","lumi-task",
  "event-selection-task","bc-selection-task","mc-generator-studies_MinBias",
  "mc-generator-studies","emcal-correction-task","skimmer-gamma-calo",
  "omega-meson-emc","omega-meson-emc_L0Triggered","omega-meson-emc_MinBias"};  


//MC gen
TString OM_ppmc_gen_hists[16] = {"hCollisionCounter","hBCCounter","Yield","Yield_Accepted",
                                  "Yield_T","Yield_TZ","Yield_TZS","Yield_TZSG","Yield_TZSGU","Yield_TZSGUE",
                                  "Yield_TZSGUE_Accepted","Yield_BC_T","Yield_BC_TC","NCollisionsMCCollisions",
                                  "NTVXCollisionsMCCollisions","hEMCollisionCounter"};

//MC pp gen fixed !Min Bias!
TString OM_ppmc_fixed_gen_hists[16] = {"hCollisionCounter","hBCCounter","Yield","Yield_Accepted",
                                  "Yield_T","Yield_TZ","Yield_TZS","Yield_TZSG","Yield_TZSGU","Yield_TZSGUE",
                                  "Yield_TZSGUE_Accepted","Yield_BC_T","Yield_BC_TC","NCollisionsMCCollisions",
                                  "NTVXCollisionsMCCollisions","hEMCollisionCounter"};

//MC OO gen fixed !Min Bias!
TString OM_OOmc_fixed_gen_hists[16] = {"hCollisionCounter","hBCCounter","Yield","Yield_Accepted",
                                  "Yield_T","Yield_TZ","Yield_TZS","Yield_TZSG","Yield_TZSGU","Yield_TZSGUE",
                                  "Yield_TZSGUE_Accepted","Yield_BC_T","Yield_BC_TC","NCollisionsMCCollisions",
                                  "NTVXCollisionsMCCollisions","hEMCollisionCounter"};

//MC pp gen large !Min Bias!
TString OM_ppmc_large_gen_hists[16] = {"hCollisionCounter","hBCCounter","Yield","Yield_Accepted",
                                  "Yield_T","Yield_TZ","Yield_TZS","Yield_TZSG","Yield_TZSGU","Yield_TZSGUE",
                                  "Yield_TZSGUE_Accepted","Yield_BC_T","Yield_BC_TC","NCollisionsMCCollisions",
                                  "NTVXCollisionsMCCollisions","hEMCollisionCounter"};


//outputnames
TString output_names[4] = {"outputs/omegagen/ppMC_gen_hists.root",
                          "outputs/omegagen/ppMC_fixed_gen_hists.root",
                          "outputs/omegagen/OOMC_fixed_gen_hists.root",
                          "outputs/omegagen/ppMC_large_gen_hists.root"};

int output_length = sizeof(output_names) / sizeof(output_names[0]);



//Function to get generated MC from ROOT files
void omegagen(){
  std::ofstream logfile("logfiles/logfile_gen.txt");
  std::streambuf* originalCoutBuffer = std::cout.rdbuf();
  std::cout.rdbuf(logfile.rdbuf());

  for(int j=3;j<output_length;++j){
    TFile *input = TFile::Open(file[j],"read");
    TFile* output = new TFile(output_names[j], "RECREATE");
    
    if(j==0){
      TDirectory* dir_OM_ppmc_gen = (TDirectory*)input->Get(dir_names_ppMC[4]);
      //MC gen
      obj = dir_OM_ppmc_gen->Get(OM_ppmc_gen_hists[9]);
      if(!obj){std::cout<<"obj not find"<<std::endl;continue;}
      TH1* h1 = dynamic_cast<TH1*>(obj);
      h1->Sumw2();
      h1->Write();
    }

    else if(j==1){
      TDirectory* dir_OM_ppmc_fixed_gen = (TDirectory*)input->Get(dir_names_ppMC_fixed[2]);
      //MC gen !Min Bias!
      obj = dir_OM_ppmc_fixed_gen->Get(OM_ppmc_fixed_gen_hists[9]);
      if(!obj){std::cout<<"obj not found"<<std::endl;continue;}
      TH1* h1 = dynamic_cast<TH1*>(obj);
      h1->Sumw2();
      h1->Write();
      continue;
    }

    else if(j==2){
      TDirectory* dir_OM_OOmc_fixed_gen = (TDirectory*)input->Get(dir_names_OOMC_fixed[4]);
      //MC gen !Min Bias!
      obj = dir_OM_OOmc_fixed_gen->Get(OM_OOmc_fixed_gen_hists[9]);
      if(!obj){std::cout<<"obj not found"<<std::endl;continue;}
      TH1* h1 = dynamic_cast<TH1*>(obj);
      h1->Sumw2();
      int nbins = h1->GetNbinsX();
      for(int i=1; i<=nbins; ++i){
          double content = h1->GetBinContent(i);
          double error = h1->GetBinError(i);
      }

      h1->Write();
      continue;
    }
    
    else if(j==3){
      TDirectory* dir_OM_ppmc_large_gen = (TDirectory*)input->Get(dir_names_ppMC_large[4]);
      //MC gen !Min Bias!
      obj = dir_OM_ppmc_large_gen->Get(OM_ppmc_large_gen_hists[9]);
      if(!obj){std::cout<<"obj not found"<<std::endl;continue;}
      TH1* h1 = dynamic_cast<TH1*>(obj);
      h1->Sumw2();
      h1->Write();
      continue;
    }

    input->Close();
    output->Close();
  }
  std::cout.rdbuf(originalCoutBuffer);
} 
