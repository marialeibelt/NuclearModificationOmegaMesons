#include <iostream>
#include <fstream>
#include <cmath>
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>



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
TString file[15] = {"Data/AnalysisResults.root",
                  "Data/AnalysisResults_OldTask_SmallppRefDataSet_Run3Alignment.root",
                  "Data/AnalysisResults_NewTask_SmallppRefDataSet_Run3Alignment.root",  
                  "Data/AnalysisResults_NewTask_MediumppRefDataSet_Run2Alignment.root",
                  "Data/AnalysisResults_NewTask_MediumppRefDataSet_Run3Alignment.root",
                  "Data/AnalysisResults_NewTask_SmallOODataSet_Run3Alignment.root",     
                  "Data/AnalysisResults_ppRefNotAligned.root",
                  "MC/AnalysisResults_MC.root",
                  "MC/AnalysisResults_fixed_ppMC.root",                                 
                  "MC/AnalysisResults_fixed_OOMC.root",                                 //9
                  "MC/AnalysisResults_pp_MC_LargeStatistics.root",                      //10
                  "Data/AnalysisResults_OOData_small_new.root",
                  "Data/AnalysisResults-ppDataQMedium.root",
                  "Data/AnalysisResults_pp_full.root",                                  //13
                  "Data/AnalysisResults-OOData.root"};

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
std::vector<double> index_list;
std::vector<double> chi2ndf_list;
std::vector<double> index_list_signalfit;
std::vector<double> chi2ndf_list_signalfit;
std::vector<double> mean_list_signalfit;
std::vector<double> width_list_signalfit;
std::vector<double> mean_error_list_signalfit;
std::vector<double> width_error_list_signalfit;
std::vector<double> omega_counts_list_signalfit;
std::vector<double> signal_integral_vec;
std::vector<double> signal_error_vec;

double hline_y2[12]={0.5*1e3,0.8*1e3,1.5*1e3,1.8*1e3,1.8*1e3,1.*1e3,0.6*1e3,0.5*1e3,0.5*1e3,0.3*1e3,0.2*1e3,50};
TObject* obj;

//Name lists
TString dir_names[11] = {"track-propagation","tof-signal","bc-selection-task",
  "lumi-task","event-selection-task","emcal-correction-task","ft0-corrected-table",
  "tof-pid-merge","photon-concersion-builder","skimmer-gamma-calo","heavy-neutral-meson"};

TString dir_names_OldTask_Smallpp[11] = {"track-propagation","emcal-correction-task",
  "lumi-task","event-selection-task","tof-signal","bc-selection-task",
  "ft0-corrected-table","tof-pid-merge","skimmer-gamma-calo",
  "photon-concersion-builder","heavy-neutral-meson-filter"};

TString dir_names_NewTask_Smallpp[13] = {"track-propagation","tof-signal","bc-selection-task",
  "emcal-correction-task","lumi-task","event-selection-task","ft0-corrected-table",
  "tof-pid-merge","skimmer-gamma-calo","strangeness-builder","omega-meson-emc",
  "tree-creator-ele-ml-dda","ml-track-qc"};

TString dir_names_NewTask_Mediumpp_Run2Alignment[7] = {"track-propagation","skimmer-gamma-calo","emcal-correction-task",
  "lumi-task","event-selection-task","bc-selection-task","omega-meson-emc"};

TString dir_names_NewTask_Mediumpp_Run3Alignment[7] = {"track-propagation","bc-selection-task","lumi-task","event-selection-task",
  "emcal-correction-task","skimmer-gamma-calo","omega-meson-emc"};

TString dir_names_NewTask_Smallpp_OO_Run3Alignment[7] = {"track-propagation","skimmer-gamma-calo","emcal-correction-task",
  "lumi-task","event-selection-task","bc-selection-task","omega-meson-emc"};

TString dir_names_pp_NoAlignment[7] = {"track-propagation","omega-meson-emc","skimmer-gamma-calo","emcal-correction-task",
  "lumi-task","event-selection-task","bc-selection-task"};

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

TString dir_names_OOData_small[9] = {"track-propagation","omega-meson-emc",
  "skimmer-gamma-calo","emcal-correction-task","lumi-task",
  "event-selection-task","bc-selection-task",
  "omega-meson-emc_L0Triggered","omega-meson-emc_MinBias"};  

TString dir_names_ppData[9] = {"track-propagation","lumi-task","event-selection-task",
  "bc-selection-task","emcal-correction-task","skimmer-gamma-calo",
  "omega-meson-emc","omega-meson-emc_L0Triggered","omega-meson-emc_MinBias"};  

TString dir_names_ppData_full[9] = {"track-propagation","omega-meson-emc","skimmer-gamma-calo",
  "emcal-correction-task","lumi-task","event-selection-task",
  "bc-selection-task","omega-meson-emc_L0Triggered","omega-meson-emc_MinBias"};

TString dir_names_OOData[9] = {"track-propagation","omega-meson-emc","skimmer-gamma-calo",
  "emcal-correction-task","lumi-task","event-selection-task",
  "bc-selection-task","omega-meson-emc_L0Triggered","omega-meson-emc_MinBias"};


//OLD original pp
TString HNM_dirs[4] = {"Event;1","GG;1","TrackCuts;1","HNM;1"};
TString HNM_folders[2] = {"Before;1","After;1"};
TString HNM_After_folders[2] = {"Omega;1","EtaPrime;1"};
TString HNM_After_Omega_folders[2] = {"PCM","EMC"};
TString HNM_After_Omega_PCM_hists[3] = {"fInvMassVsPt;1","fEta;1","fPhi;1"};
TString HNM_After_Omega_EMC_hists[3] = {"fInvMassVsPt;1","fEta;1","fPhi;1"};
TString HNM_After_EtaPrime_folders[2] = {"PCM;1","EMC;1"};
TString HNM_After_EtaPrime_PCM_hists[3] = {"fInvMassVsPt;1","fEta;1","fPhi;1"};
TString HNM_After_EtaPrime_EMC_hists[3] = {"fInvMassVsPt;1","fEta;1","fPhi;1"};

//OLD pp
TString OM_OLD_dirs[10] = {"Event","GG","TrackCuts","HNM","omegap","etaprimep","omegad","etaprimed","ppomega","ppetaprime"};
TString OM_OLD_folders[2] = {"Before","After"};
TString OM_OLD_After_folders[2] = {"Omega","EtaPrime"};
TString OM_OLD_After_Omega_folders[2] = {"PCM","EMC"};
TString OM_OLD_After_Omega_PCM_hists[3] = {"fInvMassVsPt;1","fEta;1","fPhi;1"};
TString OM_OLD_After_Omega_EMC_hists[3] = {"fInvMassVsPt;1","fEta;1","fPhi;1"};
TString OM_OLD_After_EtaPrime_folders[2] = {"PCM;1","EMC;1"};
TString OM_OLD_After_EtaPrime_PCM_hists[3] = {"fInvMassVsPt;1","fEta;1","fPhi;1"};
TString OM_OLD_After_EtaPrime_EMC_hists[3] = {"fInvMassVsPt;1","fEta;1","fPhi;1"};

//NEW small pp
TString OM_Newsmallpp_dirs[4] = {"Event","GG","TrackCuts","Omega"};
TString OM_Newsmallpp_folders[2] = {"Before","After"};
TString OM_Newsmallpp_Before_folders[1] = {"PiPlPiMi"};
TString OM_Newsmallpp_Before_PiPlPiMi_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_Newsmallpp_Before_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_Newsmallpp_After_hists[3] = {"fInvMassVsPt","fEta","fPhi"};

//NEW Run2Alignment pp
TString OM_Newmediumpp_run2al_dirs[4] = {"Event","GG","TrackCuts","Omega"};
TString OM_Newmediumpp_run2al_folders[2] = {"Before","After"};
TString OM_Newmediumpp_run2al_Before_folders[1] = {"PiPlPiMi"};
TString OM_Newmediumpp_run2al_Before_PiPlPiMi_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_Newmediumpp_run2al_Before_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_Newmediumpp_run2al_After_hists[3] = {"fInvMassVsPt","fEta","fPhi"};

//NEW Run3Alignment pp
TString OM_Newmediumpp_run3al_dirs[4] = {"Event","GG","TrackCuts","Omega"};
TString OM_Newmediumpp_run3al_folders[2] = {"Before","After"};
TString OM_Newmediumpp_run3al_Before_folders[1] = {"PiPlPiMi"};
TString OM_Newmediumpp_run3al_Before_PiPlPiMi_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_Newmediumpp_run3al_Before_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_Newmediumpp_run3al_After_hists[3] = {"fInvMassVsPt","fEta","fPhi"};

//Data OO small
TString OM_NewsmallOO_run3al_dirs[4] = {"Event","GG","TrackCuts","Omega"};
TString OM_NewsmallOO_run3al_folders[2] = {"Before","After"};
TString OM_NewsmallOO_run3al_Before_folders[1] = {"PiPlPiMi"};
TString OM_NewsmallOO_run3al_Before_PiPlPiMi_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_NewsmallOO_run3al_Before_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_NewsmallOO_run3al_After_hists[3] = {"fInvMassVsPt","fEta","fPhi"};

//No Alignment pp
TString OM_pp_noal_dirs[4] = {"Event","GG","TrackCuts","Omega"};
TString OM_pp_noal_folders[2] = {"Before","After"};
TString OM_pp_noal_Before_folders[1] = {"PiPlPiMi"};
TString OM_pp_noal_Before_PiPlPiMi_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_pp_noal_Before_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_pp_noal_After_hists[3] = {"fInvMassVsPt","fEta","fPhi"};

//MC pp
TString OM_ppmc_dirs[4] = {"Event","GG","TrackCuts","Omega"};
TString OM_ppmc_folders[2] = {"Before","After"};
TString OM_ppmc_Before_folders[1] = {"PiPlPiMi"};
TString OM_ppmc_Before_PiPlPiMi_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_ppmc_Before_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_ppmc_After_hists[3] = {"fInvMassVsPt","fEta","fPhi"};

//MC pp fixed
TString OM_ppmc_fixed_dirs[4] = {"Event","GG","TrackCuts","Omega"};
TString OM_ppmc_fixed_folders[2] = {"Before","After"};
TString OM_ppmc_fixed_Before_folders[1] = {"PiPlPiMi"};
TString OM_ppmc_fixed_Before_PiPlPiMi_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_ppmc_fixed_Before_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_ppmc_fixed_After_hists[3] = {"fInvMassVsPt","fEta","fPhi"};

//MC OO fixed
TString OM_OOmc_fixed_dirs[4] = {"Event","GG","TrackCuts","Omega"};
TString OM_OOmc_fixed_folders[2] = {"Before","After"};
TString OM_OOmc_fixed_Before_folders[1] = {"PiPlPiMi"};
TString OM_OOmc_fixed_Before_PiPlPiMi_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_OOmc_fixed_Before_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_OOmc_fixed_After_hists[3] = {"fInvMassVsPt","fEta","fPhi"};

//MC pp large
TString OM_ppmc_large_dirs[4] = {"Event","GG","TrackCuts","Omega"};
TString OM_ppmc_large_folders[2] = {"Before","After"};
TString OM_ppmc_large_Before_folders[1] = {"PiPlPiMi"};
TString OM_ppmc_large_Before_PiPlPiMi_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_ppmc_large_Before_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_ppmc_large_After_hists[3] = {"fInvMassVsPt","fEta","fPhi"};

//Data OO small
TString OM_OOData_small_dirs[4] = {"Event","GG","TrackCuts","Omega"};
TString OM_OOData_small_folders[2] = {"Before","After"};
TString OM_OOData_small_Before_folders[1] = {"PiPlPiMi"};
TString OM_OOData_small_Before_PiPlPiMi_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_OOData_small_Before_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_OOData_small_After_hists[3] = {"fInvMassVsPt","fEta","fPhi"};

//Data pp
TString OM_ppData_dirs[4] = {"Event","GG","TrackCuts","Omega"};
TString OM_ppData_folders[2] = {"Before","After"};
TString OM_ppData_Before_folders[1] = {"PiPlPiMi"};
TString OM_ppData_Before_PiPlPiMi_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_ppData_Before_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_ppData_After_hists[3] = {"fInvMassVsPt","fEta","fPhi"};

//Data pp full
TString OM_ppData_full_dirs[4] = {"Event","GG","TrackCuts","Omega"};
TString OM_ppData_full_folders[2] = {"Before","After"};
TString OM_ppData_full_Before_folders[1] = {"PiPlPiMi"};
TString OM_ppData_full_Before_PiPlPiMi_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_ppData_full_Before_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_ppData_full_After_hists[3] = {"fInvMassVsPt","fEta","fPhi"};

//Data OO
TString OM_OOData_dirs[4] = {"Event","GG","TrackCuts","Omega"};
TString OM_OOData_folders[2] = {"Before","After"};
TString OM_OOData_Before_folders[1] = {"PiPlPiMi"};
TString OM_OOData_Before_PiPlPiMi_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_OOData_Before_hists[3] = {"fInvMassVsPt","fEta","fPhi"};
TString OM_OOData_After_hists[3] = {"fInvMassVsPt","fEta","fPhi"};


//outputnames
TString output_names[15] = {"outputs/omegarec/OldOriginalTask_SmallppRefDataSet_hists.root",
  "outputs/omegarec/OldTask_SmallppRefDataSet_hists.root",
  "outputs/omegarec/NewTask_SmallppRefDataSet_hists.root",
  "outputs/omegarec/NewTask_MediumppRefDataSet_Run2Alignment_hists.root",
  "outputs/omegarec/NewTask_MediumppRefDataSet_Run3Alignment_hists.root",
  "outputs/omegarec/NewTask_SmallOODataSet_Run3Alignment_hists.root",
  "outputs/omegarec/ppRefNotAligned_hists.root",
  "outputs/omegarec/ppMC_hists.root",
  "outputs/omegarec/ppMC_fixed_hists.root",
  "outputs/omegarec/OOMC_fixed_hists.root",
  "outputs/omegarec/ppMC_large_hists.root",
  "outputs/omegarec/OOData_small_hists.root",
  "outputs/omegarec/ppData_hists.root",
  "outputs/omegarec/ppData_full_hists.root",
  "outputs/omegarec/OOData_hists.root"};

int output_length = sizeof(output_names) / sizeof(output_names[0]);



void omegarec(){
  std::ofstream logfile("logfiles/logfile_rec.txt");
  std::streambuf* originalCoutBuffer = std::cout.rdbuf();
  std::cout.rdbuf(logfile.rdbuf());

  for(int j=0;j<output_length;++j){ //j=0
    std::cout<<"j: "<<j<<std::endl;

    chi2ndf_list.clear();
    chi2ndf_list_signalfit.clear();
    mean_list_signalfit.clear();
    width_list_signalfit.clear();
    mean_error_list_signalfit.clear();
    width_error_list_signalfit.clear();
    omega_counts_list_signalfit.clear();
    signal_integral_vec.clear();
    signal_error_vec.clear();

    std::cout<<"output_names: "<<output_names[j]<<std::endl;
    TFile *input = TFile::Open(file[j],"read");
    TFile* output = new TFile(output_names[j], "RECREATE");

    if(j==0){
      TDirectory* dir_HNM = (TDirectory*)input->Get(dir_names[10]);
      TDirectory* dir_HNM1 = (TDirectory*)dir_HNM->Get(HNM_dirs[3]);
      TDirectory* dir_HNM_After = (TDirectory*)dir_HNM1->Get(HNM_folders[1]);
      TDirectory* dir_HNM_After_Omega = (TDirectory*)dir_HNM_After->Get(HNM_After_folders[0]);
      TDirectory* dir_HNM_After_Omega_EMC = (TDirectory*)dir_HNM_After_Omega->Get(HNM_After_Omega_folders[1]);
      obj = dir_HNM_After_Omega_EMC->Get(HNM_After_Omega_EMC_hists[0]);
    }

    else if(j==1){
      TDirectory* dir_OM_Old = (TDirectory*)input->Get(dir_names_OldTask_Smallpp[10]);
      TDirectory* dir_OM1_Old = (TDirectory*)dir_OM_Old->Get(OM_OLD_dirs[3]);
      TDirectory* dir_OM_OLD_After = (TDirectory*)dir_OM1_Old->Get(OM_OLD_folders[1]);
      TDirectory* dir_OM_OLD_After_Omega = (TDirectory*)dir_OM_OLD_After->Get(OM_OLD_After_folders[0]);
      TDirectory* dir_OM_OLD_After_Omega_EMC = (TDirectory*)dir_OM_OLD_After_Omega->Get(OM_OLD_After_Omega_folders[1]);
      obj = dir_OM_OLD_After_Omega_EMC->Get(OM_OLD_After_Omega_EMC_hists[0]);
    }

    else if(j==2){
      TDirectory* dir_OM_Newsmallpp = (TDirectory*)input->Get(dir_names_NewTask_Smallpp[10]);
      TDirectory* dir_OM1_Newsmallpp = (TDirectory*)dir_OM_Newsmallpp->Get(OM_Newsmallpp_dirs[3]);
      TDirectory* dir_OM_Newsmallpp_After = (TDirectory*)dir_OM1_Newsmallpp->Get(OM_Newsmallpp_folders[1]);
      obj = dir_OM_Newsmallpp_After->Get(OM_Newsmallpp_After_hists[0]);
    }

    else if(j==3){
      TDirectory* dir_OM_Newmediumpp_run2al = (TDirectory*)input->Get(dir_names_NewTask_Mediumpp_Run2Alignment[6]);
      TDirectory* dir_OM1_Newmediumpp_run2al = (TDirectory*)dir_OM_Newmediumpp_run2al->Get(OM_Newmediumpp_run2al_dirs[3]);
      TDirectory* dir_OM_Newmediumpp_run2al_After = (TDirectory*)dir_OM1_Newmediumpp_run2al->Get(OM_Newmediumpp_run2al_folders[1]);
      obj = dir_OM_Newmediumpp_run2al_After->Get(OM_Newmediumpp_run2al_After_hists[0]);
    }

    else if(j==4){
      TDirectory* dir_OM_Newmediumpp_run3al = (TDirectory*)input->Get(dir_names_NewTask_Mediumpp_Run3Alignment[6]);
      TDirectory* dir_OM1_Newmediumpp_run3al = (TDirectory*)dir_OM_Newmediumpp_run3al->Get(OM_Newmediumpp_run3al_dirs[3]);
      TDirectory* dir_OM_Newmediumpp_run3al_After = (TDirectory*)dir_OM1_Newmediumpp_run3al->Get(OM_Newmediumpp_run3al_folders[1]);
      obj = dir_OM_Newmediumpp_run3al_After->Get(OM_Newmediumpp_run3al_After_hists[0]);
    }

    else if(j==5){
      TDirectory* dir_OM_NewsmallOO_run3al = (TDirectory*)input->Get(dir_names_NewTask_Smallpp_OO_Run3Alignment[6]);
      TDirectory* dir_OM1_NewsmallOO_run3al = (TDirectory*)dir_OM_NewsmallOO_run3al->Get(OM_NewsmallOO_run3al_dirs[3]);
      TDirectory* dir_OM_NewsmallOO_run3al_After = (TDirectory*)dir_OM1_NewsmallOO_run3al->Get(OM_NewsmallOO_run3al_folders[1]);
      obj = dir_OM_NewsmallOO_run3al_After->Get(OM_NewsmallOO_run3al_After_hists[0]);
    }

    else if(j==6){
      TDirectory* dir_OM_pp_noal = (TDirectory*)input->Get(dir_names_pp_NoAlignment[1]);
      TDirectory* dir_OM1_pp_noal  = (TDirectory*)dir_OM_pp_noal->Get(OM_pp_noal_dirs[3]);
      TDirectory* dir_OM_pp_noal_After = (TDirectory*)dir_OM1_pp_noal->Get(OM_pp_noal_folders[1]);
      obj = dir_OM_pp_noal_After->Get(OM_pp_noal_After_hists[0]);
    }

    else if(j==7){
      TDirectory* dir_OM_ppmc = (TDirectory*)input->Get(dir_names_ppMC[7]);
      TDirectory* dir_OM1_ppmc = (TDirectory*)dir_OM_ppmc->Get(OM_ppmc_dirs[3]);
      TDirectory* dir_OM_ppmc_After = (TDirectory*)dir_OM1_ppmc->Get(OM_ppmc_folders[1]);
      obj = dir_OM_ppmc_After->Get(OM_ppmc_After_hists[0]);
    }

    else if(j==8){
      //Min Bias!
      TDirectory* dir_OM_ppmc_fixed = (TDirectory*)input->Get(dir_names_ppMC_fixed[10]);
      TDirectory* dir_OM1_ppmc_fixed = (TDirectory*)dir_OM_ppmc_fixed->Get(OM_ppmc_fixed_dirs[3]);
      TDirectory* dir_OM_ppmc_fixed_After = (TDirectory*)dir_OM1_ppmc_fixed->Get(OM_ppmc_fixed_folders[1]);
      obj = dir_OM_ppmc_fixed_After->Get(OM_ppmc_fixed_After_hists[0]);
    }

    else if(j==9){
      //Min Bias! OO MC
      TDirectory* dir_OM_OOmc_fixed = (TDirectory*)input->Get(dir_names_OOMC_fixed[10]);
      TDirectory* dir_OM1_OOmc_fixed = (TDirectory*)dir_OM_OOmc_fixed->Get(OM_OOmc_fixed_dirs[3]);
      TDirectory* dir_OM_OOmc_fixed_After = (TDirectory*)dir_OM1_OOmc_fixed->Get(OM_OOmc_fixed_folders[1]);
      obj = dir_OM_OOmc_fixed_After->Get(OM_OOmc_fixed_After_hists[0]);
    }

    else if(j==10){
      //Min Bias! pp large MC 
      TDirectory* dir_OM_ppmc_large = (TDirectory*)input->Get(dir_names_ppMC_large[10]);
      TDirectory* dir_OM1_ppmc_large = (TDirectory*)dir_OM_ppmc_large->Get(OM_ppmc_large_dirs[3]);
      TDirectory* dir_OM_ppmc_large_After = (TDirectory*)dir_OM1_ppmc_large->Get(OM_ppmc_large_folders[1]);
      obj = dir_OM_ppmc_large_After->Get(OM_ppmc_large_After_hists[0]);
    }

    else if(j==11){
      //Min Bias! OO small Data
      TDirectory* dir_OM_OOData_small = (TDirectory*)input->Get(dir_names_OOData_small[8]);
      TDirectory* dir_OM1_OOData_small = (TDirectory*)dir_OM_OOData_small->Get(OM_OOData_small_dirs[3]);
      TDirectory* dir_OM_OOData_small_After = (TDirectory*)dir_OM1_OOData_small->Get(OM_OOData_small_folders[1]);
      obj = dir_OM_OOData_small_After->Get(OM_OOData_small_After_hists[0]);
    }

    else if(j==12){
      //Min Bias! pp Data
      TDirectory* dir_OM_ppData = (TDirectory*)input->Get(dir_names_ppData[8]);
      TDirectory* dir_OM1_ppData = (TDirectory*)dir_OM_ppData->Get(OM_ppData_dirs[3]);
      TDirectory* dir_OM_ppData_After = (TDirectory*)dir_OM1_ppData->Get(OM_ppData_folders[1]);
      obj = dir_OM_ppData_After->Get(OM_ppData_After_hists[0]);
    }

    else if(j==13){
      //Min Bias! pp Data full dataset
      TDirectory* dir_OM_ppData_full = (TDirectory*)input->Get(dir_names_ppData_full[8]);
      TDirectory* dir_OM1_ppData_full = (TDirectory*)dir_OM_ppData_full->Get(OM_ppData_full_dirs[3]);
      TDirectory* dir_OM_ppData_full_After = (TDirectory*)dir_OM1_ppData_full->Get(OM_ppData_full_folders[1]);
      obj = dir_OM_ppData_full_After->Get(OM_ppData_full_After_hists[0]);
    }

    else if(j==14){
      //Min Bias! OO Data
      TDirectory* dir_OM_OOData = (TDirectory*)input->Get(dir_names_OOData[8]);
      TDirectory* dir_OM1_OOData = (TDirectory*)dir_OM_OOData->Get(OM_OOData_dirs[3]);
      TDirectory* dir_OM_OOData_After = (TDirectory*)dir_OM1_OOData->Get(OM_OOData_folders[1]);
      obj = dir_OM_OOData_After->Get(OM_OOData_After_hists[0]);
    }
    

    //plot fInvMassVsPt
    TH2F* h2 = dynamic_cast<TH2F*>(obj);
    if(!h2){std::cout<<"Can't find h2"<<std::endl;}
    h2->Sumw2(); 

    int n_bins = sizeof(pT_bins)/sizeof(double) - 1;

    int n_pTstart = 3;
    for (int k = n_pTstart; k < n_bins; ++k) { //k=1... leave out first bin range
      int bin_y_low = h2->GetYaxis()->FindBin(pT_bins[k]+1e-5); //excluded Grenzwert
      int bin_y_high = h2->GetYaxis()->FindBin(pT_bins[k+1]-1e-5);

      TH1D* proj = (TH1D*) h2->ProjectionX(Form("proj_pT_%d", k), bin_y_low, bin_y_high);         
      if (!proj) {std::cout<<"k:"<<k<<", Histogram proj not found!"<<std::endl;}           

      //rebin                                                                                           
      proj->Rebin(7);                                                                                      
      proj->GetXaxis()->SetRangeUser(xlow,xup);                                                            

      //Fit BG
      TF1* fitFunc = new TF1("fitFunc", FunctionBGExclusionPol3, xlow, xup, 4);
      fitFunc->SetParameters(1, 1, 1, 1);
      proj->Fit(fitFunc,"QN");                                                                              

      double chi2 = fitFunc->GetChisquare();
      int ndf = fitFunc->GetNDF();
      double chi2ndf = chi2 / ndf;
      chi2ndf_list.push_back(chi2ndf);

      //BG
      TF1* BG = new TF1("BG", FunctionBGPol3, xlow, xup, 4);
      BG->SetParameters(fitFunc->GetParameters());
      
      //Extract signal  
      proj->Sumw2();
      TH1D* h_signal = (TH1D*)proj->Clone("h_signal"); 
      h_signal->Sumw2();

      int nbins_proj = proj->GetNbinsX();
      for (int i = 1; i <= nbins_proj; ++i) {
        double sigma_BG = 0.;
        double x_proj = proj->GetBinCenter(i);
        double content_proj = proj->GetBinContent(i);
        double error_proj = proj->GetBinError(i);
        double bg = BG->Eval(x_proj);  // background estimate

        // bg error
        for (int p = 0; p < 4; ++p) {
          double par_error = BG->GetParError(p);
          double deriv = 1.0;
          if (p == 1) deriv = x_proj;
          else if (p == 2) deriv = x_proj * x_proj;
          else if (p == 3) deriv = x_proj * x_proj * x_proj;
          sigma_BG += (deriv * par_error) * (deriv * par_error);
        }
        sigma_BG = sqrt(sigma_BG);

        double content_signal = content_proj - bg;
        double error_signal = sqrt(error_proj * error_proj + sigma_BG * sigma_BG);
        h_signal->SetBinContent(i, content_signal);
        h_signal->SetBinError(i, error_signal);
      }

      // Fit signal with Gaussian
      TF1* gausFit = new TF1("gausFit", "gaus", low_cut, high_cut);
      gausFit->SetParameters(h_signal->GetMaximum(), Momega_EMCal, sigmaomega_EMCal);
      gausFit->SetParLimits(1, low_cut, high_cut);
      h_signal->Fit(gausFit, "Q");//change to R  if you want output of fitfeatures

      double chi2_signalfit = gausFit->GetChisquare();
      int ndf_signalfit = gausFit->GetNDF();
      double chi2ndf_signalfit = chi2_signalfit / ndf_signalfit;
      chi2ndf_list_signalfit.push_back(chi2ndf_signalfit);
      
      double mean_signalfit = gausFit->GetParameter(1);
      double width_signalfit = gausFit->GetParameter(2);
      double mean_error_signalfit  = gausFit->GetParError(1); // Error on mean
      double width_error_signalfit = gausFit->GetParError(2); // Error on sigma

      mean_list_signalfit.push_back(mean_signalfit);
      width_list_signalfit.push_back(width_signalfit);
      mean_error_list_signalfit.push_back(mean_error_signalfit);
      width_error_list_signalfit.push_back(width_error_signalfit);

     // Integrate signal fit
      double omega_counts_fit = gausFit->Integral(low_cut, high_cut); // keine Division durch BinWidth
      omega_counts_list_signalfit.push_back(omega_counts_fit);

      // Integrate signal (absolute counts)
      double signal_integral = 0.;
      double signal_error2 = 0.;
      //std::cout << "\n bin k=" << k << std::endl;

      for (int i = 1; i <= h_signal->GetNbinsX(); ++i) {
          double bin_center = h_signal->GetBinCenter(i);

          if (bin_center >= low_cut && bin_center <= high_cut) {
              double bin_content = h_signal->GetBinContent(i);
              double bin_error   = h_signal->GetBinError(i);

              //std::cout << "bin_content: " << bin_content << ", bin_error: " << bin_error << std::endl;

              signal_integral += bin_content;                  // absolute Counts
              signal_error2 += bin_error * bin_error;         // Fehler quadratisch summieren

              //std::cout << "signal_integral: " << signal_integral 
                //        << ", signal_error2: " << signal_error2 << std::endl;
          }
      }

      signal_integral_vec.push_back(signal_integral);
      signal_error_vec.push_back(sqrt(signal_error2));  // Fehler als sqrt(sum of squares)
      std::cout << "signal_integral: " << signal_integral<< ", signal_error: " << sqrt(signal_error2) << std::endl;

      TF1* bg_clone = (TF1*) BG->Clone(Form("BG_%d", k));
      TH1D* signal_clone = (TH1D*) h_signal->Clone(Form("signal_clone_%d", k));
      TF1* gaus_clone = (TF1*) gausFit->Clone(Form("gausFit_%d", k));
      
      signal_clone->SetDirectory(0); // wichtig!

      
      //signal_clone->GetListOfFunctions()->Add(gaus_clone);
      proj->GetListOfFunctions()->Add(bg_clone);
      signal_clone->Write();
      proj->Write();
    }

    int pTsize = npT - 1 - n_pTstart;
    //plot Chi2/NPf RAW DATA FIT
    TH1D* h_chi2ndf = new TH1D("h_chi2ndf", "fit raw data", pTsize, &pT_bins[n_pTstart]);
    for (size_t i = 0; i < pTsize; ++i){
      h_chi2ndf->SetBinContent(i+1, chi2ndf_list[i]);
    }

    //plot Chi2/NPf SIGNAL FIT
    TH1D* h_chi2ndf_signalfit = new TH1D("h_chi2ndf_signalfit", "fit signal", pTsize, &pT_bins[n_pTstart]);
    for (size_t i = 0; i < pTsize; ++i){
      h_chi2ndf_signalfit->SetBinContent(i+1, chi2ndf_list_signalfit[i]);
    }

    //plot FIT Params SIGNAL
    TH1D* h_mean_signalfit = new TH1D("h_mean_signalfit", "fit raw data", pTsize, &pT_bins[n_pTstart]);
    for (size_t i = 0; i < pTsize; ++i){
      h_mean_signalfit->SetBinContent(i+1, mean_list_signalfit[i]);
      h_mean_signalfit->SetBinError(i, mean_error_list_signalfit[i]);
    }

    TH1D* h_sigma_signalfit = new TH1D("h_sigma_signalfit", "fit raw data", pTsize, &pT_bins[n_pTstart]);
    for (size_t i = 0; i < pTsize; ++i) {
      h_sigma_signalfit->SetBinContent(i+1, width_list_signalfit[i]);
      h_sigma_signalfit->SetBinError(i, width_error_list_signalfit[i]);
    }

    TH1D* h_FWHM_signalfit = new TH1D("h_FWHM_signalfit", "fit raw data", pTsize, &pT_bins[n_pTstart]);
    for (size_t i = 0; i < pTsize; ++i) {
      double y_FWHM = FWHM_Func(width_list_signalfit[i]);
      double y_FWHM_error = FWHM_Func(width_error_list_signalfit[i]);
      h_FWHM_signalfit->SetBinContent(i+1, y_FWHM);
      h_FWHM_signalfit->SetBinError(i+1, y_FWHM_error);
    }

    //Plot Counts
    TH1D* h_omegacounts = new TH1D("h_omegacounts", "Counts Omega Mesons", pTsize, &pT_bins[n_pTstart]);
    for(int i=0;i<pTsize;++i){
      std::cout<<"signal_integral_vec[i]: "<<signal_integral_vec[i]<<std::endl;
      double pT_bin_width = pT_bins[n_pTstart + i + 1] - pT_bins[n_pTstart + i];
      h_omegacounts->SetBinContent(i+1, signal_integral_vec[i] / pT_bin_width);
      h_omegacounts->SetBinError(i+1, signal_error_vec[i] / pT_bin_width);
    }

    std::cout<<"Writing..."<<std::endl;
    h_omegacounts->Write();
    h_FWHM_signalfit->Write();
    h_sigma_signalfit->Write();
    h_mean_signalfit->Write();
    h_chi2ndf_signalfit->Write();
    h_chi2ndf->Write();
    std::cout<<"Done"<<std::endl;


    input->Close();
    output->Close();
  }
  std::cout.rdbuf(originalCoutBuffer);
} 