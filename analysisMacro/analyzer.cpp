// g++ -Wall -o analyzer `root-config --cflags --glibs` analyzer.cpp
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TROOT.h"
#include "TSystem.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TString.h"
#include "TCut.h"
#include "TMath.h"
#include "TApplication.h"
#include "TError.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPad.h"


int main(int argc, char** argv)
{
  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");

  //read input options => need to expand
  //  std::string inputFolder = argv[1];
  //  std::string nRUN = argv[2];
  //  std::cout << "----START ANALYZER: " 
  //	    << " inputFolder = " << inputFolder <<  " nRUN = " << nRUN << "-------" << std::endl;
 
  
  // need to set these configurable from cfg and change into a vector
  std::string inputFolder = "/store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTrees/";
  std::string runN = "3772";

  TChain* tree = new TChain("H4treeReco");
  tree->Add(("root://eoscms/"+inputFolder+"RECO_"+runN+".root").c_str());

  Long64_t nentries = tree->GetEntries();
  std::cout << " Tree loaded events = " << nentries << std::endl;

  //Tree variables
  UInt_t nwc;
  Float_t wc_recox[16], wc_recoy[16];
  UInt_t maxch;
  Float_t group[100],ch[100];
  Float_t pedestal[100], pedestalRMS[100], pedestalSlope[100];
  Float_t wave_max[100], wave_max_aft[100], wave_aroundmax[100][50], time_aroundmax[100][50];
  Float_t charge_integ[100], charge_integ_max[100], charge_integ_fix[100];
  Float_t charge_integ_smallw[100], charge_integ_smallw_noise[100], charge_integ_largew[100], charge_integ_largew_noise[100];
  Float_t t_max[100], t_max_frac30[100], t_max_frac50[100], t_at_threshold[100], t_over_threshold[100];

  //Read tree
  tree->SetBranchAddress("nwc",       &nwc);
  tree->SetBranchAddress("wc_recox",   wc_recox);
  tree->SetBranchAddress("wc_recoy",   wc_recoy);

  tree->SetBranchAddress("maxch",               &maxch);
  tree->SetBranchAddress("group",                group);
  tree->SetBranchAddress("ch",                   ch);
  tree->SetBranchAddress("pedestal",             pedestal);
  tree->SetBranchAddress("pedestalRMS",          pedestalRMS);
  tree->SetBranchAddress("pedestalSlope",        pedestalSlope);
  tree->SetBranchAddress("wave_max",             wave_max);
  tree->SetBranchAddress("wave_max_aft",         wave_max_aft);
  tree->SetBranchAddress("wave_aroundmax",       wave_aroundmax);
  tree->SetBranchAddress("time_aroundmax",       time_aroundmax);

  tree->SetBranchAddress("charge_integ",         charge_integ);
  tree->SetBranchAddress("charge_integ_max",     charge_integ_max);
  tree->SetBranchAddress("charge_integ_fix",     charge_integ_fix);
  tree->SetBranchAddress("charge_integ_smallw",  charge_integ_smallw);
  tree->SetBranchAddress("charge_integ_largew",  charge_integ_largew);
  tree->SetBranchAddress("charge_integ_smallw_noise",  charge_integ_smallw_noise);
  tree->SetBranchAddress("charge_integ_largew_noise",  charge_integ_largew_noise);
  tree->SetBranchAddress("t_max",                t_max);
  tree->SetBranchAddress("t_max_frac30",         t_max_frac30);
  tree->SetBranchAddress("t_max_frac50",         t_max_frac50);
  tree->SetBranchAddress("t_at_threshold",       t_at_threshold);
  tree->SetBranchAddress("t_over_threshold",     t_over_threshold);

  //Histos => needs to expand
  TH2F* WC_occupancy[2];
  for(int iw=0; iw<2; ++iw){
    WC_occupancy[iw] = new TH2F(Form("WC%d_occupancy",iw), "", 400, -20., 20., 400, -20., 20.);
  }


  //loop over entries
  for (Long64_t jentry=0; jentry<nentries; ++jentry){

    //readout the event                                                                                                                                                        
    tree->GetEntry(jentry);


    for(unsigned int iw=0; iw<nwc; ++iw){
   
      WC_occupancy[iw]->Fill(wc_recox[iw], wc_recoy[iw]);
    }
  }


  std::string outFolder = "plotFolder_"+runN;

  //Draw results 
  TCanvas* cWC[2];
  for(int iw=0; iw<2; ++iw){
    cWC[iw] = new TCanvas();
    cWC[iw]->cd();
    WC_occupancy[iw]->GetXaxis()->SetTitle(Form("wc %d x (mm)", iw));
    WC_occupancy[iw]->GetYaxis()->SetTitle(Form("wc %d y (mm)", iw));
    WC_occupancy[iw]->Draw("colz");

    cWC[iw]->Print((outFolder+"/"+std::string(WC_occupancy[iw]->GetName())+".png").c_str(), ".png");
  }

  return 0;     

}

