// g++ -Wall -o analyzer `root-config --cflags --glibs` analyzer.cpp
//   ./analyzer pType Si120_HV800_Ele150   (as example)
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
#include "TProfile.h"
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

  //read input options
  std::string inputType = argv[1];
  std::string typeRun = argv[2];

  // std::string typeRun = "Si120_HV800_Ele150";
  // std::string typeRun = "Si200_HV800_Ele150";
  // std::string typeRun = "Si300_HV800_Ele150";

  // std::string typeRun = "Si120_HV800_Ele100";
  // std::string typeRun = "Si200_HV800_Ele100";
  // std::string typeRun = "Si300_HV800_Ele100";
  
  // std::string typeRun = "Si120_HV800_Pi150"; 
  // std::string typeRun = "Si300_HV800_Pi150"; 
  // std::string typeRun = "Si200_HV800_Pi150"; 


  //  std::cout << "----START ANALYZER: " 
  //	    << " inputFolder = " << inputFolder <<  " nRUN = " << nRUN << "-------" << std::endl;
 
  
  std::string inputFolder = "/store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTrees_v1/";
  if(inputType == "pType") inputFolder = "/store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/June2016/recoTrees_v1/";
  
  int ADCcutValue1 = 700;
  int ADCcutValueN = 700;


  std::vector<std::string> runN;
  if(typeRun == "Si120_HV800_Ele150" && inputType == "nType"){
    runN.push_back("3713");
    runN.push_back("3714");
    runN.push_back("3717");
    ADCcutValue1 = 340;
    ADCcutValueN = 255;
  }
  if(typeRun == "Si120_HV800_Ele100" && inputType == "pType"){
    runN.push_back("4174");
    runN.push_back("4175");
    runN.push_back("4176");
    runN.push_back("4177");
    runN.push_back("4178");
    runN.push_back("4179");
    runN.push_back("4180");
    runN.push_back("4181");
    runN.push_back("4182");
    runN.push_back("4183");
    runN.push_back("4184");
    ADCcutValue1 = 340;
    ADCcutValueN = 255;
  }

  if(typeRun == "Si200_HV800_Ele150" && inputType == "nType"){
    runN.push_back("3772");
    runN.push_back("3773");
    ADCcutValue1 = 460;
    ADCcutValueN = 345;
  }
  if(typeRun == "Si200_HV800_Ele100" && inputType == "pType"){
    runN.push_back("4187");
    runN.push_back("4188");
    runN.push_back("4189");
    runN.push_back("4190");
    runN.push_back("4191");
    runN.push_back("4195");
    runN.push_back("4196");
    runN.push_back("4197");
    runN.push_back("4198_1");
    runN.push_back("4199");
    runN.push_back("4200");
    runN.push_back("4202");
    ADCcutValue1 = 460;
    ADCcutValueN = 345;
  }

  if(typeRun == "Si300_HV800_Ele150" && inputType == "nType"){
    runN.push_back("3774");
    runN.push_back("3775");
    ADCcutValue1 = 700;
    ADCcutValueN = 525;
  }
  if(typeRun == "Si300_HV800_Ele100" && inputType == "pType"){
    runN.push_back("4203");
    runN.push_back("4204");
    runN.push_back("4205");
    runN.push_back("4206");
    runN.push_back("4207");
    runN.push_back("4208");
    runN.push_back("4209");
    runN.push_back("4210");
    runN.push_back("4211");
    ADCcutValue1 = 700;
    ADCcutValueN = 525;
  }


  // noise
  if(typeRun == "Si120_HV800_Pi150"){
    runN.push_back("3697");
  }
  if(typeRun == "Si200_HV800_Pi150"){
    runN.push_back("3663");
    runN.push_back("3670");
  }
  if(typeRun == "Si300_HV800_Pi150"){
    runN.push_back("3679");
  }


 
  TChain* tree = new TChain("H4treeReco");
  for(unsigned int iR=0; iR<runN.size(); ++iR){
    tree->Add(("root://eoscms/"+inputFolder+"RECO_"+runN.at(iR)+".root").c_str());
    std::cout << " Add: " << ("root://eoscms/"+inputFolder+"RECO_"+runN.at(iR)+".root").c_str() << std::endl;
  }

  Long64_t nentries = tree->GetEntries();
  std::cout << " Tree loaded events = " << nentries << std::endl;

  //Tree variables
  UInt_t nwc;
  UInt_t evtNumber;
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
  tree->SetBranchAddress("evtNumber", &evtNumber);
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


  //define reference diode and diodes order
  int ChRef = 4;
  int ChPadMapping[6] = {6, 5, 4, 1, 2, 3};
  if(inputType == "pType"){
    for(int iP=0; iP<6; ++iP) ChPadMapping[iP] = iP+1;
    ChRef = 1;
  }


  //Histos => needs to expand
  TH2F* WC_occupancy[2];
  for(int iw=0; iw<2; ++iw){
    WC_occupancy[iw] = new TH2F(Form("WC%d_occupancy",iw), "", 200, -20., 20., 200, -20., 20.);
  }

  TH2F* WC_diff = new TH2F("WC_diff", "", 20, -5.5, 1.5, 20, -3., 1.);
  TH2F* WC_diff_Pad1 = new TH2F("WC_diff_Pad1", "", 50, -20., 20., 50, -20., 20.);



  TH1F* h_TimeDiff[6];
  TProfile* tp_TimeDiff_12 = new TProfile("tp_TimeDiff_12", "", 10, 0., 30);
  TH2F* h2_TimeDiff_12 = new TH2F("h2_TimeDiff_12", "", 15, 0., 30., 1500, -2., 2.);
  TH2F* h2_TimeDiff_vsSoN[6];
  for(int iCh=0; iCh<6; ++iCh){
    h_TimeDiff[iCh] = new TH1F(Form("h_TimeDiff_Pad%d",ChPadMapping[iCh]), "", 1500, -2., 2.);
    h2_TimeDiff_vsSoN[iCh] = new TH2F(Form("h2_TimeDiff_vsSoN_Pad%d",ChPadMapping[iCh]),  "", 10, 0., 100., 1500, -2., 2.);
  }


  //loop over entries
  for (Long64_t jentry=0; jentry<nentries; ++jentry){
    
    //readout the event                                                               
    tree->GetEntry(jentry);
    

    //look at the position on the wire chambers
    WC_diff->Fill(wc_recox[0]-wc_recox[1], wc_recoy[0]-wc_recoy[1]);
    //FIXME    cut 
    for(unsigned int iw=0; iw<nwc; ++iw){  
      WC_occupancy[iw]->Fill(wc_recox[iw], wc_recoy[iw]);
    }

    
    for(int iPad=1; iPad<7; ++iPad){
      
      //wire chamber      
      if(abs(wc_recox[0]-wc_recox[1]+3.4) < 0.8 && abs(wc_recoy[0]-wc_recoy[1]+0.8) < 0.8 &&
	 wave_max[iPad] > ADCcutValueN && wave_max[ChRef] > ADCcutValue1){
	WC_diff_Pad1->Fill(wc_recox[0], wc_recoy[0]);
      }

      
      //timing
      if( wave_max[iPad] > ADCcutValueN && wave_max[ChRef] > ADCcutValue1){
	h_TimeDiff[iPad-1]->Fill(t_max_frac50[iPad] - t_max_frac50[ChRef]);
      }
      
      if(wave_max[ChRef] > ADCcutValue1)
	h2_TimeDiff_vsSoN[iPad-1]->Fill(wave_max[iPad]/pedestalRMS[iPad], t_max_frac50[iPad] - t_max_frac50[ChRef]);
    }//loop over pads
    
    
    //timing between un-irradiated pads
    // valid only for the 300mum thick
    if(wave_max[ChRef] > ADCcutValue1/20.*3. && wave_max[ChRef+1] > ADCcutValue1/20.*3.){
      if(wave_max[4] > wave_max[ChRef+1]) {
	
	tp_TimeDiff_12->Fill(wave_max[ChRef+1]/35., t_max_frac50[ChRef+1] - t_max_frac50[ChRef]);
	h2_TimeDiff_12->Fill(wave_max[ChRef+1]/35., t_max_frac50[ChRef+1] - t_max_frac50[ChRef]);
      }
      else{ 
	tp_TimeDiff_12->Fill(wave_max[ChRef]/35., t_max_frac50[ChRef+1] - t_max_frac50[ChRef]);
	h2_TimeDiff_12->Fill(wave_max[ChRef]/35., t_max_frac50[ChRef+1] - t_max_frac50[ChRef]);	
      }
    }
    

  }// loop over entries
  

  std::string outFolder = "plotFolder_"+typeRun;



  //Draw results occupancy plots
  TCanvas* cWC[2];
  for(int iw=0; iw<2; ++iw){
    cWC[iw] = new TCanvas();
    cWC[iw]->cd();
    WC_occupancy[iw]->GetXaxis()->SetTitle(Form("wc %d x (mm)", iw));
    WC_occupancy[iw]->GetYaxis()->SetTitle(Form("wc %d y (mm)", iw));
    WC_occupancy[iw]->Draw("colz");

    cWC[iw]->Print((outFolder+"/"+std::string(WC_occupancy[iw]->GetName())+".png").c_str(), ".png");
  }

  TCanvas* cWCdiff = new TCanvas();
  cWCdiff->cd();
  WC_diff->GetXaxis()->SetTitle("wc diff x (mm)");
  WC_diff->GetYaxis()->SetTitle("wc diff y (mm)");
  WC_diff->Draw("colz");
  cWCdiff->Print((outFolder+"/"+std::string(WC_diff->GetName())+".png").c_str(), ".png");

  TCanvas* cWCdiff_Pad1 = new TCanvas();
  cWCdiff_Pad1->cd();
  WC_diff_Pad1->GetXaxis()->SetTitle("wc x (mm)");
  WC_diff_Pad1->GetXaxis()->SetRangeUser(0., 15.);
  WC_diff_Pad1->GetYaxis()->SetTitle("wc y (mm)");
  WC_diff_Pad1->GetYaxis()->SetRangeUser(-20., -5.);
  //  WC_diff_Pad1->GetZaxis()->SetRangeUser(0., 200.);
  WC_diff_Pad1->Draw("colz");
  cWCdiff_Pad1->Print((outFolder+"/WC0_acceptance_"+typeRun+".png").c_str(), ".png");

  

  // build TGraph with results for timing
  //////////////////////////////////////////////////////////////////////////////
  double fluenceSet3S[6] = {0.00, 0.00, 4.00E+14, 6.00E+14, 6.00E+14, 9.00E+14};
  double fluenceSet2S[6] = {0.00, 0.00, 1.50E+15, 2.50E+15, 2.50E+15, 4.00E+15};
  double fluenceSet1S[6] = {0.00, 0.00, 6.25E+15, 6.25E+15, 1.00E+16, 1.60E+16};
  double updatedFluence[6];

  //graph for timing
  TGraphErrors* hTime = new TGraphErrors();
  hTime->SetName("hTime");
  for(int iC=0; iC<6; ++iC){
    if(typeRun == "Si120_HV800_Ele150" || typeRun == "Si120_HV800_Pi150") updatedFluence[iC] = fluenceSet1S[ChPadMapping[iC]-1];
    if(typeRun == "Si200_HV800_Ele150" || typeRun == "Si200_HV800_Pi150") updatedFluence[iC] = fluenceSet2S[ChPadMapping[iC]-1];
    if(typeRun == "Si300_HV800_Ele150" || typeRun == "Si300_HV800_Pi150") updatedFluence[iC] = fluenceSet3S[ChPadMapping[iC]-1];
    }

  hTime->SetPoint(0, 1.e+14, 0.001);
  TF1* unIr1_time = new TF1("unIr1_time", "pol0", 1., 1.e+17);
  TF1* tFit;
  TCanvas* cx;

  for(int iCh=0; iCh<6; ++iCh){
    if(iCh == ChRef-1 || h_TimeDiff[iCh]->GetEntries() < 30) continue;
    tFit = new TF1("tFit", "gaus", h_TimeDiff[iCh]->GetMean()-0.2, h_TimeDiff[iCh]->GetMean()+0.2);
    tFit->SetLineColor(kRed);
    tFit->SetLineWidth(2);
    tFit->SetRange(h_TimeDiff[iCh]->GetMean()-0.2, h_TimeDiff[iCh]->GetMean()+0.2);
    tFit->SetParameter(1, h_TimeDiff[iCh]->GetMean());
    tFit->SetParameter(2, 0.04);
    //    h_TimeDiff[iCh]->Rebin(2);
    
    cx = new TCanvas();
    cx->cd();
    //    h_TimeDiff[iCh]->GetXaxis()->SetRangeUser(-0.5, 0.5);
    h_TimeDiff[iCh]->GetXaxis()->SetRangeUser(h_TimeDiff[iCh]->GetMean()-0.15, h_TimeDiff[iCh]->GetMean()+0.15);
    h_TimeDiff[iCh]->Draw();

    h_TimeDiff[iCh]->Fit("tFit");
    if(iCh == ChRef) unIr1_time->SetParameter(0, tFit->GetParameter(2));
    else{
      hTime->SetPoint(ChPadMapping[iCh]-2, updatedFluence[iCh], tFit->GetParameter(2));
      hTime->SetPointError(ChPadMapping[iCh]-2, 0., tFit->GetParError(2));
    }
    // std::cout << " name = " << ChPadMapping[iCh] << std::endl;
    // std::cout << "iCh = " << iCh << " x = " << updatedFluence[iCh] << " >>> val = " 
    // 	      << tFit->GetParameter(2) << " err = " << tFit->GetParError(2) 
    // 	      << " 1/sqrt = " << 1./sqrt(h_TimeDiff[iCh]->GetEntries()) << std::endl;

    tFit->Draw("same");
    cx->Print(Form(("plotsTiming/"+inputType+"/fitTime_Pad%d_type_"+typeRun+".png").c_str(), ChPadMapping[iCh]), "png");
    cx->Delete();
    tFit->Delete();
  } 
  hTime->SetPoint(hTime->GetN()+1, 2.e+16, 2.);


  ///// time vs MIP  fro unirradiated diodes
  TGraphErrors* hTimevsMIP = new TGraphErrors();
  hTimevsMIP->SetPoint(0, 0., -1.);

  for(int iBin=1; iBin<h2_TimeDiff_12->GetNbinsX(); ++iBin){
    TH1D* hist12 = (TH1D*)h2_TimeDiff_12->ProjectionY(Form("hist12_bin%d",iBin), iBin, iBin+1);
    TF1* tFit12 = new TF1("tFit12", "gaus", hist12->GetMean()-0.2, hist12->GetMean()+0.2);
    TCanvas* cx12 = new TCanvas();
    cx12->cd();
    hist12->Fit("tFit12");
    hist12->Draw();
    tFit12->Draw("same");
    cx12->Print(Form(("plotsTiming12/"+inputType+"/fitTime_Pad12_bin%d_type_"+typeRun+".png").c_str(), iBin), "png");
    //    std::cout << " >>>  hTimevsMIP->GetN() = " << hTimevsMIP->GetN() << std::endl;
    hTimevsMIP->SetPoint(iBin, 2*iBin-1, tFit12->GetParameter(2));
    hTimevsMIP->SetPointError(iBin, 1., tFit12->GetParError(2));
    // std::cout << " >> x = " << 2*iBin-1 
    // 	      << " >> y = " << tFit12->GetParameter(2) << std::endl;
  }
  hTimevsMIP->SetPoint(hTimevsMIP->GetN()+1, 100., 2.);


  //time versus signal/Noise
  TGraphErrors* hTimevsSoN[6];
  for(int iCh=0; iCh<6; ++iCh){
    //    if(iCh == 3 ) continue;

    hTimevsSoN[iCh] = new TGraphErrors();
    hTimevsSoN[iCh]->SetName(Form("hTimevsSoN_pad%d", ChPadMapping[iCh]));
    hTimevsSoN[iCh]->SetPoint(0, 0., -1.);    

    for(int iBin=1; iBin<h2_TimeDiff_vsSoN[iCh]->GetNbinsX(); ++iBin){
      TH1D* hist12 = (TH1D*)h2_TimeDiff_vsSoN[iCh]->ProjectionY(Form(("histSoN_bin%d_Pad%d_type_"+typeRun).c_str(),iBin, ChPadMapping[iCh]), iBin, iBin+1);
      if(hist12->GetEntries() < 200)  {hTimevsSoN[iCh]->SetPoint(iBin, 5*(iBin+iBin-1), 2.); continue;}
      TF1* tFit12 = new TF1("tFit12", "gaus", hist12->GetMean()-0.2, hist12->GetMean()+0.2);
      TCanvas* cx12 = new TCanvas();
      cx12->cd();
      hist12->Fit("tFit12");
      hist12->Draw();
      tFit12->Draw("same");
      cx12->Print(Form(("plotsTimingSoN/"+inputType+"/fitTime_Pad%d_bin%d_type_"+typeRun+".png").c_str(), ChPadMapping[iCh], iBin), "png");
      // if(iBin == 9 && (ChPadMapping[iCh] == 2 || ChPadMapping[iCh] == 3 || ChPadMapping[iCh] == 4))
      // 	cx12->Print(Form(("plotsTimingSoN/fitTime_Pad%d_bin%d_type_"+typeRun+".root").c_str(), ChPadMapping[iCh], iBin), "root");
      hTimevsSoN[iCh]->SetPoint(iBin, 5*(iBin+iBin-1), tFit12->GetParameter(2));
      hTimevsSoN[iCh]->SetPointError(iBin, 1., tFit12->GetParError(2));      
    }
    hTimevsSoN[iCh]->SetPoint(hTimevsSoN[iCh]->GetN()+1, 100., 2.);
  }




  //save histos in out file (root format)
  TFile outF(("outFile_Timing_"+inputType+"_"+typeRun+".root").c_str(), "recreate");
  outF.cd();
  for(int iw=0; iw<2; ++iw){
    WC_occupancy[iw]->Write();
  }
  WC_diff->Write();
  WC_diff_Pad1->Write();

  for(int iCh=0; iCh<6; ++iCh){
    h_TimeDiff[iCh]->Write(h_TimeDiff[iCh]->GetName());

    hTimevsSoN[iCh]->Write(hTimevsSoN[iCh]->GetName());

    h2_TimeDiff_vsSoN[iCh]->Write();
  }

  hTime->Write();

  unIr1_time->Write();

  tp_TimeDiff_12->Write();
  hTimevsMIP->Write("hTimevsMIP");
  h2_TimeDiff_12->Write();

  outF.Close();

  return 0;     

}

