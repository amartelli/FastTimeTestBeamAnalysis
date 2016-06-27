#include "TLegend.h"
#include "TFile.h"
#include "TProfile.h"
#include "TLatex.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include "TTree.h"
#include "TChain.h"

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

#include <vector>
#include <fstream>
#include <string>

void doPlot_vsSoN()
{
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  std::string type = "nType";

  TFile* inF[3];
  inF[0] = TFile::Open(("outFile_Timing_"+type+"_Si120_HV800_Ele150.root").c_str());
  inF[1] = TFile::Open(("outFile_Timing_"+type+"_Si200_HV800_Ele150.root").c_str());
  inF[2] = TFile::Open(("outFile_Timing_"+type+"_Si300_HV800_Ele150.root").c_str());


  TGraphErrors* gE100[6];
  TGraphErrors* gE200[6];
  TGraphErrors* gE300[6];


  for(iP=0; iP<6; ++iP){
    gE100[iP] = (TGraphErrors*)inF[0]->Get(Form("hTimevsSoN_pad%d",iP+1));
    gE200[iP] = (TGraphErrors*)inF[1]->Get(Form("hTimevsSoN_pad%d",iP+1));
    gE300[iP] = (TGraphErrors*)inF[2]->Get(Form("hTimevsSoN_pad%d",iP+1));
  }


  std::string fluenceSet3[6] = {"0.00E+00",  "0.00E+00", "4.00E+14", "6.00E+14", "6.00E+14", "9.00E+14"};
  std::string fluenceSet2[6] = {"0.00E+00",  "0.00E+00", "1.50E+15", "2.50E+15", "2.50E+15", "4.00E+15"};
  std::string fluenceSet1[6] = {"0.00E+00",  "0.00E+00", "6.25E+15", "6.25E+15", "1.00E+16", "1.60E+16"};

  int colorsSet1[6] = {kBlack, kBlack, kRed, kRed, kBlue, kCyan};
  int colorsSet2[6] = {kBlack, kBlack, kRed, kBlue, kBlue, kCyan};
  int colorsSet3[6] = {kBlack, kBlack, kRed, kBlue, kBlue, kCyan};
  int colorsExt[6];
  int styleG[6] = {20, 23, 20, 23, 20, 23};
  float sizeG[6] = {1.1, 1.4, 1.1, 1.4, 1.1, 1.4};


  for(int iC=0; iC<3; ++iC){
    for(int cc=0; cc<6; ++cc){
      if(iC == 0) colorsExt[cc] = colorsSet1[cc];
      if(iC == 1) colorsExt[cc] = colorsSet2[cc];
      if(iC == 2) colorsExt[cc] = colorsSet3[cc];
      //    std::cout << " >>> iC = " << iC << " colorsExt[cc] = " << colorsExt[cc] << std::endl;
      
    if(iC == 0){
      gE100[cc]->SetMarkerColor(colorsExt[cc]);
      gE100[cc]->SetMarkerStyle(styleG[cc]);
      gE100[cc]->SetMarkerSize(sizeG[cc]);
    }
    if(iC == 1){
      gE200[cc]->SetMarkerColor(colorsExt[cc]);
      gE200[cc]->SetMarkerStyle(styleG[cc]);
      gE200[cc]->SetMarkerSize(sizeG[cc]);
    }
    if(iC == 2){
      gE300[cc]->SetMarkerColor(colorsExt[cc]);
      gE300[cc]->SetMarkerStyle(styleG[cc]);
      gE300[cc]->SetMarkerSize(sizeG[cc]);
    }
    }
  }



  std::string folderPlot = "plotsTimingSoN/"+type;

  // plots
  TLegend *legG1 = new TLegend(0.75,0.70,0.95,0.95,NULL,"brNDC");
  legG1->SetTextFont(42);
  legG1->SetTextSize(0.04);
  legG1->SetFillColor(kWhite);
  legG1->SetLineColor(kWhite);
  legG1->SetShadowColor(kWhite);
  legG1->SetHeader("Si P-on-N 120#mum");
  legG1->AddEntry(gE100[1], fluenceSet1[0].c_str(), "pl");
  legG1->AddEntry(gE100[2], fluenceSet1[2].c_str(), "pl");
  legG1->AddEntry(gE100[3], fluenceSet1[3].c_str(), "pl");
  legG1->AddEntry(gE100[4], fluenceSet1[4].c_str(), "pl");
  legG1->AddEntry(gE100[5], fluenceSet1[5].c_str(), "pl");


  TCanvas* graph1 = new TCanvas();
  graph1->cd();

  gE100[1]->GetXaxis()->SetTitle("S/N");
  gE100[1]->GetXaxis()->SetRangeUser(0, 100.);
  gE100[1]->GetYaxis()->SetTitle("#sigma(t_{PN} - t_{P1}) (ns)");
  gE100[1]->GetYaxis()->SetRangeUser(0.001, 0.2);
  gE100[1]->Draw("ap");
  for(int iC=2; iC<6; ++iC)
    gE100[iC]->Draw("p, same");
  legG1->Draw("same");
  graph1->Print((folderPlot+"/graph/timeResolutionSoN_N100.png").c_str(), "png");
  graph1->Print((folderPlot+"/graph/timeResolutionSoN_N100.root").c_str(), "root");


  //G2
  TLegend *legG2 = new TLegend(0.75,0.70,0.95,0.95,NULL,"brNDC");
  legG2->SetTextFont(42);
  legG2->SetTextSize(0.04);
  legG2->SetFillColor(kWhite);
  legG2->SetLineColor(kWhite);
  legG2->SetShadowColor(kWhite);
  legG2->SetHeader("Si P-on-N 200#mum");
  //  legG2->AddEntry(gE200[1], fluenceSet2[0].c_str(), "pl");
  legG2->AddEntry(gE200[2], fluenceSet2[2].c_str(), "pl");
  legG2->AddEntry(gE200[3], fluenceSet2[3].c_str(), "pl");
  legG2->AddEntry(gE200[4], fluenceSet2[4].c_str(), "pl");
  legG2->AddEntry(gE200[5], fluenceSet2[5].c_str(), "pl");


  TCanvas* graph2 = new TCanvas();
  gE200[2]->GetXaxis()->SetTitle("S/N");
  gE200[2]->GetXaxis()->SetRangeUser(0, 100.);
  gE200[2]->GetYaxis()->SetTitle("#sigma(t_{PN} - t_{P1}) (ns)");
  gE200[2]->GetYaxis()->SetRangeUser(0.001, 0.2);
  gE200[2]->Draw("ap");
  for(int iC=3; iC<6; ++iC)
    gE200[iC]->Draw("p, same");
  legG2->Draw("same");
  graph2->Print((folderPlot+"/graph/timeResolutionSoN_N200.png").c_str(), "png");
  graph2->Print((folderPlot+"/graph/timeResolutionSoN_N200.root").c_str(), "root");


  // G3
  TLegend *legG3 = new TLegend(0.75,0.70,0.95,0.95,NULL,"brNDC");
  legG3->SetTextFont(42);
  legG3->SetTextSize(0.04);
  legG3->SetFillColor(kWhite);
  legG3->SetLineColor(kWhite);
  legG3->SetShadowColor(kWhite);
  legG3->SetHeader("Si P-on-N 320#mum");
  legG3->AddEntry(gE300[1], fluenceSet3[0].c_str(), "pl");
  legG3->AddEntry(gE300[2], fluenceSet3[2].c_str(), "pl");
  legG3->AddEntry(gE300[3], fluenceSet3[3].c_str(), "pl");
  legG3->AddEntry(gE300[4], fluenceSet3[4].c_str(), "pl");
  legG3->AddEntry(gE300[5], fluenceSet3[5].c_str(), "pl");


  TCanvas* graph3 = new TCanvas();
  gE300[1]->GetXaxis()->SetTitle("S/N");
  gE300[1]->GetXaxis()->SetRangeUser(0, 100.);
  gE300[1]->GetYaxis()->SetTitle("#sigma(t_{PN} - t_{P1}) (ns)");
  gE300[1]->GetYaxis()->SetRangeUser(0.001, 0.2);
  gE300[1]->Draw("ap");
  for(int iC=2; iC<6; ++iC)
    gE300[iC]->Draw("p, same");
  legG3->Draw("same");
  graph3->Print((folderPlot+"/graph/timeResolutionSoN_N300.png").c_str(), "png");
  graph3->Print((folderPlot+"/graph/timeResolutionSoN_N300.root").c_str(), "root");


}
