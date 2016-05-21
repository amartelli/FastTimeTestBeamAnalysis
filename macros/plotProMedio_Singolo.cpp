// g++ -Wall -o plotProMedio_Singolo `root-config --cflags --glibs`  plotProMedio_Singolo.cpp

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
#include <vector>
#include <fstream>
#include <string>

//#include "../interface/Waveform.hpp"
#include "waveFormUtils.h"


int main(){
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* inF[3];
  inF[0] = TFile::Open("Above200Thr/waveTemplates2016_120um.root");
  inF[1] = TFile::Open("Above200Thr/waveTemplates2016_200um.root");
  inF[2] = TFile::Open("Above200Thr/waveTemplates2016_300um.root");

  //vs gen                                                                                                                                                                      

  TProfile* tp_Pad1[3];
  TProfile* tp_Pad2[3];
  TProfile* tp_Pad3[3];
  TProfile* tp_Pad4[3];
  TProfile* tp_Pad5[3];
  TProfile* tp_Pad6[3];

  //set3
  std::string fluenceSet3[6] = {"0.00E+00",  "0.00E+00", "4.00E+14", "6.00E+14", "6.00E+14", "9.00E+14"};
  int colorsSet3[6] = {kBlack, kBlack, kRed, kBlue, kBlue, kCyan};
  int colorsSetG3[6] = {kBlack, kBlack, kRed, kBlue, kBlue, kCyan};
  //set2
  std::string fluenceSet2[6] = {"0.00E+00",  "0.00E+00", "1.50E+15", "2.50E+15", "2.50E+15", "4.00E+15"};
  //0.00E+00   0.00E+00   1.50E+15   2.50E+15   2.50E+15   4.00E+15
  int colorsSet2[6] = {kBlack, kBlack, kRed, kBlue, kBlue, kCyan};
  int colorsSetG2[6] = {kBlack, kBlack, kGreen, kViolet, kViolet, kMagenta};
  //set1
  std::string fluenceSet1[6] = {"0.00E+00",  "0.00E+00", "6.25E+15", "6.25E+15", "1.00E+16", "1.60E+16"};
  //0.00E+00   0.00E+00   6.25E+15   6.25E+15   1.00E+16   1.60E+16   
  int colorsSet1[6] = {kBlack, kBlack, kRed, kRed, kBlue, kCyan};
  int colorsSetG1[6] = {kBlack, kBlack, kOrange+1, kOrange+1, kGreen+2, kRed+2};
  int colorsExt[6];

  int styleM[6] = {7, 7, 7, 7, 7, 7};
  int styleG[6] = {20, 23, 20, 23, 20, 23};
  float sizeG[6] = {1.1, 1.4, 1.1, 1.4, 1.1, 1.4};

  for(int iC=0; iC<3; ++iC){
    std::cout << "iC = " << iC << std::endl;
    tp_Pad1[iC] = (TProfile*)(inF[iC]->Get("SiPad1_waveProfile"));
    tp_Pad2[iC] = (TProfile*)(inF[iC]->Get("SiPad2_waveProfile"));
    tp_Pad3[iC] = (TProfile*)(inF[iC]->Get("SiPad3_waveProfile"));
    tp_Pad4[iC] = (TProfile*)(inF[iC]->Get("SiPad4_waveProfile"));
    tp_Pad5[iC] = (TProfile*)(inF[iC]->Get("SiPad5_waveProfile"));
    tp_Pad6[iC] = (TProfile*)(inF[iC]->Get("SiPad6_waveProfile"));

    for(int cc=0; cc<6; ++cc){
    if(iC == 0) colorsExt[cc] = colorsSet1[cc];
    if(iC == 1) colorsExt[cc] = colorsSet2[cc];
    if(iC == 2) colorsExt[cc] = colorsSet3[cc];
    std::cout << " >>> iC = " << iC << " colorsExt[cc] = " << colorsExt[cc] << std::endl;
    }
    tp_Pad1[iC]->SetLineColor(colorsExt[0]);
    tp_Pad2[iC]->SetLineColor(colorsExt[1]);
    tp_Pad3[iC]->SetLineColor(colorsExt[2]);
    tp_Pad4[iC]->SetLineColor(colorsExt[3]);
    tp_Pad5[iC]->SetLineColor(colorsExt[4]);
    tp_Pad6[iC]->SetLineColor(colorsExt[5]);

    tp_Pad1[iC]->SetMarkerColor(colorsExt[0]);
    tp_Pad2[iC]->SetMarkerColor(colorsExt[1]);
    tp_Pad3[iC]->SetMarkerColor(colorsExt[2]);
    tp_Pad4[iC]->SetMarkerColor(colorsExt[3]);
    tp_Pad5[iC]->SetMarkerColor(colorsExt[4]);
    tp_Pad6[iC]->SetMarkerColor(colorsExt[5]);

    tp_Pad1[iC]->SetMarkerStyle(styleM[0]);
    tp_Pad2[iC]->SetMarkerStyle(styleM[1]);
    tp_Pad3[iC]->SetMarkerStyle(styleM[2]);
    tp_Pad4[iC]->SetMarkerStyle(styleM[3]);
    tp_Pad5[iC]->SetMarkerStyle(styleM[4]);
    tp_Pad6[iC]->SetMarkerStyle(styleM[5]);

    tp_Pad1[iC]->SetLineWidth(2);
    tp_Pad2[iC]->SetLineWidth(2);
    tp_Pad3[iC]->SetLineWidth(2);
    tp_Pad4[iC]->SetLineWidth(2);
    tp_Pad5[iC]->SetLineWidth(2);
    tp_Pad6[iC]->SetLineWidth(2);
    
    /*
    tp_Pad1[iC]->Rebin(3);
    tp_Pad2[iC]->Rebin(3);
    tp_Pad3[iC]->Rebin(3);
    tp_Pad4[iC]->Rebin(3);
    tp_Pad5[iC]->Rebin(3);
    tp_Pad6[iC]->Rebin(3);
    */
  }


  std::vector<float> max_Pad1;
  std::vector<float> max_Pad2;
  std::vector<float> max_Pad3;
  std::vector<float> max_Pad4;
  std::vector<float> max_Pad5;
  std::vector<float> max_Pad6;

  std::vector<float> at30_Pad1;
  std::vector<float> at30_Pad2;
  std::vector<float> at30_Pad3;
  std::vector<float> at30_Pad4;
  std::vector<float> at30_Pad5;
  std::vector<float> at30_Pad6;
  std::vector<float> atMax_Pad1;
  std::vector<float> atMax_Pad2;
  std::vector<float> atMax_Pad3;
  std::vector<float> atMax_Pad4;
  std::vector<float> atMax_Pad5;
  std::vector<float> atMax_Pad6;


  std::vector<float> times_P1;
  std::vector<float> times_P2;
  std::vector<float> times_P3;
  std::vector<float> times_P4;
  std::vector<float> times_P5;
  std::vector<float> times_P6;
  std::vector<float> samples_P1;
  std::vector<float> samples_P2;
  std::vector<float> samples_P3;
  std::vector<float> samples_P4;
  std::vector<float> samples_P5;
  std::vector<float> samples_P6;

  TProfile* templateP1[3];
  TProfile* templateP2[3];
  TProfile* templateP3[3];
  TProfile* templateP4[3];
  TProfile* templateP5[3];
  TProfile* templateP6[3];
  for(int iSet=0; iSet<3; ++iSet){
    templateP1[iSet] = new TProfile(Form("SiPad1_waveProfile_Set%d", iSet+1), "", 800, -50., 50.);
    templateP2[iSet] = new TProfile(Form("SiPad2_waveProfile_Set%d", iSet+1), "", 800, -50., 50.);
    templateP3[iSet] = new TProfile(Form("SiPad3_waveProfile_Set%d", iSet+1), "", 800, -50., 50.);
    templateP4[iSet] = new TProfile(Form("SiPad4_waveProfile_Set%d", iSet+1), "", 800, -50., 50.);
    templateP5[iSet] = new TProfile(Form("SiPad5_waveProfile_Set%d", iSet+1), "", 800, -50., 50.);
    templateP6[iSet] = new TProfile(Form("SiPad6_waveProfile_Set%d", iSet+1), "", 800, -50., 50.);
  }

  for(int iC=0; iC<3; ++iC){
    times_P1.clear();
    times_P2.clear();
    times_P3.clear();
    times_P4.clear();
    times_P5.clear();
    times_P6.clear();
    samples_P1.clear();
    samples_P2.clear();
    samples_P3.clear();
    samples_P4.clear();
    samples_P5.clear();
    samples_P6.clear();
    for(int i=0; i<tp_Pad1[0]->GetNbinsX(); ++i){
      times_P1.push_back(tp_Pad1[iC]->GetBinCenter(i+1));
      samples_P1.push_back(tp_Pad1[iC]->GetBinContent(i+1));
      times_P2.push_back(tp_Pad2[iC]->GetBinCenter(i+1));
      samples_P2.push_back(tp_Pad2[iC]->GetBinContent(i+1));
      times_P3.push_back(tp_Pad3[iC]->GetBinCenter(i+1));
      samples_P3.push_back(tp_Pad3[iC]->GetBinContent(i+1));
      times_P4.push_back(tp_Pad3[iC]->GetBinCenter(i+1));
      samples_P4.push_back(tp_Pad4[iC]->GetBinContent(i+1));
      times_P5.push_back(tp_Pad5[iC]->GetBinCenter(i+1));
      samples_P5.push_back(tp_Pad5[iC]->GetBinContent(i+1));
      times_P6.push_back(tp_Pad6[iC]->GetBinCenter(i+1));
      samples_P6.push_back(tp_Pad6[iC]->GetBinContent(i+1));
    }
    Waveform templatePS_P1(times_P1, samples_P1);
    Waveform::max_amplitude_informations wave_max_P1 = templatePS_P1.max_amplitude(240, 560, 5);
    at30_Pad1.push_back(templatePS_P1.time_at_frac(240, 560, 0.1, wave_max_P1, 5));
    atMax_Pad1.push_back(templatePS_P1.time_at_frac(240, 560, 0.9, wave_max_P1, 5));
    // at30_Pad1.push_back(templatePS_P1.time_at_frac(240, 560, 0.3, wave_max_P1, 5));
    // atMax_Pad1.push_back(wave_max_P1.time_at_max);
    max_Pad1.push_back(wave_max_P1.max_amplitude);


    Waveform templatePS_P2(times_P2, samples_P2);
    Waveform::max_amplitude_informations wave_max_P2 = templatePS_P2.max_amplitude(240, 560, 5);
    at30_Pad2.push_back(templatePS_P2.time_at_frac(240, 560, 0.1, wave_max_P2, 5));
    atMax_Pad2.push_back(templatePS_P2.time_at_frac(240, 560, 0.9, wave_max_P2, 5));
    //at30_Pad2.push_back(templatePS_P2.time_at_frac(240, 560, 0.3, wave_max_P2, 5));
    //atMax_Pad2.push_back(wave_max_P2.time_at_max);
    max_Pad2.push_back(wave_max_P2.max_amplitude);

    Waveform templatePS_P3(times_P3, samples_P3);
    Waveform::max_amplitude_informations wave_max_P3 = templatePS_P3.max_amplitude(240, 560, 5);
    at30_Pad3.push_back(templatePS_P3.time_at_frac(240, 560, 0.1, wave_max_P3, 5));
    atMax_Pad3.push_back(templatePS_P3.time_at_frac(240, 560, 0.9, wave_max_P3, 5));
    //at30_Pad3.push_back(templatePS_P3.time_at_frac(240, 560, 0.3, wave_max_P3, 5));
    //atMax_Pad3.push_back(wave_max_P3.time_at_max);
    max_Pad3.push_back(wave_max_P3.max_amplitude);

    Waveform templatePS_P4(times_P4, samples_P4);
    Waveform::max_amplitude_informations wave_max_P4 = templatePS_P4.max_amplitude(240, 560, 5);
    at30_Pad4.push_back(templatePS_P4.time_at_frac(240, 560, 0.1, wave_max_P4, 5));
    atMax_Pad4.push_back(templatePS_P4.time_at_frac(240, 560, 0.9, wave_max_P4, 5));
    //at30_Pad4.push_back(templatePS_P4.time_at_frac(240, 560, 0.3, wave_max_P4, 5));
    //atMax_Pad4.push_back(wave_max_P4.time_at_max);
    max_Pad4.push_back(wave_max_P4.max_amplitude);

    Waveform templatePS_P5(times_P5, samples_P5);
    Waveform::max_amplitude_informations wave_max_P5 = templatePS_P5.max_amplitude(240, 560, 5);
    at30_Pad5.push_back(templatePS_P5.time_at_frac(240, 560, 0.1, wave_max_P5, 5));
    atMax_Pad5.push_back(templatePS_P5.time_at_frac(240, 560, 0.9, wave_max_P5, 5));
    //at30_Pad5.push_back(templatePS_P5.time_at_frac(240, 560, 0.3, wave_max_P5, 5));
    //atMax_Pad5.push_back(wave_max_P5.time_at_max);
    max_Pad5.push_back(wave_max_P5.max_amplitude);

    Waveform templatePS_P6(times_P6, samples_P6);
    Waveform::max_amplitude_informations wave_max_P6 = templatePS_P6.max_amplitude(240, 560, 5);
    at30_Pad6.push_back(templatePS_P6.time_at_frac(240, 560, 0.1, wave_max_P6, 5));
    atMax_Pad6.push_back(templatePS_P6.time_at_frac(240, 560, 0.9, wave_max_P6, 5));
    //at30_Pad6.push_back(templatePS_P6.time_at_frac(240, 560, 0.3, wave_max_P6, 5));
    //atMax_Pad6.push_back(wave_max_P6.time_at_max);
    max_Pad6.push_back(wave_max_P6.max_amplitude);

    for(int i=0; i<tp_Pad1[0]->GetNbinsX(); ++i){
      templateP1[iC]->Fill(tp_Pad1[iC]->GetBinCenter(i+1)/0.2, tp_Pad1[iC]->GetBinContent(i+1)/wave_max_P1.max_amplitude);
      templateP2[iC]->Fill(tp_Pad2[iC]->GetBinCenter(i+1)/0.2, tp_Pad2[iC]->GetBinContent(i+1)/wave_max_P2.max_amplitude);
      templateP3[iC]->Fill(tp_Pad3[iC]->GetBinCenter(i+1)/0.2, tp_Pad3[iC]->GetBinContent(i+1)/wave_max_P3.max_amplitude);
      templateP4[iC]->Fill(tp_Pad4[iC]->GetBinCenter(i+1)/0.2, tp_Pad4[iC]->GetBinContent(i+1)/wave_max_P4.max_amplitude);
      templateP5[iC]->Fill(tp_Pad5[iC]->GetBinCenter(i+1)/0.2, tp_Pad5[iC]->GetBinContent(i+1)/wave_max_P5.max_amplitude);
      templateP6[iC]->Fill(tp_Pad6[iC]->GetBinCenter(i+1)/0.2, tp_Pad6[iC]->GetBinContent(i+1)/wave_max_P6.max_amplitude);
    }

  }

  TFile outTemplate("AllwaveTemplate.root", "recreate");
  outTemplate.cd();
  for(int iSet=0; iSet<3; ++iSet){
    templateP1[iSet]->Write(templateP1[iSet]->GetName());
    templateP2[iSet]->Write(templateP2[iSet]->GetName());
    templateP3[iSet]->Write(templateP3[iSet]->GetName());
    templateP4[iSet]->Write(templateP4[iSet]->GetName());
    templateP5[iSet]->Write(templateP5[iSet]->GetName());
    templateP6[iSet]->Write(templateP6[iSet]->GetName());
  }
  outTemplate.Close();


  TGraph* PS = new TGraph();
  PS->SetPoint(0, -1, -1);
  PS->SetPoint(1, 0, atMax_Pad1.at(0) - at30_Pad1.at(0));
  PS->SetPoint(2, 6.25e+15, atMax_Pad3.at(0) - at30_Pad3.at(0));
  PS->SetPoint(3, 6.25e+15, atMax_Pad4.at(0) - at30_Pad4.at(0));
  PS->SetPoint(4, 10e+15, atMax_Pad5.at(0) - at30_Pad5.at(0));
  PS->SetPoint(5, 16e+15, atMax_Pad6.at(0) - at30_Pad6.at(0));
  PS->SetPoint(6, 30e+15, 2);
  PS->SetMarkerStyle(20);

  TGraph* P1[3];
  TGraph* P2[3];
  TGraph* P3[3];
  TGraph* P4[3];
  TGraph* P5[3];
  TGraph* P6[3];
  for(int iC=0; iC<3; ++iC){
    P1[iC] = new TGraph();
    P2[iC] = new TGraph();
    P3[iC] = new TGraph();
    P4[iC] = new TGraph();
    P5[iC] = new TGraph();
    P6[iC] = new TGraph();
  }
  for(int iC=0; iC<3; ++iC){
    P1[iC]->SetPoint(0, 0, 0);
    P2[iC]->SetPoint(0, 0, 0);
    P3[iC]->SetPoint(0, 0, 0);
    P4[iC]->SetPoint(0, 0, 0);
    P5[iC]->SetPoint(0, 0, 0);
    P6[iC]->SetPoint(0, 0, 0);
  }
  for(int iC=0; iC<3; ++iC){
    P1[iC]->SetPoint(iC+1, (iC+1)*100, atMax_Pad1.at(iC) - at30_Pad1.at(iC));
    P2[iC]->SetPoint(iC+1, (iC+1)*100, atMax_Pad2.at(iC) - at30_Pad2.at(iC));
    P3[iC]->SetPoint(iC+1, (iC+1)*100, atMax_Pad3.at(iC) - at30_Pad3.at(iC));
    P4[iC]->SetPoint(iC+1, (iC+1)*100, atMax_Pad4.at(iC) - at30_Pad4.at(iC));
    P5[iC]->SetPoint(iC+1, (iC+1)*100, atMax_Pad5.at(iC) - at30_Pad5.at(iC));
    P6[iC]->SetPoint(iC+1, (iC+1)*100, atMax_Pad6.at(iC) - at30_Pad6.at(iC));
  }
  for(int iC=0; iC<3; ++iC){
    P1[iC]->SetPoint(4, 500, 2);
    P2[iC]->SetPoint(4, 500, 2);
    P3[iC]->SetPoint(4, 500, 2);
    P4[iC]->SetPoint(4, 500, 2);
    P5[iC]->SetPoint(4, 500, 2);
    P6[iC]->SetPoint(4, 500, 2);
  }
  
  for(int iC=0; iC<3; ++iC){
    for(int cc=0; cc<6; ++cc){
      if(iC == 0) colorsExt[cc] = colorsSetG1[cc];
      if(iC == 1) colorsExt[cc] = colorsSetG2[cc];
      if(iC == 2) colorsExt[cc] = colorsSetG3[cc];
      std::cout << " >>> iC = " << iC << " colorsExt[cc] = " << colorsExt[cc] << std::endl;
    }
    P1[iC]->SetMarkerColor(colorsExt[0]);
    P2[iC]->SetMarkerColor(colorsExt[1]);
    P3[iC]->SetMarkerColor(colorsExt[2]);
    P4[iC]->SetMarkerColor(colorsExt[3]);
    P5[iC]->SetMarkerColor(colorsExt[4]);
    P6[iC]->SetMarkerColor(colorsExt[5]);
    P1[iC]->SetMarkerStyle(styleG[0]);
    P2[iC]->SetMarkerStyle(styleG[1]);
    P3[iC]->SetMarkerStyle(styleG[2]);
    P4[iC]->SetMarkerStyle(styleG[3]);
    P5[iC]->SetMarkerStyle(styleG[4]);
    P6[iC]->SetMarkerStyle(styleG[5]);
    P1[iC]->SetMarkerSize(sizeG[0]);
    P2[iC]->SetMarkerSize(sizeG[1]);
    P3[iC]->SetMarkerSize(sizeG[2]);
    P4[iC]->SetMarkerSize(sizeG[3]);
    P5[iC]->SetMarkerSize(sizeG[4]);
    P6[iC]->SetMarkerSize(sizeG[5]);
  }
  

  TLegend *legG1 = new TLegend(0.75,0.70,0.95,0.95,NULL,"brNDC");
  legG1->SetTextFont(42);
  legG1->SetTextSize(0.04);
  legG1->SetFillColor(kWhite);
  legG1->SetLineColor(kWhite);
  legG1->SetShadowColor(kWhite);
  legG1->SetHeader("Si N-on-P 120#mum");
  legG1->AddEntry(P1[0], fluenceSet1[0].c_str(), "pl");
  legG1->AddEntry(P3[0], fluenceSet1[2].c_str(), "pl");
  legG1->AddEntry(P4[0], fluenceSet1[3].c_str(), "pl");
  legG1->AddEntry(P5[0], fluenceSet1[4].c_str(), "pl");
  legG1->AddEntry(P6[0], fluenceSet1[5].c_str(), "pl");

  //  legG->AddEntry(tp_Pad2[0], fluenceSet1[1].c_str(), "pl");
  TLegend *legG2 = new TLegend(0.75,0.43,0.95,0.68,NULL,"brNDC");
  legG2->SetTextFont(42);
  legG2->SetTextSize(0.04);
  legG2->SetFillColor(kWhite);
  legG2->SetLineColor(kWhite);
  legG2->SetShadowColor(kWhite);
  legG2->SetHeader("Si N-on-P 200#mum");
  legG2->AddEntry(P1[1], fluenceSet2[0].c_str(), "pl");
  //  legG->AddEntry(P2[0], fluenceSet1[1].c_str(), "pl");
  legG2->AddEntry(P3[1], fluenceSet2[2].c_str(), "pl");
  legG2->AddEntry(P4[1], fluenceSet2[3].c_str(), "pl");
  legG2->AddEntry(P5[1], fluenceSet2[4].c_str(), "pl");
  legG2->AddEntry(P6[1], fluenceSet2[5].c_str(), "pl");

  //
  TLegend *legG3 = new TLegend(0.75,0.16,0.95,0.41,NULL,"brNDC");
  legG3->SetTextFont(42);
  legG3->SetTextSize(0.04);
  legG3->SetFillColor(kWhite);
  legG3->SetLineColor(kWhite);
  legG3->SetShadowColor(kWhite);
  legG3->SetHeader("Si N-on-P 320#mum");
  //
  legG3->AddEntry(P1[2], fluenceSet3[0].c_str(), "pl");
  //  legG->AddEntry(P2[0], fluenceSet1[1].c_str(), "pl");
  legG3->AddEntry(P3[2], fluenceSet3[2].c_str(), "pl");
  legG3->AddEntry(P4[2], fluenceSet3[3].c_str(), "pl");
  legG3->AddEntry(P5[2], fluenceSet3[4].c_str(), "pl");
  legG3->AddEntry(P6[2], fluenceSet3[5].c_str(), "pl");


  TCanvas* graph = new TCanvas();
  P1[0]->GetXaxis()->SetTitle("Si thickness (#mum)");
  P1[0]->GetXaxis()->SetRangeUser(80, 380.);
  P1[0]->GetYaxis()->SetTitle("t_{max} - t_{30%} (ns)");
  P1[0]->GetYaxis()->SetRangeUser(0.95, 1.15);
  P1[0]->Draw("ap");
  //  P2[0]->Draw("p, same");
  P3[0]->Draw("p, same");
  P4[0]->Draw("p, same");
  P5[0]->Draw("p, same");
  P6[0]->Draw("p, same");
  for(int iC=1; iC<3; ++iC){
    P1[iC]->Draw("p, same");
    //    P2[iC]->Draw("p, same");
    P3[iC]->Draw("p, same");
    P4[iC]->Draw("p, same");
    P5[iC]->Draw("p, same");
    P6[iC]->Draw("p, same");
  }
  legG1->Draw("same");
  legG2->Draw("same");
  legG3->Draw("same");
  // graph->Print("plots_ProMediAbove200_Single/timeRise_pulseShapes_10_90.png", "png");
  // graph->Print("plots_ProMediAbove200_Single/timeRise_pulseShapes_10_90.root", ".root");
  graph->Print("plots_ProMediAbove200_Single/timeRise_pulseShapes_30_100.png", "png");
  graph->Print("plots_ProMediAbove200_Single/timeRise_pulseShapes_30_100.root", ".root");

  ///singolo
  TLegend *legS = new TLegend(0.6,0.65,0.8,0.9,NULL,"brNDC");
  legS->SetTextFont(42);
  legS->SetTextSize(0.04);
  legS->SetFillColor(kWhite);
  legS->SetLineColor(kWhite);
  legS->SetShadowColor(kWhite);
  legS->SetHeader("Si N-on-P 120#mum");
  
  TCanvas* graphS = new TCanvas();
  //  PS->GetXaxis()->SetTitle("fluence (n/cm^{2} x E+15)");
  PS->GetXaxis()->SetTitle("fluence (n/cm^{2})");
  PS->GetXaxis()->SetRangeUser(-1, 20.);
  PS->GetYaxis()->SetTitle("t_{90%} - t_{10%} (ns)");
  PS->GetYaxis()->SetRangeUser(0.9, 1.1);
  PS->Draw("ap");
  legS->Draw("same");
  graphS->Print("plots_ProMediAbove200_Single/timeRise_pulseShapes_120_90_10.png", "png");
  graphS->Print("plots_ProMediAbove200_Single/timeRise_pulseShapes_120_90_10.root", ".root");





  TFile outRoot("outRoot.root", "recreate");
  outRoot.cd();
  for(int iC=0; iC<3; ++iC){
    P1[iC]->Write();
    P3[iC]->Write();
    P4[iC]->Write();
    P5[iC]->Write();
    P6[iC]->Write();
  }
  outRoot.Close();
  //  return 10;
  std::cout << " ci sono " << std::endl;

  //gStyle->SetLegendTextSize(20);
  TLegend *leg1 = new TLegend(0.6,0.65,0.8,0.90,NULL,"brNDC");
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.04);
  leg1->SetFillColor(kWhite);
  leg1->SetLineColor(kWhite);
  leg1->SetShadowColor(kWhite);
  leg1->SetHeader("Si N-on-P 120#mum");
  leg1->AddEntry(tp_Pad1[0], fluenceSet1[0].c_str(), "pl");
  //  leg1->AddEntry(tp_Pad2[0], fluenceSet1[1].c_str(), "pl");
  leg1->AddEntry(tp_Pad3[0], fluenceSet1[2].c_str(), "pl");
  leg1->AddEntry(tp_Pad4[0], fluenceSet1[3].c_str(), "pl");
  leg1->AddEntry(tp_Pad5[0], fluenceSet1[4].c_str(), "pl");
  leg1->AddEntry(tp_Pad6[0], fluenceSet1[5].c_str(), "pl");

  TLegend *leg2 = new TLegend(0.6,0.65,0.8,0.90,NULL,"brNDC");
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.04);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);
  leg2->SetShadowColor(kWhite);
  leg2->SetHeader("Si N-on-P 200#mum");
  leg2->AddEntry(tp_Pad1[1], fluenceSet2[0].c_str(), "pl");
  //  leg2->AddEntry(tp_Pad2[1], fluenceSet2[1].c_str(), "pl");
  leg2->AddEntry(tp_Pad3[1], fluenceSet2[2].c_str(), "pl");
  leg2->AddEntry(tp_Pad4[1], fluenceSet2[3].c_str(), "pl");
  leg2->AddEntry(tp_Pad5[1], fluenceSet2[4].c_str(), "pl");
  leg2->AddEntry(tp_Pad6[1], fluenceSet2[5].c_str(), "pl");

  TLegend *leg3 = new TLegend(0.6,0.65,0.8,0.90,NULL,"brNDC");
  leg3->SetTextFont(42);
  leg3->SetTextSize(0.04);
  leg3->SetFillColor(kWhite);
  leg3->SetLineColor(kWhite);
  leg3->SetShadowColor(kWhite);
  leg3->SetHeader("Si N-on-P 320#mum");
  leg3->AddEntry(tp_Pad1[2], fluenceSet3[0].c_str(), "pl");
  //  leg3->AddEntry(tp_Pad2[2], fluenceSet3[1].c_str(), "pl");
  leg3->AddEntry(tp_Pad3[2], fluenceSet3[2].c_str(), "pl");
  leg3->AddEntry(tp_Pad4[2], fluenceSet3[3].c_str(), "pl");
  leg3->AddEntry(tp_Pad5[2], fluenceSet3[4].c_str(), "pl");
  leg3->AddEntry(tp_Pad6[2], fluenceSet3[5].c_str(), "pl");

  std::string folder = "plots_ProMediAbove200_Single";

  TCanvas* ch_100 = new TCanvas();
  ch_100->cd();
  tp_Pad1[0]->GetXaxis()->SetRangeUser(-4, +5);
  tp_Pad1[0]->GetXaxis()->SetTitle("t - t_{max}^{0.00E+00} (ns)");
  tp_Pad1[0]->GetYaxis()->SetTitle("amp / amp_{max} ");
  tp_Pad1[0]->Draw("h");
  //  tp_Pad2[0]->Draw("h,same");
  tp_Pad3[0]->Draw("h,same");
  tp_Pad4[0]->Draw("h,same");
  tp_Pad5[0]->Draw("h,same");
  tp_Pad6[0]->Draw("h,same");
  leg1->Draw("same");
  ch_100->Print((folder+"/Si_N120um.png").c_str(), "png");

  TCanvas* ch_200 = new TCanvas();
  ch_200->cd();
  tp_Pad1[1]->GetXaxis()->SetRangeUser(-4, +5);
  tp_Pad1[1]->GetXaxis()->SetTitle("t - t_{max}^{0.00E+00} (ns)");
  tp_Pad1[1]->GetYaxis()->SetTitle("amp / amp_{max} ");
  tp_Pad1[1]->Draw("h");
  //  tp_Pad2[1]->Draw("h,same");
  tp_Pad3[1]->Draw("h,same");
  tp_Pad4[1]->Draw("h,same");
  tp_Pad5[1]->Draw("h,same");
  tp_Pad6[1]->Draw("h,same");
  leg2->Draw("same");
  ch_200->Print((folder+"/Si_N200um.png").c_str(), "png");


  TCanvas* ch_300 = new TCanvas();
  ch_300->cd();
  tp_Pad1[2]->GetXaxis()->SetRangeUser(-4, +5);
  tp_Pad1[2]->GetXaxis()->SetTitle("t - t_{max}^{0.00E+00} (ns)");
  tp_Pad1[2]->GetYaxis()->SetTitle("amp / amp_{max} ");
  tp_Pad1[2]->Draw("h");
  //tp_Pad2[2]->Draw("h,same");
  tp_Pad3[2]->Draw("h,same");
  tp_Pad4[2]->Draw("h,same");
  tp_Pad5[2]->Draw("h,same");
  tp_Pad6[2]->Draw("h,same");
  leg3->Draw("same");
  ch_300->Print((folder+"/Si_N320um.png").c_str(), "png");

  /*
  tp_Pad1[1]->SetLineColor(kGreen+2);
  tp_Pad1[1]->SetMarkerColor(kGreen+2);
  tp_Pad1[2]->SetLineColor(kBlue);
  tp_Pad1[2]->SetMarkerColor(kBlue);

  tp_Pad1[0]->SetLineWidth(2);
  tp_Pad1[1]->SetLineWidth(2);
  tp_Pad1[2]->SetLineWidth(2);

  TLegend *leg2 = new TLegend(0.7,0.70,1.,0.9,NULL,"brNDC");
  leg2->SetTextFont(42);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);
  leg2->SetShadowColor(kWhite);
  leg2->AddEntry(tp_Pad1[0], "Pad1 120#mum", "pl");
  leg2->AddEntry(tp_Pad1[1], "Pad1 200#mum", "pl");
  leg2->AddEntry(tp_Pad1[2], "Pad1 320#mum", "pl");

  TCanvas* ch_pad1 = new TCanvas();
  ch_pad1->cd();
  tp_Pad1[0]->GetXaxis()->SetTitle("time wrt max (ns)");
  tp_Pad1[0]->GetYaxis()->SetTitle("normalized amplitude");
  tp_Pad1[0]->Draw("h");
  tp_Pad1[1]->Draw("h,same");
  tp_Pad1[2]->Draw("h,same");
  leg2->Draw("same");
  ch_300->Print((folder+"/h_Pad1.png").c_str(), "png");
  */
  return 0;
}
