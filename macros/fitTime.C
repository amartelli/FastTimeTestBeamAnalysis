#include "TFile.h"
#include "TChain.h"
#include "TROOT.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMath.h"
#include "TString.h"

void fitTime(float deltaT1,float deltaT2, float MIPcalib, float MIPcut, TString tag, TString title)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111111);
  TCanvas* c1=new TCanvas("c1","c1",900,700);
  TChain *c=new TChain("H4treeReco");
  // c->Add("root://eoscms///eos/cms/store/caf/user/meridian/test/RECO_3353.root");
  // c->Add("root://eoscms///eos/cms/store/caf/user/meridian/test/RECO_3354.root");
  // c->Add("root://eoscms///eos/cms/store/caf/user/meridian/test/RECO_3355.root");
  // c->Add("root://eoscms///eos/cms/store/caf/user/meridian/test/RECO_3356.root");
  c->Add("root://eoscms///eos/cms/store/caf/user/meridian/test/RECO_3367.root");

  c->Draw(Form("(t_max_frac50[2]-t_max_frac50[1]+%4.3f)>>h1(100,-0.2,0.2)",deltaT1),Form("wave_max[1]>200 && wave_fit_smallw_ampl[2]>%3.1f*%3.1f",MIPcalib,MIPcut));
  c->Draw(Form("(t_max_frac50[3]-t_max_frac50[1]+%4.3f)>>h2(100,-0.2,0.2)",deltaT2),Form("wave_max[1]>200 && wave_fit_smallw_ampl[3]>%3.1f*%3.1f",MIPcalib,MIPcut));
  // c->Draw("(t_max_frac50[3]-t_max_frac50[1]+5.41)>>h2(100,-0.2,0.2)","wave_max[1]>200 && wave_fit_smallw_ampl[3]>41*20");
  // c->Draw("(t_max_frac50[2]-t_max_frac50[3]+0.155)>>h3(100,-0.2,0.2)","wave_fit_smallw_ampl[2]>41*20 && wave_fit_smallw_ampl[3]>41*20");
  c->Draw(Form("(t_max_frac50[2]-t_max_frac50[3]+%4.3f)>>h3(100,-0.2,0.2)",deltaT1-deltaT2),Form("wave_max[1]>200 && wave_fit_smallw_ampl[3]>%3.1f*%3.1f && wave_fit_smallw_ampl[2]>%3.1f*%3.1f",MIPcalib,MIPcut,MIPcalib,MIPcut));

  TLatex t;
  t.SetTextSize(0.04);
  TH1F* h1=(TH1F*)gROOT->FindObject("h1");
  h1->Sumw2();
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.3);
  h1->SetMarkerColor(kMagenta);
  h1->SetLineColor(kBlack);  
  h1->Fit("gaus");
  h1->GetXaxis()->SetTitle("#Delta_{t} Si_{1} - MCP (ns)");
  t.DrawLatexNDC(0.07,0.92,title);
  t.DrawLatexNDC(0.16,0.8,Form("Amplitude > %2.0f MIP",MIPcut));
  c1->SaveAs(Form("DeltaT_ch%d_MCP_%s.png",2,tag.Data()));

  TH1F* h2=(TH1F*)gROOT->FindObject("h2");
  h2->Sumw2();
  h2->SetMarkerStyle(20);
  h2->SetMarkerSize(0.3);
  h2->SetMarkerColor(kMagenta);
  h2->GetXaxis()->SetTitle("#Delta_{t} Si_{2} - MCP (ns)");
  h2->SetLineColor(kBlack);  
  h2->Fit("gaus");
  t.DrawLatexNDC(0.07,0.92,title);
  t.DrawLatexNDC(0.16,0.8,Form("Amplitude > %2.0f MIP",MIPcut));
  c1->SaveAs(Form("DeltaT_ch%d_MCP_%s.png",3,tag.Data()));

  TH1F* h3=(TH1F*)gROOT->FindObject("h3");
  h3->Sumw2();
  h3->SetMarkerStyle(20);
  h3->SetMarkerSize(0.3);
  h3->SetMarkerColor(kMagenta);
  h3->GetXaxis()->SetTitle("#Delta_{t} Si_{1} - Si_{2} (ns)");
  h3->SetLineColor(kBlack);  
  h3->Fit("gaus");
  t.DrawLatexNDC(0.07,0.92,"e- 50 GeV Si 320#mum p-type. 4 X_{0} lead");
  t.DrawLatexNDC(0.16,0.8,"Amplitude > 20 MIP");
  t.DrawLatexNDC(0.07,0.92,title);
  t.DrawLatexNDC(0.16,0.8,Form("Amplitude > %2.0f MIP",MIPcut));
  c1->SaveAs(Form("DeltaT_ch%d_ch%d_%s.png",2,3,tag.Data()));

  TF1* fRes=new TF1("fRes","TMath::Sqrt(([0]*[0])/(x*x)+([1]*[1]))",1,100);
  fRes->SetLineWidth(2);
  
  // c->Draw("(t_max_frac50[2]-t_max_frac50[1]+5.57):wave_fit_smallw_ampl[2]/41.>>h1_2(30,1,41,50,-0.2,0.2)","wave_max[1]>200 && wave_fit_smallw_ampl[2]>20");
  c->Draw(Form("(t_max_frac50[2]-t_max_frac50[1]+%4.3f):wave_fit_smallw_ampl[2]/%3.1f>>h1_2(30,1,41,50,-0.2,0.2)",deltaT1,MIPcalib),"wave_max[1]>200 && wave_fit_smallw_ampl[2]>20");
  TH2F* h1_2=(TH2F*)gROOT->FindObject("h1_2");
  h1_2->FitSlicesY();
  TH1F* h1_2_2=(TH1F*)gROOT->FindObject("h1_2_2");
  TH1F* h1_2_2_corr=(TH1F*)h1_2_2->Clone(h1_2_2->GetName()+TString("_corr"));
  for (int i=1;i<=h1_2_2_corr->GetNbinsX();++i)
    {
      if (h1_2_2->GetBinContent(i)>0.0221)
	{
	  h1_2_2_corr->SetBinContent(i,TMath::Sqrt(h1_2_2->GetBinContent(i)*h1_2_2->GetBinContent(i)-0.0221*0.0221));
	  h1_2_2_corr->SetBinError(i,TMath::Sqrt(h1_2_2->GetBinError(i)*h1_2_2->GetBinError(i)+0.003*0.003)); 
	}
      else
	{
	  h1_2_2_corr->SetBinContent(i,0.);
	  h1_2_2_corr->SetBinError(i,0.008);
	}
      //     std::cout << i << "," << h1_2_2_corr->GetBinContent(i) << "," << h1_2_2_corr->GetBinError(i) << std::endl;
    }
  h1_2_2_corr->Draw();
  h1_2_2_corr->Fit("fRes","R");
  h1_2_2_corr->SetMinimum(0.);
  h1_2_2_corr->SetMaximum(0.3);
  h1_2_2_corr->GetXaxis()->SetTitle("Signal amplitude (# MIP)");
  h1_2_2_corr->GetYaxis()->SetTitle("#sigma_{t} (ns)");
  h1_2_2_corr->SetMarkerStyle(20);
  h1_2_2_corr->SetMarkerSize(1.2);
  h1_2_2_corr->SetMarkerColor(kBlue);
  h1_2_2_corr->SetLineColor(kBlack);
  t.DrawLatexNDC(0.12,0.92,title);
  t.DrawLatexNDC(0.12,0.82,"#sigma_{MCP}=22.1ps subtracted");
  c1->SaveAs(Form("DeltaT_vs_Ampl_ch%d_MCP_%s.png",2,tag.Data()));


  c->Draw(Form("(t_max_frac50[3]-t_max_frac50[1]+%4.3f):wave_fit_smallw_ampl[3]/%3.1f>>h2_2(30,1,41,50,-0.2,0.2)",deltaT2,MIPcalib),"wave_max[1]>200 && wave_fit_smallw_ampl[3]>20");
  TH2F* h2_2=(TH2F*)gROOT->FindObject("h2_2");
  h2_2->FitSlicesY();
  TH1F* h2_2_2=(TH1F*)gROOT->FindObject("h2_2_2");
  TH1F* h2_2_2_corr=(TH1F*)h2_2_2->Clone(h2_2_2->GetName()+TString("_corr"));
  for (int i=1;i<=h2_2_2_corr->GetNbinsX();++i)
    {
      if (h2_2_2->GetBinContent(i)>0.0221)
	{
	  h2_2_2_corr->SetBinContent(i,TMath::Sqrt(h2_2_2->GetBinContent(i)*h2_2_2->GetBinContent(i)-0.0221*0.0221));
	  h2_2_2_corr->SetBinError(i,TMath::Sqrt(h2_2_2->GetBinError(i)*h2_2_2->GetBinError(i)+0.003*0.003)); 
	}
      else
	{
	  h2_2_2_corr->SetBinContent(i,0.);
	  h2_2_2_corr->SetBinError(i,0.008);
	}
      //     std::cout << i << "," << h2_2_2_corr->GetBinContent(i) << "," << h2_2_2_corr->GetBinError(i) << std::endl;
    }

  h2_2_2_corr->Draw();
  h2_2_2_corr->Fit("fRes","R");
  h2_2_2_corr->SetMinimum(0.);
  h2_2_2_corr->SetMaximum(0.3);
  h2_2_2_corr->GetXaxis()->SetTitle("Signal amplitude (# MIP)");
  h2_2_2_corr->GetYaxis()->SetTitle("#sigma_{t} (ns)");
  h2_2_2_corr->SetMarkerStyle(20);
  h2_2_2_corr->SetMarkerSize(1.2);
  h2_2_2_corr->SetMarkerColor(kBlue);
  h2_2_2_corr->SetLineColor(kBlack);
  t.DrawLatexNDC(0.12,0.92,title);
  t.DrawLatexNDC(0.12,0.82,"#sigma_{MCP}=22.1ps subtracted");
  c1->SaveAs(Form("DeltaT_vs_Ampl_ch%d_MCP_%s.png",3,tag.Data()));

}

