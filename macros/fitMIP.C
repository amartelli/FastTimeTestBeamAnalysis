#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TString.h"

#include "fitFunctions.C"

void fitMIP(int ch, int run, float xc, float yc,TString tag)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(111111);
  TFile *_file0 = TFile::Open(Form("root://eoscms//eos/cms/store/caf/user/meridian/test/RECO_%d.root",run));
  TTree* H4treeReco=(TTree*)_file0->Get("H4treeReco");
  TCanvas* c1=new TCanvas("c1","c1",900,900);
  H4treeReco->Draw(Form("wave_fit_smallw_ampl[%d]>>pedestal(250,-50,50)",ch),"wc_xl_hits[0]==0 && wc_xr_hits[0]==0 && wc_xl_hits[1]==0 && wc_xr_hits[1]==0 && wc_yu_hits[1]==0 && wc_yd_hits[1]==0");
  TH1F* pedestal=(TH1F*)gROOT->FindObject("pedestal");
  pedestal->Fit("gaus");
  c1->SaveAs(Form("noise_ch%d_run%d.png",ch,run));
  H4treeReco->Draw(Form("wave_fit_smallw_ampl[%d]-%4.1f>>noise_corrected(250,-50,50)",ch,pedestal->GetFunction("gaus")->GetParameter(1)),"wc_xl_hits[0]==0 && wc_xr_hits[0]==0 && wc_xl_hits[1]==0 && wc_xr_hits[1]==0 && wc_yu_hits[1]==0 && wc_yd_hits[1]==0");
  TH1F* noise_corrected=(TH1F*)gROOT->FindObject("noise_corrected");
  noise_corrected->Fit("gaus");
  c1->SaveAs(Form("noise_corrected_ch%d_run%d.png",ch,run));

  int cha;
  float mcpdelta;

  if (ch==2)
    {
      cha=3;
      mcpdelta=5.1;
    }
  else if (ch==3)
    {
      cha=2;
      mcpdelta=5.6;
    }
  else
    {
      cha=-999;
      mcpdelta=-999;
    }
    
  H4treeReco->Draw(Form("wave_fit_smallw_ampl[%d]-%4.1f>>signal1(50,-50,200)",ch,pedestal->GetFunction("gaus")->GetParameter(1)),Form("wave_max[1]>100 && abs(wc_recox[1]-%3.1f)<2.5 && abs(wc_recoy[1]-%3.1f)<2.5 && abs(wc_recox[0]-wc_recox[1]-24.22)<3 && wave_fit_smallw_ampl[%d]<80",xc,yc,cha));
  //  H4treeReco->Draw(Form("wave_fit_smallw_ampl[%d]-%4.1f>>signal1(100,-50,200)",ch,pedestal->GetFunction("gaus")->GetParameter(1)),Form("wave_max[%d]>10 && wave_max[%d]<70 && abs(t_max[%d]-t_max[1]+%2.1f)<2. && wave_max[1]>100",yc,xc,cha,cha,cha,mcpdelta));
  // H4treeReco->Draw(Form("wave_fit_smallw_ampl[3]-%4.1f>>signal2(100,-50,200)",pedestal->GetFunction("gaus")->GetParameter(1)),"wave_max[1]>200 && abs(t_max[2]-t_max[1]+5.6)<1. && wave_max[2]>10 && wave_max[2]<70"); 
  // H4treeReco->Draw(Form("wave_fit_smallw_ampl[2]-%4.1f>>signal1(80,-200,1000)",pedestal->GetFunction("gaus")->GetParameter(1)),"wave_max[1]>200 && abs(t_max[3]-t_max[1]+5.7)<.5 && wave_max[3]>10 && wave_max[3]<200 && abs(wc_recox[1]-wc_recox[0]+24.2)<4 && abs(((wc_recox[1]+(wc_recox[0]-24.2)/2)-12.5)-17)<10 && abs((wc_recoy[1]+35)+1)<10"); 
  // H4treeReco->Draw(Form("wave_fit_smallw_ampl[3]-%4.1f>>signal2(80,-200,1000)",pedestal->GetFunction("gaus")->GetParameter(1)),"wave_max[1]>200 && abs(t_max[2]-t_max[1]+5.7)<.5 && wave_max[2]>10 && wave_max[2]<200 && abs(wc_recox[1]-wc_recox[0]+24.2)<4 && abs(((wc_recox[1]+(wc_recox[0]-24.2)/2)-12.5)-17)<10 && abs((wc_recoy[1]+35)+1)<10"); 
  // H4treeReco->Draw(Form("wave_fit_smallw_ampl[2]-%4.1f>>signal1(250,-1000,5000)",pedestal->GetFunction("gaus")->GetParameter(1)),"wave_max[1]<100 && abs(t_max[3]-t_max[1]+5.7)>4 && wave_max[3]<100"); 
  // H4treeReco->Draw(Form("wave_fit_smallw_ampl[3]-%4.1f>>signal2(250,-1000,5000)",pedestal->GetFunction("gaus")->GetParameter(1)),"wave_max[1]<100 && abs(t_max[2]-t_max[1]+5.7)>4 && wave_max[2]<100"); 
  TH1F* signal1=(TH1F*)gROOT->FindObject("signal1");
  // TH1F* signal2=(TH1F*)gROOT->FindObject("signal2");
  TH1F* signal=(TH1F*)signal1->Clone("signal");
  signal->Sumw2();
  // signal->Add(signal1);
  //H4treeReco->Draw(Form("wave_fit_smallw_ampl[3]-%4.1f>>signal(250,-500,5000)",pedestal->GetFunction("gaus")->GetParameter(1)),"wave_max[1]>200"); 
  //  gROOT->ProcessLine(".L ~/work/H2_2015/FastTimeTestBeamAnalysis/macros/fitFunctions.C++");
  TF1* fitFunc=new TF1("fitFunc",sigFunc,-50,200,7);
  fitFunc->FixParameter(1,0.);
  fitFunc->FixParameter(2,pedestal->GetFunction("gaus")->GetParameter(2));
  fitFunc->SetParLimits(0,0,100000);
  fitFunc->SetParLimits(3,1,10);
  fitFunc->SetParName(0,"NoiseNorm");
  fitFunc->SetParName(3,"LandauWidth");
  fitFunc->SetParName(4,"LandauMPV");
  fitFunc->SetParName(5,"LandauNorm");
  fitFunc->SetParLimits(4,10,100);
  //  fitFunc->SetParLimits(6,50,3000);
  fitFunc->SetParameter(3,4);
  fitFunc->SetParameter(4,40);
  fitFunc->FixParameter(6,pedestal->GetFunction("gaus")->GetParameter(2));
  // fitFunc->FixParameter(3,100.);
  //  fitFunc->FixParameter(4,230.);
  // fitFunc->FixParameter(5,0.);
  signal->Fit("fitFunc","BR");
  signal->SetMarkerStyle(20);
  signal->SetMarkerSize(0.8);
  signal->GetXaxis()->SetTitle("Signal Amplitude [ADC]");
  signal->Draw("PE");
  fitFunc->SetLineWidth(3);
  fitFunc->Draw("SAME");
  Double_t fitparams[7];
  fitFunc->GetParameters(fitparams);
  TF1 *landau= new TF1("landau",lanconvgau,-50,200,4);
  landau->SetParameters(&fitparams[3]);
  landau->SetLineColor(kBlue);
  landau->SetLineStyle(2);
  landau->Draw("SAME");

  TF1 *noi=new TF1("noi","gaus",-50,200);
  noi->SetParameter(0,fitparams[0]);
  noi->SetParameter(1,fitparams[1]);
  noi->SetParameter(2,fitparams[2]);
  noi->SetLineColor(kMagenta);
  noi->SetLineStyle(2);
  noi->Draw("SAME");
  c1->SaveAs(Form("MIP_ch%d_run%d_%s.png",ch,run,tag.Data())); 

}
