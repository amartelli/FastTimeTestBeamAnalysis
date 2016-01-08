#include "TSystem.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TLine.h"
#include "TStyle.h"
#include "TString.h"
#include "TLatex.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TRandom.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TLine.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include "makeRunSummary.C"

using namespace std;

Int_t ifitType=1;
TString referenceWorkspaceUrl("~/www/HGCal/FastTimeTB/test/Run0/workspace.root");
TString referenceCalibration("~/www/HGCal/FastTimeTB/test/runsummary.root");
 
//
void produceCalibratedEnergySpectrumFor(std::vector<Int_t> runs,Float_t siWidth=120,TString inDir="/store/cmst3/group/hgcal/TimingTB_H2_Jul2015/RECO/6000b73")
{
  setROOTstyle();

  //readout the adc2mip calibration
  TFile *fIn=TFile::Open(referenceCalibration);
  Float_t adc2mip[2]={1.0,1.0};
  for(Int_t siChargeEstIdx=0; siChargeEstIdx<2; siChargeEstIdx++)
    {
      TF1 *adc2mipFunc=(TF1 *)fIn->Get( Form("adc2mip_chargeEst%d_siPad-1_fit%d",siChargeEstIdx,ifitType) );
      adc2mip[siChargeEstIdx]=adc2mipFunc->Eval(siWidth);
    }
  fIn->Close();

  //readout the workspace to get pedestals, noise and fiducial beamspot region definitions
  fIn=TFile::Open(referenceWorkspaceUrl);
  RooWorkspace *w=(RooWorkspace *) fIn->Get("w");
  
  //beam spot 95%CL coordinates reconstructed from wire chambers
  Float_t beamSpotX005_wc=w->var("x005_inc_1")->getVal();
  Float_t beamSpotX095_wc=w->var("x095_inc_1")->getVal();
  Float_t beamSpotY005_wc=w->var("y005_inc_1")->getVal();
  Float_t beamSpotY095_wc=w->var("y095_inc_1")->getVal();

  //beam spot 95%CL coordinates reconstructed from wire chambers requiring Si signals above noise
  Float_t beamSpotX005_siFid=w->var(Form("x005_fid_%s_1",siChargeEst[0].Data()))->getVal();
  Float_t beamSpotX095_siFid=w->var(Form("x095_fid_%s_1",siChargeEst[0].Data()))->getVal();
  Float_t beamSpotY005_siFid=w->var(Form("y005_fid_%s_1",siChargeEst[0].Data()))->getVal();
  Float_t beamSpotY095_siFid=w->var(Form("y095_fid_%s_1",siChargeEst[0].Data()))->getVal();

  //MCP data
  Float_t mcp_pedestal=w->var("wave_max_pedestal_1")->getVal();
  Float_t mcp_noise=w->var("wave_max_noise_1")->getVal();

  //Si data
  std::map<Int_t,std::map<Int_t,Float_t> > siNoise,siPedestal;
  std::map<Int_t,Float_t> templateCalibMap;
  for(Int_t ich=2; ich<4; ich++)
    {
      siNoise[ich]=templateCalibMap;
      siPedestal[ich]=templateCalibMap;
      for(Int_t siChargeEstIdx=0; siChargeEstIdx<2; siChargeEstIdx++)
	{
	  //get the noise estimate for this pad, estimator and selection type
	  TString noiseName(Form("%s_noise_%d",siChargeEst[siChargeEstIdx].Data(),ich));
	  if(ifitType==2)  { noiseName+= "_bias0.0"; }
	  if(ifitType==3)  { noiseName+= "_bias1.0"; }
	  if(ifitType==4)  { noiseName+= "_bias1.0"; }
	  RooRealVar *noiseVar=(RooRealVar *)w->var(noiseName);
	  siNoise[ich][siChargeEstIdx]=noiseVar->getVal();
	  
	  //get the pedestal
	  TString pedestalName(noiseName); pedestalName.ReplaceAll("noise","pedestal");
	  RooRealVar *pedestalVar=(RooRealVar *)w->var(pedestalName);
	  siPedestal[ich][siChargeEstIdx]=pedestalVar->getVal();
	}
    }
  fIn->Close();

  //report reco conditions
  cout << "***************************************************************" << endl
       << "[Reconstruction conditions]" << endl
       << "Beamspot fiducial\t"
       << beamSpotX005_wc << " <x<" << beamSpotX095_wc << "\t"
       << beamSpotY005_wc << " <y<" << beamSpotY095_wc
       << endl
       << "Beamspot Si fiducial\t"
       << beamSpotX005_siFid << " <x<" << beamSpotX095_siFid << "\t"
       << beamSpotY005_siFid << " <y<" << beamSpotY095_siFid
       << endl
       << "MCP\t pedestal=" << mcp_pedestal << " noise=" << mcp_noise << std::endl; 
  for(Int_t siChargeEstIdx=0; siChargeEstIdx<2; siChargeEstIdx++)
    { 
      cout << "Calibration for charge estimator for " << siChargeEst[siChargeEstIdx].Data() << endl;
      for(Int_t ich=2; ich<4; ich++)
	std::cout << "\tSi #" << ich << "\t pedestal=" << siPedestal[ich][siChargeEstIdx] << " noise=" << siNoise[ich][siChargeEstIdx] << std::endl;
      cout << "ADC2MIP conversion=" << adc2mip[siChargeEstIdx] << endl;
    }
  cout << "***************************************************************" << endl;

  //prepare to read the tree
  TChain *H4treeReco=new TChain("H4treeReco");
  for(size_t irun=0; irun<runs.size(); irun++)
    {
      TString url( Form("root://eoscms//eos/cms/%s/RECO_%d.root",inDir.Data(),runs[irun]) );
      H4treeReco->Add(url);
    }

  UInt_t runNumber,spillNumber,evtNumber;
  H4treeReco->SetBranchAddress("runNumber",    &runNumber);
  H4treeReco->SetBranchAddress("spillNumber",  &spillNumber);
  H4treeReco->SetBranchAddress("evtNumber",    &evtNumber);
  Float_t wc_recox[16], wc_recoy[16];
  UInt_t wc_xl_hits[16], wc_xr_hits[16], wc_yu_hits[16], wc_yd_hits[16]; 
  UInt_t nwc;
  H4treeReco->SetBranchAddress("nwc",       &nwc);
  H4treeReco->SetBranchAddress("wc_recox",   wc_recox);
  H4treeReco->SetBranchAddress("wc_recoy",   wc_recoy);
  H4treeReco->SetBranchAddress("wc_xl_hits",   wc_xl_hits);
  H4treeReco->SetBranchAddress("wc_xr_hits",   wc_xr_hits);
  H4treeReco->SetBranchAddress("wc_yu_hits",   wc_yu_hits);
  H4treeReco->SetBranchAddress("wc_yd_hits",   wc_yd_hits);
 
  UInt_t maxch;
  Float_t wave_max[100],charge_integ_smallw_mcp[100], wave_fit_smallw_ampl[100];
  H4treeReco->SetBranchAddress("maxch",                   &maxch);
  H4treeReco->SetBranchAddress("wave_max",                 wave_max);
  H4treeReco->SetBranchAddress("charge_integ_smallw_mcp",  charge_integ_smallw_mcp);
  H4treeReco->SetBranchAddress("wave_fit_smallw_ampl",     wave_fit_smallw_ampl);
  Float_t t_max[100],  t_max_frac30[100], t_max_frac50[100];
  H4treeReco->SetBranchAddress("t_max",                t_max);
  H4treeReco->SetBranchAddress("t_max_frac30",         t_max_frac30);
  H4treeReco->SetBranchAddress("t_max_frac50",         t_max_frac50);
  Float_t wave_aroundmax[100][50], time_aroundmax[100][50];
  H4treeReco->SetBranchAddress("wave_aroundmax",       wave_aroundmax);
  H4treeReco->SetBranchAddress("time_aroundmax",       time_aroundmax);
 
  //prepare output summary tree
  TFile *fOut=TFile::Open(Form("timesummary_%d.root",Int_t(siWidth)),"RECREATE");
  fOut->cd();
  TTree *ttree=new TTree("ttree","ttree");
  ttree->SetDirectory(fOut);
  Float_t calibEn[3],rawEn[3],AoverN[3];
  Float_t out_t_max[3],  out_t_max_frac30[3], out_t_max_frac50[3],out_t_mip10[3],out_t_mip20[3];
  ttree->Branch("rawEn",             rawEn,            "rawEn[3]/F");
  ttree->Branch("AoverN",            AoverN,           "AoverN[3]/F");
  ttree->Branch("calibEn",           calibEn,          "calibEn[3]/F");
  ttree->Branch("t_max",             out_t_max,            "t_max[3]/F");
  ttree->Branch("t_max_frac30",      out_t_max_frac30,     "t_max_frac30[3]/F");
  ttree->Branch("t_max_frac50",      out_t_max_frac50,     "t_max_frac50[3]/F");
  ttree->Branch("t_mip10",           out_t_mip10,          "t_mip10[3]/F");
  ttree->Branch("t_mip20",           out_t_mip20,          "t_mip20[3]/F");
  ttree->Branch("wave_aroundmax",    wave_aroundmax,   "wave_aroundmax[3][50]/F");
  ttree->Branch("time_aroundmax",    time_aroundmax,   "time_aroundmax[3][50]/F");

  //canvas to hold wave forms
  TCanvas *c=new TCanvas("c","c",500,500);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.1);

  //loop over events
  Int_t nSel(0),nWaveFormPlots(0);
  for(int i=0; i<H4treeReco->GetEntries(); i++)
    {
      H4treeReco->GetEntry(i);
      
      //require hits in the wire chambers
      Bool_t isEmpty(wc_xl_hits[0]==0 && wc_xr_hits[0]==0 && 
		     wc_xl_hits[1]==0 && wc_xr_hits[1]==0 && 
		     wc_yu_hits[1]==0 && wc_yd_hits[1]==0);
      if(isEmpty) continue;
      
      //require reconstructed beamspot to be meaningful
      Bool_t isInclusiveBeamspot(wc_recox[1]>beamSpotX005_wc && 
				 wc_recox[1]<beamSpotX095_wc && 
				 wc_recoy[1]>beamSpotY005_wc && 
				 wc_recoy[1]<beamSpotY095_wc);
      if(!isInclusiveBeamspot) continue;
      
      //fiducial in Si
      Bool_t isFiducialBeamspot(wc_recox[1]>beamSpotX005_siFid && 
				wc_recox[1]<beamSpotX095_siFid &&
				wc_recoy[1]>beamSpotY005_siFid && 
				wc_recoy[1]<beamSpotY095_siFid);
      //if(!isFiducialBeamspot) continue;

      //has MCP trigger
      Bool_t hasMCPTrigger(wave_max[1]-mcp_pedestal>50*mcp_noise);
      if(!hasMCPTrigger) continue;

      nSel++;
      
      //loop over channels
      Bool_t doWaveFormPlot(false);
      std::vector<TGraph *> waveFormGrs;
      for(int ich=1; ich<4; ich++)
	{
	  int chIdx(ich==1 ? 2 : ich-2);

	  //save pre-computed time estimates
	  out_t_max[chIdx]            = t_max[ich];
	  out_t_max_frac30[chIdx]     = t_max_frac30[ich];
	  out_t_max_frac50[chIdx]     = t_max_frac50[ich];

	  //get charge and calibrate it to the energy
	  rawEn[chIdx]=TMath::Max(Float_t(0.),Float_t(charge_integ_smallw_mcp[ich]-siPedestal[ich][0]));
	  if(chIdx<2)
	    {
	      AoverN[chIdx]  = rawEn[chIdx]/siNoise[ich][0];
	      calibEn[chIdx] = rawEn[chIdx]/adc2mip[0];
	    }

	  //save waveform
	  if(ich==1) continue;
	  if(ich==2 && calibEn[chIdx]>5 && nWaveFormPlots<30 && gRandom->Uniform()>0.5 ) doWaveFormPlot=true;
	  if(!doWaveFormPlot) continue;

	  TGraph *timeGr=new TGraph();
	  timeGr->SetFillStyle(0);
	  timeGr->SetMarkerStyle(20+ich);
	  timeGr->SetName( Form("ch%d",ich) );
	  timeGr->SetTitle( Form("#scale[0.8]{E(%d) = %3.2f MIP}",ich,calibEn[chIdx]) );
	  for(Int_t itime=0; itime<50; itime++)
	    timeGr->SetPoint(timeGr->GetN(),time_aroundmax[ich][itime]*1e9,wave_aroundmax[ich][itime]);
	  waveFormGrs.push_back(timeGr);	  
	}
      
      //save tree
      for(int ich=2; ich<4; ich++)
	{
	  Int_t chIdx=ich-2;
	  out_t_mip10[chIdx]=999999.;
	  out_t_mip20[chIdx]=999999.;
	  for(Int_t itime=0; itime<50; itime++)
	    {
	      Float_t dt=t_max[ich]-time_aroundmax[ich][itime]*1e9;
	      Float_t ampl=(wave_aroundmax[ich][itime]-siPedestal[ich][1])/adc2mip[1];
	      if(dt<0)
		{
		  if(ampl>10)  out_t_mip10[chIdx] = TMath::Min(time_aroundmax[ich][itime],out_t_mip10[chIdx]);
		  if(ampl>20)  out_t_mip20[chIdx] = TMath::Min(time_aroundmax[ich][itime],out_t_mip20[chIdx]);
		}

	      //update with calibrated values
	      time_aroundmax[chIdx][itime]=dt;	      
	      wave_aroundmax[chIdx][itime]=ampl;
	    }
	  out_t_mip10[chIdx] = out_t_mip10[chIdx]>1e3 ? -1 : out_t_mip10[chIdx]*1e9;
	  out_t_mip20[chIdx] = out_t_mip20[chIdx]>1e3 ? -1 : out_t_mip20[chIdx]*1e9;
	}
      ttree->Fill();

      //save graph if required
      if(!doWaveFormPlot) continue;

      nWaveFormPlots++;
      c->Clear();
      TLegend *leg=new TLegend(0.15,0.85,0.4,0.75);
      leg->SetTextFont(42);
      leg->SetTextSize(0.035);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      
      Float_t ymin(9999999999.),ymax(-ymin);
      for(size_t igr=0; igr<waveFormGrs.size(); igr++)
	{
	  waveFormGrs[igr]->SetLineColor(igr+1);
	  waveFormGrs[igr]->SetMarkerColor(igr+1);
	  waveFormGrs[igr]->SetFillStyle(0);
	  waveFormGrs[igr]->Draw(igr==0? "ap" : "p");
	  waveFormGrs[igr]->GetXaxis()->SetTitle("Time [ns]");
	  waveFormGrs[igr]->GetYaxis()->SetTitle("ADC counts");
	  waveFormGrs[igr]->GetYaxis()->SetTitleOffset(1.4);
	  leg->AddEntry(waveFormGrs[igr],waveFormGrs[igr]->GetTitle(),"p");
	  if(igr==0){
	    ymin=TMath::Min((Float_t)waveFormGrs[igr]->GetYaxis()->GetXmin(),(Float_t)ymin);
	    ymax=TMath::Max((Float_t)waveFormGrs[igr]->GetYaxis()->GetXmax(),(Float_t)ymax);
	  }
	}

      for(size_t igr=0; igr<waveFormGrs.size(); igr++)
	{
	  TLine *l=new TLine;
	  l->SetLineColor(igr+1);
	  l->SetLineStyle(2);
	  l->DrawLine(out_t_max[igr],ymin,out_t_max[igr],ymax);
	  l->DrawLine(out_t_max_frac30[igr],ymin,out_t_max_frac30[igr],ymax);
	  //l->DrawLine(ymin,t_max_frac50[igr],ymax,t_max_frac50[igr]);	  
	}
      
      TLatex *txt=new TLatex;
      txt->SetTextFont(42);
      txt->SetTextSize(0.035);
      txt->SetNDC();
      txt->DrawLatex(0.15,0.90,"#bf{CMS}");
      txt->DrawLatex(0.15,0.85,"#scale[0.8]{#it{Fast-timing test beam preliminary}}");
      txt->DrawLatex(0.6,0.90,Form("#scale[0.8]{Run %d Spill %d Event %d}",runNumber,spillNumber,evtNumber));
      
      leg->Draw();
      
      c->SaveAs(Form("waveform_%d_si%d.png",nWaveFormPlots,(int)siWidth));
    }

  cout << nSel << " events selected out of " << H4treeReco->GetEntries() << " initial events" << endl;

  fOut->cd();
  ttree->Write();
  fOut->Close();
}

//
void getCalibratedEnergySpectrum(Float_t siWidth=285,TString inDir="/store/cmst3/group/hgcal/TimingTB_H2_Jul2015/RECO/6000b73")
{
  std::vector<Int_t> runs;
  if(siWidth==285)
    {

      referenceWorkspaceUrl="~/www/HGCal/FastTimeTB/test/Run3346/workspace.root";

      //3X0
      runs.push_back(3330);
      runs.push_back(3333);
      runs.push_back(3349);
      runs.push_back(3350);
      runs.push_back(3352);
      runs.push_back(3358);

      //2X0
      runs.push_back(3334);
      runs.push_back(3335);
      runs.push_back(3336);
      runs.push_back(3337);
      runs.push_back(3338);
      runs.push_back(3339);
      runs.push_back(3340);
      runs.push_back(3341);
      runs.push_back(3342);
      runs.push_back(3348);

      //4X0
      runs.push_back(3353);
      runs.push_back(3355);
      runs.push_back(3356);
    }
  if(siWidth==133)
    {
      referenceWorkspaceUrl="~/www/HGCal/FastTimeTB/test/Run3363/workspace.root";

      //4X0
      runs.push_back(3358);
      runs.push_back(3362);
    }
  if(siWidth==211.5)
    {
      //4X0
      referenceWorkspaceUrl="~/www/HGCal/FastTimeTB/test/Run0/workspace.root";
      runs.push_back(3365);
      runs.push_back(3367);
    }



  produceCalibratedEnergySpectrumFor(runs,siWidth,inDir);
};
