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
#include "TFitResult.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#ifndef __CINT__
#include "RooCFunction1Binding.h" 
#include "RooCFunction3Binding.h"
#endif
#include "RooTFnBinding.h" 

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include "fitFunctions.C"

using namespace std;

Bool_t doPedestalByQuantile=true;
TString mcpChargeEst="wave_max";
TString siChargeEst[]={"charge_integ_smallw_mcp","wave_fit_smallw_ampl"};
TString CMSLABEL="#splitline{#bf{CMS} #it{work in progress}}{#scale[0.8]{#it{Fast-timing testbeam}}}";
struct MIPFitSummary_t
{
  Int_t chargeEstimator, siPad, fitType,ndof;
  Float_t chi2,prob,adc2mip,adc2mipUnc,mip1frac,mip1fracUnc,mip2frac,mip2fracUnc,relWidth,relWidthUnc;
};



//
std::pair<Float_t,Float_t> getParameterRatio(float num,float numUnc,float den, float denUnc, float covNumDen)
{
  float val=  num/den;
  float valUnc=fabs(val)*TMath::Sqrt( fabs(pow(numUnc/num,2)  + pow(denUnc/den,2) -2*covNumDen/(num*den)) );
  return std::pair<Float_t,Float_t>(val,valUnc);
}

//
void setROOTstyle()
{
  //configure root
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111111);
  gStyle->SetStatW(0.20);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.9);
  gStyle->SetStatH(0.15);
  gROOT->SetBatch(true);
}

//pedestal control in empty events
void doPedestals(TChain* H4treeReco,RooWorkspace *w, Float_t sigmaBias, TCanvas *c,TString outDir)
{
  Int_t runNumber=w->var("runNumber")->getVal();

  TString emptyEvtCut("wc_xl_hits[0]==0 && wc_xr_hits[0]==0 && wc_xl_hits[1]==0 && wc_xr_hits[1]==0 && wc_yu_hits[1]==0 && wc_yd_hits[1]==0");


  for(Int_t ich=1; ich<4; ich++)
    {
      for(size_t iest=0; iest<sizeof(siChargeEst)/sizeof(TString); iest++)
	{
	  
	  TString chargeEst(siChargeEst[iest]);
	  if(ich==1) 
	    {
	      chargeEst=mcpChargeEst;
	      if(iest>0) continue;
	    }

	  TString biasPostFix("");
	  TString biasCut("");
	  if(ich>1 && sigmaBias>=0)
	    {
	      biasPostFix=Form("_bias%3.1f",sigmaBias);
	      Int_t otherSiPadCh(ich==2? 3 : 2);
	      RooRealVar *otherPedestal=w->var(Form("%s_pedestal_%d",chargeEst.Data(),otherSiPadCh));
	      RooRealVar *otherNoise=w->var(Form("%s_noise_%d",chargeEst.Data(),otherSiPadCh));
	      biasCut = Form(" && %s[%d]-%f>%f*%f",
			     chargeEst.Data(),
			     otherSiPadCh,
			     otherPedestal?otherPedestal->getVal():0.0,
			     sigmaBias,
			     otherNoise?otherNoise->getVal():0.0);
	    }
	  
	  
	  //pedestal evolution
	  TGraphErrors *pedestalEvolGr=new TGraphErrors();
	  pedestalEvolGr->SetName(Form("pedestalevol_%s_%d%s",chargeEst.Data(),ich,biasPostFix.Data()));
	  pedestalEvolGr->SetMarkerStyle(20);
	  pedestalEvolGr->SetFillStyle(0);

	  //project pedestal vs spill to trace the evolution throughout the run, estimate from gaussian fit or quantiles
	  Float_t baseline(0);
	  RooRealVar *baselineVar=w->var(Form("%s_pedestal_%d%s",chargeEst.Data(),ich,biasPostFix.Data()));
	  if(baselineVar) baseline=baselineVar->getVal();
	  H4treeReco->Draw(Form("%s[%d]-%f:spillNumber>>pedestalvsspill(500,0,500,250,-250,250)",
				chargeEst.Data(),
				ich,
				baselineVar? baselineVar->getVal():0.0),
			   emptyEvtCut+biasCut);
	  
	  TH2* pedestalvsspill=(TH2*)gROOT->FindObject("pedestalvsspill");
	  Int_t minSpill(99999),maxSpill(-1);
	  Float_t pedestalAvg(0),pedestalSigma(0);
	  Float_t pedestalAvgUnc(0),pedestalSigmaUnc(0);
	  for(Int_t xbin=0; xbin<=pedestalvsspill->GetNbinsX(); xbin++)
	    {
	      Int_t firstxbin(0),lastxbin(-1);
	      if(xbin>0) {firstxbin=xbin; lastxbin=xbin;}
	      TH1 *proj=pedestalvsspill->ProjectionY("_py",firstxbin,lastxbin);
	      
	      //at least 10 events to do the control
	      if(proj->Integral()<10) continue;
	      proj->Fit("gaus","MQ+");
	      TF1 *gaus=proj->GetFunction("gaus"); 
	      Float_t avg   = gaus->GetParameter(1);
	      Float_t sigma = gaus->GetParameter(2);
	      gaus->SetParLimits(1,avg-sigma,avg+sigma);
	      proj->Fit(gaus,"MQ+");
	      avg   = gaus->GetParameter(1);
	      sigma = gaus->GetParameter(2);

	      if(doPedestalByQuantile)
		{
		  Double_t bsq_x[3]   ={0.16, 0.50, 0.84};
		  Double_t bsq_val[3] ={0.00, 0.00, 0.00};
		  proj->GetQuantiles(3,bsq_val,bsq_x);
		  avg=bsq_val[1];
		  sigma=0.5*(bsq_val[2]-bsq_val[0]);
		}

	      //
	      if(xbin==0) 
		{ 
		  pedestalAvg=avg;                     pedestalSigma=sigma; 
		  pedestalAvgUnc=gaus->GetParError(1); pedestalSigmaUnc=gaus->GetParError(2);
		  
		  c->Clear();
		  if(ich==1) c->SetLogy();
		  else c->SetLogy(false);
		  proj->Draw("PE");
		  proj->GetXaxis()->SetTitle("Integrated charge [ADC]");
		  proj->GetYaxis()->SetTitle("Events");
		  
		  TLatex txt;
		  txt.SetNDC();
		  txt.SetTextFont(42);
		  txt.SetTextSize(0.03);
		  txt.DrawLatex(0.15,0.93,CMSLABEL);
		  txt.DrawLatex(0.8,0.93,Form("Run %d", runNumber)); 
		  if(baselineVar)
		    c->SaveAs(outDir+Form("/residualincnoisefit_%s_ch%d%s.png",chargeEst.Data(),ich,biasPostFix.Data()));
		  else
		    c->SaveAs(outDir+Form("/incnoisefit_%s_ch%d%s.png",chargeEst.Data(),ich,biasPostFix.Data()));
		}
	      else
		{
		  Int_t spillNb=pedestalvsspill->GetXaxis()->GetBinLowEdge(xbin);
		  if(spillNb<minSpill) minSpill=spillNb;
		  if(spillNb>maxSpill) maxSpill=spillNb;
		  
		  Int_t np=pedestalEvolGr->GetN();
		  pedestalEvolGr->SetPoint(np,spillNb,avg);
		  pedestalEvolGr->SetPointError(np,0,sigma);
		}
	      
	      proj->Delete();
	    }
	  
	  TGraph *pedestalGr=new TGraph();
	  pedestalGr->SetName(Form("pedestal_%s_%d%s",chargeEst.Data(),ich,biasPostFix.Data()));
	  pedestalGr->SetFillColor(kGray);
	  pedestalGr->SetFillStyle(1001);
	  pedestalGr->SetLineColor(kGray);
	  pedestalGr->SetMarkerColor(kGray);
	  pedestalGr->SetPoint(0,minSpill,pedestalAvg-pedestalSigma);
	  pedestalGr->SetPoint(1,minSpill,pedestalAvg+pedestalSigma);
	  pedestalGr->SetPoint(2,maxSpill,pedestalAvg+pedestalSigma);
	  pedestalGr->SetPoint(3,maxSpill,pedestalAvg-pedestalSigma);
	  pedestalGr->SetPoint(4,maxSpill,pedestalAvg-pedestalSigma);
     
	  //draw summary
	  c->Clear();
	  c->SetLogy(false);
	  pedestalGr->Draw("af");
	  pedestalGr->GetYaxis()->SetRangeUser(pedestalGr->GetYaxis()->GetXmin()*1.5,pedestalGr->GetYaxis()->GetXmax()*1.5);
	  pedestalGr->GetXaxis()->SetTitle("Spill number");
	  pedestalGr->GetYaxis()->SetTitle("<Noise>");
	  TLine *line = new TLine(minSpill,pedestalAvg,maxSpill,pedestalAvg);
	  line->Draw();
	  pedestalEvolGr->Draw("p");
	  TLatex txt;
	  txt.SetNDC();
	  txt.SetTextFont(42);
	  txt.SetTextSize(0.03);
	  txt.DrawLatex(0.15,0.93,CMSLABEL);
	  txt.DrawLatex(0.8,0.93,Form("Run %d", runNumber)); 
	  if(ich>1 && sigmaBias>=0) txt.DrawLatex(0.5,0.93,Form("Noise bias > %3.1f#sigma",sigmaBias));

	  if(baselineVar)
	    {
	      c->SaveAs(outDir+Form("/residualnoise_%s_ch%d%s.png",chargeEst.Data(),ich,biasPostFix.Data()));
	    }
	  else
	    {
	      c->SaveAs(outDir+Form("/noise_%s_ch%d%s.png",chargeEst.Data(),ich,biasPostFix.Data()));
	      w->factory(Form("%s_pedestal_%d%s[%f]",chargeEst.Data(),ich,biasPostFix.Data(),pedestalAvg));
	      w->var(Form("%s_pedestal_%d%s",chargeEst.Data(),ich,biasPostFix.Data()))->setError(pedestalAvgUnc);
	      w->factory(Form("%s_noise_%d%s[%f]",chargeEst.Data(),ich,biasPostFix.Data(),pedestalSigma));
	      w->var(Form("%s_noise_%d%s",chargeEst.Data(),ich,biasPostFix.Data()))->setError(pedestalSigmaUnc);
	    }
	}

    }
}

//pedestal summary
void doPedestalSummary(RooWorkspace *w,TCanvas *c,TString outDir)
{

  Int_t runNumber=w->var("runNumber")->getVal();
  for(size_t iest=0; iest<sizeof(siChargeEst)/sizeof(TString); iest++)
    {
      
      TString chargeEst(siChargeEst[iest]);
      for(Int_t ich=2; ich<4; ich++)
	{
	  Float_t pedestalAvg(w->var(Form("%s_pedestal_%d",chargeEst.Data(),ich))->getVal());
	  Float_t pedestalSigma(w->var(Form("%s_noise_%d",chargeEst.Data(),ich))->getVal());
	 
	  TGraph *pedestalGr=new TGraph();
	  pedestalGr->SetName(Form("pedestal_%s_%d",chargeEst.Data(),ich));
	  pedestalGr->SetFillColor(kGray);
	  pedestalGr->SetFillStyle(1001);
	  pedestalGr->SetLineColor(kGray);
	  pedestalGr->SetMarkerColor(kGray);
	  pedestalGr->SetPoint(0,-0.5,pedestalAvg-pedestalSigma);
	  pedestalGr->SetPoint(1,-0.5,pedestalAvg+pedestalSigma);
	  pedestalGr->SetPoint(2,2.5,pedestalAvg+pedestalSigma);
	  pedestalGr->SetPoint(3,2.5,pedestalAvg-pedestalSigma);
	  pedestalGr->SetPoint(4,-0.5,pedestalAvg-pedestalSigma);
	  TGraphErrors *pedestalEvolGr=new TGraphErrors();
	  pedestalEvolGr->SetName(Form("pedestalevol_%s_%d",chargeEst.Data(),ich));
	  pedestalEvolGr->SetMarkerStyle(20);
	  for(Float_t sigmaBias=0; sigmaBias<=2; sigmaBias+=1.0)
	    {
	      TString biasPostFix(Form("_bias%3.1f",sigmaBias));
	      Float_t pedestal(w->var(Form("%s_pedestal_%d%s",chargeEst.Data(),ich,biasPostFix.Data()))->getVal());
	      Float_t noise(w->var(Form("%s_noise_%d%s",chargeEst.Data(),ich,biasPostFix.Data()))->getVal());
	      Int_t np=pedestalEvolGr->GetN();
	      pedestalEvolGr->SetPoint(np,sigmaBias,pedestal);
	      pedestalEvolGr->SetPointError(np,0,noise);
	    }

	  c->Clear();
	  pedestalGr->Draw("af");
	  c->SetLogy(false);
	  pedestalGr->GetYaxis()->SetRangeUser(pedestalGr->GetYaxis()->GetXmin()*1.5,pedestalGr->GetYaxis()->GetXmax()*1.5);
	  pedestalGr->GetXaxis()->SetTitle("Bias applied to other Si pad (Nx#sigma_{noise})");
	  pedestalGr->GetYaxis()->SetTitle("Pedestal");
	  TLine *line = new TLine(-0.5,pedestalAvg,2.5,pedestalAvg);
	  line->Draw();
	  pedestalEvolGr->Draw("p");
	  TLatex txt;
	  txt.SetNDC();
	  txt.SetTextFont(42);
	  txt.SetTextSize(0.03);
	  txt.DrawLatex(0.15,0.93,CMSLABEL);
	  txt.DrawLatex(0.8,0.93,Form("Run %d", runNumber)); 
	  c->SaveAs(outDir+Form("/noiseSummary_%s_ch%d.png",chargeEst.Data(),ich));
	}
    }
}



//beamspot position
void determineFiducialBeamSpotFromWireChambers(TChain* H4treeReco,RooWorkspace *w,TCanvas *c,TString outDir,bool doSiFiducial)
{

  Int_t runNumber=w->var("runNumber")->getVal();

  Double_t bsq_x[4]   ={0.05, 0.50, 0.90, 0.95};
  Double_t bsq_val[4] ={0.00, 0.00, 0.00, 0.00};
  for(Int_t iwc=0; iwc<2; iwc++)
    {
      TString signalInWCCut(Form("wc_xl_hits[%d]>0 && wc_xr_hits[%d]>0 && wc_yu_hits[%d]>0 && wc_yd_hits[%d]>0",iwc,iwc,iwc,iwc));
      
      //require a signal in the MCP
      RooRealVar *baselineMCPVar=w->var(Form("%s_pedestal_1",mcpChargeEst.Data()));
      RooRealVar *baselineWidthMCPVar=w->var(Form("%s_noise_1",mcpChargeEst.Data()));
      Float_t baselineMCP(baselineMCPVar ?  baselineMCPVar->getVal() : 0.0);
      Float_t baselineWidthMCP(baselineWidthMCPVar ?  baselineWidthMCPVar->getVal() : 0.0);
      TString signalInMCPCut(Form("%s[1]-%f>50*%f",mcpChargeEst.Data(),baselineMCP,baselineWidthMCP));
    
      for(size_t iest=0; iest<sizeof(siChargeEst)/sizeof(TString); iest++)
        {
	  if(!doSiFiducial && iest>0) continue;
	  
	  TString finalCut(signalInWCCut + " && " + signalInMCPCut);

	  //add requirement for signal in the Si      
	  TString tag("inc");
	  if(doSiFiducial)
	    {
	      TString chargeEst(siChargeEst[iest]);
	      tag="fid_"+chargeEst;

	      finalCut +=" && ";
	      RooRealVar *baselineSiVar=w->var(Form("%s_pedestal_2",chargeEst.Data()));
	      RooRealVar *baselineWidthSiVar=w->var(Form("%s_noise_2",chargeEst.Data()));
	      Float_t baselineSi1(baselineSiVar ?  baselineSiVar->getVal() : 0.0);
	      Float_t baselineWidthSi1(baselineWidthSiVar ? baselineWidthSiVar->getVal():0.0);
	      finalCut += Form("%s[2]-%f>2*%f",chargeEst.Data(),baselineSi1,baselineWidthSi1);
	      
	      finalCut +=" && ";
	      baselineSiVar=w->var(Form("%s_pedestal_3",chargeEst.Data()));
	      baselineWidthSiVar=w->var(Form("%s_noise_3",chargeEst.Data()));
	      Float_t baselineSi2(baselineSiVar ?  baselineSiVar->getVal() : 0.0);
	      Float_t baselineWidthSi2(baselineWidthSiVar ? baselineWidthSiVar->getVal():0.0);
	      finalCut += Form("%s[3]-%f>2*%f",chargeEst.Data(),baselineSi2,baselineWidthSi2);
	    }
	  
	  TGraphAsymmErrors *bsxEvolGr=new TGraphAsymmErrors();
	  bsxEvolGr->SetName(Form("beamspotxevol_%s_%d",tag.Data(),iwc));
	  bsxEvolGr->SetMarkerStyle(20);
	  bsxEvolGr->SetFillStyle(0);
	  
	  TGraphAsymmErrors *bsyEvolGr=(TGraphAsymmErrors *)bsxEvolGr->Clone(Form("beamspotxevol_%s_%d",tag.Data(),iwc));
      
	  //project beamspot radius
	  H4treeReco->SetAlias(Form("wc_recorho%d",iwc),Form("sqrt(pow(wc_recoy[%d],2)+pow(wc_recox[%d],2))",iwc,iwc));
	  H4treeReco->Draw(Form("wc_recorho%d >> beamspotradius",iwc),finalCut);
	  TH1* beamspotradius=(TH1*)gROOT->FindObject("beamspotradius");
	  if(beamspotradius->Integral()<10) continue;

	  //determine 95% quantile to cut on the radius
	  beamspotradius->GetQuantiles(4,bsq_val,bsq_x);
	  Float_t maxRadius=bsq_val[3];
	  
	  signalInWCCut += Form(" && wc_recorho%d<%f",iwc,maxRadius);
	  H4treeReco->Draw(Form("wc_recoy[%d]:spillNumber>>beamspotyvsspill(500,0,500,200,-100,100)",iwc),finalCut);
	  TH2* beamspotyvsspill=(TH2*)gROOT->FindObject("beamspotyvsspill");
	  H4treeReco->Draw(Form("wc_recox[%d]:spillNumber>>beamspotxvsspill(500,0,500,200,-100,100)",iwc),finalCut);
	  TH2* beamspotxvsspill=(TH2*)gROOT->FindObject("beamspotxvsspill");

	  //profile as function of the spill number
	  Int_t minSpill(99999),maxSpill(-1);      
	  for(Int_t xbin=0; xbin<=beamspotyvsspill->GetNbinsX(); xbin++)
	    {
	      Int_t firstxbin(0),lastxbin(-1);
	      if(xbin>0) {firstxbin=xbin; lastxbin=xbin;}
	      TH1 *projy=beamspotyvsspill->ProjectionY("_py",firstxbin,lastxbin);
	      TH1 *projx=beamspotxvsspill->ProjectionY("_px",firstxbin,lastxbin);

	      Int_t spillNb=beamspotyvsspill->GetXaxis()->GetBinLowEdge(xbin);

	      //at least 10 events to do the control
	      if(projy->Integral()>10)
		{
		  if(spillNb<minSpill) minSpill=spillNb;
		  if(spillNb>maxSpill) maxSpill=spillNb;
		  
		  //y
		  projy->GetQuantiles(4,bsq_val,bsq_x);
		  if(xbin==0) 
		    {
		      w->factory(Form("y005_%s_%d[%f]",tag.Data(),iwc,bsq_val[0]));
		      w->factory(Form("y050_%s_%d[%f]",tag.Data(),iwc,bsq_val[1]));
		      w->factory(Form("y095_%s_%d[%f]",tag.Data(),iwc,bsq_val[3]));
		    }
		  else
		    {
		      Int_t np=bsyEvolGr->GetN();
		      bsyEvolGr->SetPoint(np,spillNb,bsq_val[1]);
		      bsyEvolGr->SetPointError(np,0,0,bsq_val[1]-bsq_val[0],bsq_val[3]-bsq_val[1]);
		    }
		}
	      projy->Delete();

	      //x
	      if(projx->Integral()>0)
		{
		  if(spillNb<minSpill) minSpill=spillNb;
		  if(spillNb>maxSpill) maxSpill=spillNb;
		  
		  projx->GetQuantiles(4,bsq_val,bsq_x);
		  if(xbin==0) 
		    {
		      w->factory(Form("x005_%s_%d[%f]",tag.Data(),iwc,bsq_val[0]));
		      w->factory(Form("x050_%s_%d[%f]",tag.Data(),iwc,bsq_val[1]));
		      w->factory(Form("x095_%s_%d[%f]",tag.Data(),iwc,bsq_val[3]));
		    }
		  else
		    {
		      Int_t np=bsxEvolGr->GetN();
		      bsxEvolGr->SetPoint(np,spillNb,bsq_val[1]);
		      bsxEvolGr->SetPointError(np,0,0,bsq_val[1]-bsq_val[0],bsq_val[3]-bsq_val[1]);
		    }
		}
	      projx->Delete();
	    }

	  //y summary
	  c->Clear();
	  TGraph *yrangeGr=new TGraph();
	  yrangeGr->SetName(Form("yrange_%s_%d",tag.Data(),iwc));
	  yrangeGr->SetFillColor(kGray);
	  yrangeGr->SetFillStyle(1001);
	  yrangeGr->SetLineColor(kGray);
	  yrangeGr->SetMarkerColor(kGray);
	  yrangeGr->SetPoint(0,minSpill,w->var(Form("y005_inc_%d",iwc))->getVal());
	  yrangeGr->SetPoint(1,minSpill,w->var(Form("y095_inc_%d",iwc))->getVal());
	  yrangeGr->SetPoint(2,maxSpill,w->var(Form("y095_inc_%d",iwc))->getVal());
	  yrangeGr->SetPoint(3,maxSpill,w->var(Form("y095_inc_%d",iwc))->getVal());
	  yrangeGr->SetPoint(4,maxSpill,w->var(Form("y005_inc_%d",iwc))->getVal());
	  yrangeGr->Draw("af");
	  //      yrangeGr->GetYaxis()->SetRangeUser(yrangeGr->GetYaxis()->GetXmin()*1.5,yrangeGr->GetYaxis()->GetXmax()*1.5);
	  yrangeGr->GetXaxis()->SetTitle("Spill number");
	  yrangeGr->GetYaxis()->SetTitle("y");
	  TLine *line = new TLine();
	  line->SetLineStyle(9);
	  line->DrawLine(minSpill,w->var(Form("y050_%s_%d",tag.Data(),iwc))->getVal(),maxSpill,w->var(Form("y050_%s_%d",tag.Data(),iwc))->getVal());	 
	  line->DrawLine(minSpill,w->var(Form("y005_%s_%d",tag.Data(),iwc))->getVal(),maxSpill,w->var(Form("y005_%s_%d",tag.Data(),iwc))->getVal());	 
	  line->DrawLine(minSpill,w->var(Form("y095_%s_%d",tag.Data(),iwc))->getVal(),maxSpill,w->var(Form("y095_%s_%d",tag.Data(),iwc))->getVal());	 
	  bsyEvolGr->Draw("p");
	  TLatex txt;
	  txt.SetNDC();
	  txt.SetTextFont(42);
	  txt.SetTextSize(0.03);
	  txt.DrawLatex(0.15,0.93,CMSLABEL);
	  txt.DrawLatex(0.8,0.93,Form("Run %d", runNumber)); 
	  if(doSiFiducial)
	    c->SaveAs(outDir+Form("/y_%s_wc%d.png",tag.Data(),iwc));
	  else
	    c->SaveAs(outDir+Form("/yinc_%s_wc%d.png",tag.Data(),iwc));
	  
	  //x summary
	  c->Clear();
	  TGraph *xrangeGr=(TGraph *) yrangeGr->Clone(Form("xrange_%s_%d",tag.Data(),iwc));
	  xrangeGr->SetPoint(0,minSpill,w->var(Form("x005_inc_%d",iwc))->getVal());
	  xrangeGr->SetPoint(1,minSpill,w->var(Form("x095_inc_%d",iwc))->getVal());
	  xrangeGr->SetPoint(2,maxSpill,w->var(Form("x095_inc_%d",iwc))->getVal());
	  xrangeGr->SetPoint(3,maxSpill,w->var(Form("x095_inc_%d",iwc))->getVal());
	  xrangeGr->SetPoint(4,maxSpill,w->var(Form("x005_inc_%d",iwc))->getVal());
	  xrangeGr->Draw("af");
	  //xrangeGr->GetYaxis()->SetRangeUser(xrangeGr->GetYaxis()->GetXmin()*1.5,xrangeGr->GetYaxis()->GetXmax()*1.5);
	  xrangeGr->GetXaxis()->SetTitle("Spill number");
	  xrangeGr->GetYaxis()->SetTitle("x");
	  line->DrawLine(minSpill,w->var(Form("x005_%s_%d",tag.Data(),iwc))->getVal(),maxSpill,w->var(Form("x005_%s_%d",tag.Data(),iwc))->getVal());
	  line->DrawLine(minSpill,w->var(Form("x050_%s_%d",tag.Data(),iwc))->getVal(),maxSpill,w->var(Form("x050_%s_%d",tag.Data(),iwc))->getVal());
	  line->DrawLine(minSpill,w->var(Form("x095_%s_%d",tag.Data(),iwc))->getVal(),maxSpill,w->var(Form("x095_%s_%d",tag.Data(),iwc))->getVal());
	  bsxEvolGr->Draw("p");
	  txt.DrawLatex(0.15,0.93,CMSLABEL);
	  txt.DrawLatex(0.8,0.93,Form("Run %d", runNumber)); 
	  if(doSiFiducial)
	    c->SaveAs(outDir+Form("/x_%s_wc%d.png",tag.Data(),iwc));
	  else
	    c->SaveAs(outDir+Form("/xinc_%s_wc%d.png",tag.Data(),iwc));
	}
    }
}

//
void runMIPFits(Int_t run,TString inDir,Bool_t allowSecondPeak)
{
  setROOTstyle();

  TString outDir=Form("%s/Run%d",inDir.Data(),run);
  TFile *fIn=TFile::Open(Form("%s/workspace.root",outDir.Data()));
  RooWorkspace *w=(RooWorkspace *) fIn->Get("w");
  w->Print("v");

  TCanvas* c=new TCanvas("c","c",500,500);
  c->SetTopMargin(0.02);
  c->SetRightMargin(0.02);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);

  //prepare summary tree
  TFile *fOut=TFile::Open(outDir+Form("/fitsummary.root"),"UPDATE");
  MIPFitSummary_t fitSummary;
  TTree *tree=new TTree("mipcalib","mipcalib");
  tree->Branch("runNumber",       &run,                        "runNumber/I");
  tree->Branch("chargeEstimator", &fitSummary.chargeEstimator, "chargeEstimator/I");
  tree->Branch("siPad",           &fitSummary.siPad,           "siPad/I");
  tree->Branch("fitType",         &fitSummary.fitType,         "fitType/I");
  tree->Branch("ndof",            &fitSummary.ndof,            "ndof/I");
  tree->Branch("chi2",            &fitSummary.chi2,            "chi2/F");
  tree->Branch("prob",            &fitSummary.prob,            "prob/F");
  tree->Branch("adc2mip",         &fitSummary.adc2mip,         "adc2mip/F");
  tree->Branch("adc2mipUnc",      &fitSummary.adc2mipUnc,      "adc2mipUnc/F");
  tree->Branch("mip1frac",        &fitSummary.mip1frac,        "mip1frac/F");
  tree->Branch("mip1fracUnc",     &fitSummary.mip1fracUnc,     "mip1fracUnc/F");
  tree->Branch("mip2frac",        &fitSummary.mip2frac,        "mip2frac/F");
  tree->Branch("mip2fracUnc",     &fitSummary.mip2fracUnc,     "mip2fracUnc/F");
  tree->Branch("relWidth",        &fitSummary.relWidth,        "relWidth/F");
  tree->Branch("relWidthUnc",     &fitSummary.relWidthUnc,     "relWidthUnc/F");

  //loop over charge estimators
  for(size_t iest=0; iest<sizeof(siChargeEst)/sizeof(TString); iest++)
    {
      RooDataSet *data=(RooDataSet *)w->data("data_"+siChargeEst[iest]);
      fitSummary.chargeEstimator=iest;

      //loop over pads
      for(Int_t ich=2; ich<4; ich++)
	{
	  fitSummary.siPad=ich;

	  for(Int_t ifit=0; ifit<5; ifit++)
	    {
	      fitSummary.fitType=ifit;

	      //get variable to be fit and associated noise
	      TString noiseName(Form("%s_noise_%d",siChargeEst[iest].Data(),ich));
	      TString varName(Form("charge_si%d",ich));
	      if(ifit==2)  { varName+="_bias0.0"; noiseName+= "_bias0.0"; }
	      if(ifit==3)  { varName+="_bias1.0"; noiseName+= "_bias1.0"; }
	      if(ifit==4)  { varName+="_bias2.0"; noiseName+= "_bias1.0"; }
	      RooRealVar *var=(RooRealVar *) w->var(varName);
  	      RooRealVar *noiseVar=(RooRealVar *)w->var(noiseName);

	      //selection cut
	      TString redCut(Form("%s>-1000",varName.Data()));
	      if(ifit>0)   redCut += " && isFiducial==1";

	      //reduce the original dataset
	      RooDataSet *redData=(RooDataSet *)data->reduce(RooArgSet(*var),redCut);
	      Float_t minCharge(iest==0?-100:-50), maxCharge(iest==0?1000:250);
	      TH1 *h=redData->createHistogram("binneddata",
					      *var,
					      RooFit::Binning(100,minCharge,maxCharge));
	      h->SetLineColor(1);
	      h->SetMarkerStyle(20);
	      h->SetMarkerColor(1);
	      TF1* fitFunc=new TF1("fitFunc",sigFunc,minCharge,maxCharge,11);
	      
	      fitFunc->SetParName(0,"N_{0}");
	      fitFunc->SetParName(1,"pedestal");
	      fitFunc->SetParName(2,"#sigma_{noise}");
	      fitFunc->SetParName(3,"#xi_{1}");
	      fitFunc->SetParName(4,"MPV_{1}");
	      fitFunc->SetParName(5,"N_{1}");
	      fitFunc->SetParName(6,"#sigma_{noise}");
	      if(allowSecondPeak)
		{
		  fitFunc->SetParName(7,"#xi_{2}");
		  fitFunc->SetParName(8,"MPV_{2}");
		  fitFunc->SetParName(9,"N_{2}");
		  fitFunc->SetParName(10,"#sigma_{noise}");
		}

	      fitFunc->SetParLimits(0,0.,h->Integral()); //h->GetXaxis()->FindBin(0.),2*h->Integral(1,h->GetXaxis()->FindBin(0.)));
	      fitFunc->FixParameter(1,0.);
	      if(ifit==4) fitFunc->SetParLimits(1,-noiseVar->getVal(),noiseVar->getVal());
	      fitFunc->FixParameter(2,noiseVar->getVal());

	      fitFunc->SetParLimits(3,(iest==0 ? 9 : 1),  (iest==0 ? 35 : 20));
	      fitFunc->SetParLimits(4,(iest==0 ? 90 : 5), (iest==0 ? 350 : 100));
	      fitFunc->SetParLimits(5,0,h->Integral()*1000);
	      fitFunc->FixParameter(6,noiseVar->getVal());
	      
	      fitFunc->FixParameter(7,0.1);
	      fitFunc->FixParameter(8,0.);
	      if(allowSecondPeak) fitFunc->SetParLimits(9,0,h->Integral()*10);
	      else                fitFunc->FixParameter(9,0.);
	      fitFunc->FixParameter(10,noiseVar->getVal());

	      TFitResultPtr fitRes=h->Fit(fitFunc,"SERB+");
	      //Int_t status = (Int_t)fitRes;

	      //show the result of the fit
	      c->Clear();

	      fitFunc->Draw();
	      fitFunc->SetLineColor(kBlue);

	      Double_t fitparams[11];
	      fitFunc->GetParameters(fitparams);
	      TF1 *mip1= new TF1("mip1",lanconvgau,minCharge,maxCharge,4);
	      mip1->SetParameters(&fitparams[3]);
	      mip1->SetLineColor(kGray);
	      mip1->SetFillColor(kGray);
	      mip1->SetFillStyle(1001);
	      mip1->Draw("fcsame");
	      
	      if(allowSecondPeak)
		{
		  TF1 *mip2= new TF1("mip2",lanconvgau,minCharge,maxCharge,4);
		  Double_t par[4]={1.66*TMath::Sqrt(3.0)*fitparams[3],1.22*3*fitparams[4],fitparams[9],fitparams[2]};
		  mip2->SetParameters(&par[0]);
		  mip2->SetLineColor(kGray+1);
		  mip2->SetFillColor(kGray+1);
		  mip2->SetFillStyle(3001);
		  mip2->Draw("fcsame");
		}
	      
	      h->Draw("same");
	      fitFunc->GetYaxis()->SetTitle("Events");
	      fitFunc->GetXaxis()->SetTitle(iest==0 ? 
					    "Charge sum [ADC]" :
					    "Wave amplitude [ADC]");
	      
	      if(ifit<2) 
		{
		  c->SetLogy(true);
		  fitFunc->GetYaxis()->SetRangeUser(0.1,h->GetMaximum()*10);
		}
	      else
		{
		  c->SetLogy(false);
		  fitFunc->GetYaxis()->SetRangeUser(0.1,h->GetMaximum()*1.3);
		}

	      TLatex txt;
	      txt.SetNDC();
	      txt.SetTextFont(42);
	      txt.SetTextSize(0.03);
	      txt.DrawLatex(0.15,0.93,CMSLABEL);
	      txt.DrawLatex(0.8,0.93,Form("Run %d",run));
	      TString fitName("inclusive");
	      if(ifit>0)  fitName = "fiducial,";
	      if(ifit>=2) fitName += Form(">%d#sigma_{noise}",ifit-2);
	      txt.DrawLatex(0.15,0.88,Form("#scale[0.8]{Si%d %s}",ich,fitName.Data()));
	      	      
	      Float_t chi2=fitFunc->GetChisquare();
	      Int_t ndof=fitFunc->GetNDF();
	      txt.DrawLatex(0.6,0.88,Form("#chi^{2}/ndof = %3.3f/%d",chi2,ndof));
	      Float_t prob=fitFunc->GetProb();
	      txt.DrawLatex(0.6,0.83,Form("p-value = %3.3f",prob));
	      Float_t adc2mip(fitFunc->GetParameter(4)),adc2mipUnc(fitFunc->GetParError(4));
	      txt.DrawLatex(0.6,0.78,Form("MIP/ADC = %3.3f#pm%3.3f",adc2mip,adc2mipUnc));

	      TMatrixDSym cov=fitRes->GetCovarianceMatrix();
	      std::pair<Float_t,Float_t> relWidth = getParameterRatio(fitFunc->GetParameter(3),fitFunc->GetParError(3),
								      fitFunc->GetParameter(4),fitFunc->GetParError(4),
								      cov(3,4));
	      txt.DrawLatex(0.6,0.73,Form("#xi/MIP = %3.3f#pm%3.3f",relWidth.first,relWidth.second));


	      Float_t normFactor( noiseVar->getVal()*TMath::Sqrt(2*TMath::Pi()) );
	      std::pair<Float_t,Float_t> mip1frac = getParameterRatio(fitFunc->GetParameter(5),fitFunc->GetParError(5),
								      normFactor*fitFunc->GetParameter(0),normFactor*fitFunc->GetParError(0),
								      cov(0,5));
	      txt.DrawLatex(0.6,0.68,Form("S(1)/N = %3.3f#pm%3.3f",mip1frac.first,mip1frac.second));
	      std::pair<Float_t,Float_t> mip2frac(-1,0);
	      if(allowSecondPeak)
		{
		  mip2frac = getParameterRatio(fitFunc->GetParameter(9),fitFunc->GetParError(9),
					       fitFunc->GetParameter(5),fitFunc->GetParError(5),
					       cov(0,9));
		  txt.DrawLatex(0.6,0.63,Form("S(3)/S(1) = %3.3f#pm%3.3f",mip2frac.first,mip2frac.second));
		}
	      
	      //if(status!=0) txt.DrawLatex(0.6,0.53,Form("#scale[0.8]{#color[2]{Warning status=%d}}",status));
	      c->SaveAs(outDir+Form("/mipfit_pad%d_fit%d_%s.png",ich,ifit,siChargeEst[iest].Data()));

	      //save summary
	      fitSummary.chi2         = chi2;
	      fitSummary.ndof         = ndof;
	      fitSummary.prob         = prob;
	      fitSummary.adc2mip      = adc2mip;
	      fitSummary.adc2mipUnc   = adc2mipUnc;
	      fitSummary.relWidth     = relWidth.first;
	      fitSummary.relWidthUnc  = relWidth.second;
	      fitSummary.mip1frac     = mip1frac.first;
	      fitSummary.mip1fracUnc  = mip1frac.second;
	      fitSummary.mip2frac     = mip2frac.first;
	      fitSummary.mip2fracUnc  = mip2frac.second;
	      tree->Fill();
	    }
	}
    }

  //close output file
  tree->Write();
  fOut->Close();
}

//
void prepareWorkspaceForRun(TString url="root://eoscms//eos/cms/store/cmst3/group/hgcal/TimingTB_H2_Jul2015/RECO/3eb93c8/RECO_3373.root",TString outDir="~/www/HGCal")
{

  setROOTstyle();

  //open file
  TChain *H4treeReco=new TChain("H4treeReco");
  H4treeReco->Add(url);
  UInt_t runNumber;
  H4treeReco->SetBranchAddress("runNumber",    &runNumber);
  H4treeReco->GetEntry(0);

  //prepare workspace
  RooWorkspace *w=new RooWorkspace("w");
  w->factory(Form("runNumber[%d]",runNumber));

  outDir+=Form("/Run%d/",runNumber);
  gSystem->Exec("mkdir -p "+outDir);

  TCanvas* c=new TCanvas("c","c",500,500);
  c->SetTopMargin(0.02);
  c->SetRightMargin(0.02);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);

  //pedestals (2nd time computes residuals as they'll be stored in the workspace already)
  for(Float_t sigmaBias=-1; sigmaBias<=2; sigmaBias+=1.0)
    {
      doPedestals(H4treeReco,w,sigmaBias,c,outDir);
      doPedestals(H4treeReco,w,sigmaBias,c,outDir);
    }
  doPedestalSummary(w,c,outDir);

  //beamspot
  determineFiducialBeamSpotFromWireChambers(H4treeReco,w,c,outDir,false);
  determineFiducialBeamSpotFromWireChambers(H4treeReco,w,c,outDir,true);

  //prepare to read the tree
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

  //prepare variables to store in the datasets
  RooArgSet dataVars;
  dataVars.add( *((RooRealVar *)w->factory("isFiducial[0,0,1]")) );
  for(Float_t sigmaBias=-1; sigmaBias<=2; sigmaBias+=1.0)
    {
      TString biasPostFix(sigmaBias>=0 ? Form("_bias%3.1f",sigmaBias) : "");
      for(int ich=2; ich<4; ich++)
	{
	  TString name(Form("charge_si%d%s",ich,biasPostFix.Data()));
	  dataVars.add( *((RooRealVar *)w->factory(Form("%s[0,-999999999,999999999]",name.Data()))) );
	}
    }

  //loop over charge estimators
  for(size_t iest=0; iest<sizeof(siChargeEst)/sizeof(TString); iest++)
    {
      RooDataSet data(Form("data_%s",siChargeEst[iest].Data()),Form("data_%s",siChargeEst[iest].Data()),dataVars);

      //loop over events
      for(int i=0; i<H4treeReco->GetEntries(); i++)
	{
	  H4treeReco->GetEntry(i);

	  //require hits in the wire chambers
	  Bool_t isEmpty(wc_xl_hits[0]==0 && wc_xr_hits[0]==0 && wc_xl_hits[1]==0 && wc_xr_hits[1]==0 && wc_yu_hits[1]==0 && wc_yd_hits[1]==0);
	  if(isEmpty) continue;

	  //require reconstructed beamspot to be meaningful
	  Bool_t isInclusiveBeamspot(wc_recox[1]>w->var("x005_inc_1")->getVal() && 
				     wc_recox[1]<w->var("x095_inc_1")->getVal() && 
				     wc_recoy[1]>w->var("y005_inc_1")->getVal() && 
				     wc_recoy[1]<w->var("y095_inc_1")->getVal());
	  if(!isInclusiveBeamspot) continue;

	  //has MCP trigger
	  Bool_t hasMCPTrigger(wave_max[1]-w->var("wave_max_pedestal_1")->getVal()>50*w->var("wave_max_noise_1")->getVal());
	  if(!hasMCPTrigger) continue;

	  //fiducial in Si
	  Bool_t isFiducialBeamspot(wc_recox[1]>w->var(Form("x005_fid_%s_1",siChargeEst[iest].Data()))->getVal() && 
				    wc_recox[1]<w->var(Form("x095_fid_%s_1",siChargeEst[iest].Data()))->getVal() && 
				    wc_recoy[1]>w->var(Form("y005_fid_%s_1",siChargeEst[iest].Data()))->getVal() && 
				    wc_recoy[1]<w->var(Form("y095_fid_%s_1",siChargeEst[iest].Data()))->getVal());
	  ((RooRealVar *)dataVars.find("isFiducial"))->setVal( isFiducialBeamspot );
	  
	  //loop over channels
	  for(int ich=2; ich<4; ich++)
	    {
	      //cout << "Channel=" << ich << " ";

	      //check charge on the other Si pad
	      Int_t otherSiPadCh(ich==2? 3 : 2);
	      Float_t otherRawCharge(charge_integ_smallw_mcp[otherSiPadCh]);
	      if(iest==1) otherRawCharge=wave_fit_smallw_ampl[otherSiPadCh];
	      RooRealVar *baselineOtherVar=w->var(Form("%s_pedestal_%d",siChargeEst[iest].Data(),otherSiPadCh));
	      RooRealVar *baselineOtherWidthVar=w->var(Form("%s_noise_%d",siChargeEst[iest].Data(),otherSiPadCh));
	      Float_t otherCharge(otherRawCharge-baselineOtherVar->getVal());

	      //get charge
	      Float_t rawCharge(charge_integ_smallw_mcp[ich]);
	      if(iest==1) rawCharge=wave_fit_smallw_ampl[ich];
	      for(Float_t sigmaBias=-1; sigmaBias<=2; sigmaBias+=1.0)
		{
		  TString biasPostFix("");
		  Bool_t passBiasCut(true);
		  if(sigmaBias>=0)  
		    {
		      biasPostFix=Form("_bias%3.1f",sigmaBias);
		      passBiasCut=(otherCharge>baselineOtherWidthVar->getVal());
		    }

		  RooRealVar *baselineVar=w->var(Form("%s_pedestal_%d%s",siChargeEst[iest].Data(),ich,biasPostFix.Data()));
		  Float_t charge(passBiasCut ? rawCharge-baselineVar->getVal() : -999999.);
		  TString name(Form("charge_si%d%s",ich,biasPostFix.Data()));
		  ((RooRealVar *)dataVars.find(name.Data()))->setVal(charge);	
		}
	    }

	  data.add(dataVars);
	}

      //import dataset
      w->import(data);
      data.Print("v");
    }

  w->Print("v");
  w->writeToFile(outDir+"/workspace.root");
}


//wrapper for all the MIP runs
void makeRunSummary(TString outDir="~/www/HGCal/FastTimeTB/MIPCalibrationRuns",
		    TString inDir="/store/cmst3/group/hgcal/TimingTB_H2_Jul2015/RECO/3eb93c8",
		    Bool_t redoWorkspace=true)
{
  Int_t runs[]    = {3370,3371,3372,3373,3374,3375,3376,3363,3369,3346};
  for(size_t i=0; i<sizeof(runs)/sizeof(Int_t); i++)
    {
      TString url(Form("root://eoscms//eos/cms/%s/RECO_%d.root",inDir.Data(),runs[i]));
      if(redoWorkspace) prepareWorkspaceForRun(url,outDir);
      
      bool doSecondPeak(runs[i]==3363 || runs[i]==3369 || runs[i]==3346); //these are electron runs 
      runMIPFits(runs[i],outDir,doSecondPeak);

      gSystem->Exec(Form("cp /afs/cern.ch/user/p/psilva/work/HGCal/TestBeam/Analysis/FastTimeTestBeamAnalysis/scripts/index.php %s/Run%d",outDir.Data(),runs[i]));
    }
}


