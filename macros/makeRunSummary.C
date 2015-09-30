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

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include "fitFunctions.C"

using namespace std;

struct RunSummary_t
{
  Int_t runNumber;
  std::map<Int_t, Float_t> pedestals,pedestalsSigma;
  std::map<Int_t, Float_t> pedestalsUnc,pedestalsSigmaUnc;
  std::map<Int_t, Float_t> x005,x050,x095;
  std::map<Int_t, Float_t> y005,y050,y095;
  std::map<Int_t, Float_t> xinc005,xinc050,xinc095;
  std::map<Int_t, Float_t> yinc005,yinc050,yinc095;
  std::map<Int_t, Float_t> mip1,mip1Unc,mip2,mip2Unc,chi2ndof1,chi2ndof2;
};

TString CMSLABEL="#splitline{#bf{CMS} #it{work in progress}}{#scale[0.8]{#it{Fast-timing testbeam}}}";
RunSummary_t runSummary;
enum MIPFitMode_t { FIXEDNOISE=0, CORRELATEDNOISE=1, INDEPENDENTNOISE=2 };


//pedestal control in empty events
void doPedestals(TTree* H4treeReco,TCanvas *c,TString outDir,TString siChargeEst="charge_integ_smallw_mcp")
{

  std::map<Int_t,Float_t> pedestalsPerChannel;
  TString emptyEvtCut("wc_xl_hits[0]==0 && wc_xr_hits[0]==0 && wc_xl_hits[1]==0 && wc_xr_hits[1]==0 && wc_yu_hits[1]==0 && wc_yd_hits[1]==0");
  for(Int_t ich=1; ich<4; ich++)
    {
      TGraphErrors *pedestalEvolGr=new TGraphErrors();
      pedestalEvolGr->SetName(Form("pedestalevol_%d",ich));
      pedestalEvolGr->SetMarkerStyle(20);
      pedestalEvolGr->SetFillStyle(0);

      TString chargeEst(siChargeEst);
      if(ich==1) chargeEst="wave_max";

      //project pedestal vs spill to trace the evolution throughout the run, estimate from gaussian fit
      Float_t baseline(0);
      if(runSummary.pedestals.find(ich) != runSummary.pedestals.end()) baseline=runSummary.pedestals[ich];
      H4treeReco->Draw(Form("%s[%d]-%f:spillNumber>>pedestalvsspill(500,0,500,250,-250,250)",chargeEst.Data(),ich,baseline),emptyEvtCut);
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
	      txt.DrawLatex(0.8,0.93,Form("Run %d", runSummary.runNumber)); 
	      if(runSummary.pedestals.find(ich)!=runSummary.pedestals.end())
		c->SaveAs(outDir+Form("/residualincnoisefit_ch%d.png",ich));
	      else
		c->SaveAs(outDir+Form("/innoisefit_ch%d.png",ich));
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
      pedestalGr->SetName(Form("pedestal_%d",ich));
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
      txt.DrawLatex(0.8,0.93,Form("Run %d", runSummary.runNumber)); 
      if(runSummary.pedestals.find(ich)!=runSummary.pedestals.end())
	c->SaveAs(outDir+Form("/residualnoise_ch%d.png",ich));
      else
	{
	  c->SaveAs(outDir+Form("/noise_ch%d.png",ich));
	  runSummary.pedestals[ich]=pedestalAvg;
	  runSummary.pedestalsSigma[ich]=pedestalSigma;
	  runSummary.pedestalsUnc[ich]=pedestalAvgUnc;
	  runSummary.pedestalsSigmaUnc[ich]=pedestalSigmaUnc;
	}
   }

}

//beamspot position
void determineFiducialBeamSpotFromWireChambers(TTree* H4treeReco,TCanvas *c,TString outDir,bool doSiFiducial,TString siChargeEst="charge_integ_smallw_mcp")
{

  Double_t bsq_x[4]   ={0.05, 0.50, 0.90, 0.95};
  Double_t bsq_val[4] ={0.00, 0.00, 0.00, 0.00};
  for(Int_t iwc=0; iwc<2; iwc++)
    {
      TString signalInWCCut(Form("wc_xl_hits[%d]>0 && wc_xr_hits[%d]>0 && wc_yu_hits[%d]>0 && wc_yd_hits[%d]>0",iwc,iwc,iwc,iwc));
      
      //require a signal in the MCP
      Float_t baselineMCP(runSummary.pedestals.find(1)!=runSummary.pedestals.end() ? runSummary.pedestals[1] : 0.0);
      Float_t baselineWidthMCP(runSummary.pedestalsSigma.find(1)!=runSummary.pedestalsSigma.end() ? runSummary.pedestalsSigma[1] : 0.0);
      TString signalInMCPCut(Form("wave_max[1]-%f>50*%f",baselineMCP,baselineWidthMCP));

      TString finalCut(signalInWCCut + " && " + signalInMCPCut);

      //add requirement for signal in the Si      
      if(doSiFiducial)
	{
	  finalCut +=" && ";

	  Float_t baselineSi1(runSummary.pedestals.find(2)!=runSummary.pedestals.end() ? runSummary.pedestals[2] : 0.0);
	  Float_t baselineWidthSi1(runSummary.pedestalsSigma.find(2)!=runSummary.pedestalsSigma.end() ? runSummary.pedestalsSigma[2] : 0.0);
	  finalCut += Form("%s[2]-%f>2*%f",siChargeEst.Data(),baselineSi1,baselineWidthSi1);

	  finalCut +=" && ";

	  Float_t baselineSi2(runSummary.pedestals.find(3)!=runSummary.pedestals.end() ? runSummary.pedestals[3] : 0.0);
	  Float_t baselineWidthSi2(runSummary.pedestalsSigma.find(3)!=runSummary.pedestalsSigma.end() ? runSummary.pedestalsSigma[3] : 0.0);
	  finalCut += Form("%s[3]-%f>2*%f",siChargeEst.Data(),baselineSi2,baselineWidthSi2);
	}

      TGraphAsymmErrors *bsxEvolGr=new TGraphAsymmErrors();
      bsxEvolGr->SetName(Form("beamspotxevol_%d",iwc));
      bsxEvolGr->SetMarkerStyle(20);
      bsxEvolGr->SetFillStyle(0);
      
      TGraphAsymmErrors *bsyEvolGr=(TGraphAsymmErrors *)bsxEvolGr->Clone(Form("beamspotxevol_%d",iwc));
      
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
		  if(doSiFiducial)
		    {
		      runSummary.y005[iwc]=bsq_val[0]; runSummary.y050[iwc]=bsq_val[1]; runSummary.y095[iwc]=bsq_val[3]; 
		    }
		  else
		    {
		      runSummary.yinc005[iwc]=bsq_val[0]; runSummary.yinc050[iwc]=bsq_val[1]; runSummary.yinc095[iwc]=bsq_val[3]; 
		    }
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
		  if(doSiFiducial)
		    {
		      runSummary.x005[iwc]=bsq_val[0]; runSummary.x050[iwc]=bsq_val[1]; runSummary.x095[iwc]=bsq_val[3]; 
		    }
		  else
		    {
		      runSummary.xinc005[iwc]=bsq_val[0]; runSummary.xinc050[iwc]=bsq_val[1]; runSummary.xinc095[iwc]=bsq_val[3]; 
		    }
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
      yrangeGr->SetName(Form("yrange_%d",iwc));
      yrangeGr->SetFillColor(kGray);
      yrangeGr->SetFillStyle(1001);
      yrangeGr->SetLineColor(kGray);
      yrangeGr->SetMarkerColor(kGray);
      yrangeGr->SetPoint(0,minSpill,runSummary.yinc005[iwc]);
      yrangeGr->SetPoint(1,minSpill,runSummary.yinc095[iwc]);
      yrangeGr->SetPoint(2,maxSpill,runSummary.yinc095[iwc]);
      yrangeGr->SetPoint(3,maxSpill,runSummary.yinc095[iwc]);
      yrangeGr->SetPoint(4,maxSpill,runSummary.yinc005[iwc]);
      yrangeGr->Draw("af");
      //      yrangeGr->GetYaxis()->SetRangeUser(yrangeGr->GetYaxis()->GetXmin()*1.5,yrangeGr->GetYaxis()->GetXmax()*1.5);
      yrangeGr->GetXaxis()->SetTitle("Spill number");
      yrangeGr->GetYaxis()->SetTitle("y");
      TLine *line = new TLine();
      if(doSiFiducial)
	line->DrawLine(minSpill,runSummary.y050[iwc],maxSpill,runSummary.y050[iwc]);
      else
	line->DrawLine(minSpill,runSummary.yinc050[iwc],maxSpill,runSummary.yinc050[iwc]);
      bsyEvolGr->Draw("p");
      TLatex txt;
      txt.SetNDC();
      txt.SetTextFont(42);
      txt.SetTextSize(0.03);
      txt.DrawLatex(0.15,0.93,CMSLABEL);
      txt.DrawLatex(0.8,0.93,Form("Run %d", runSummary.runNumber)); 
      if(doSiFiducial)
	c->SaveAs(outDir+Form("/y_wc%d.png",iwc));
      else
	c->SaveAs(outDir+Form("/yinc_wc%d.png",iwc));

      //x summary
      c->Clear();
      TGraph *xrangeGr=(TGraph *) yrangeGr->Clone(Form("xrange_%d",iwc));
      xrangeGr->SetPoint(0,minSpill,runSummary.xinc005[iwc]);
      xrangeGr->SetPoint(1,minSpill,runSummary.xinc095[iwc]);
      xrangeGr->SetPoint(2,maxSpill,runSummary.xinc095[iwc]);
      xrangeGr->SetPoint(3,maxSpill,runSummary.xinc095[iwc]);
      xrangeGr->SetPoint(4,maxSpill,runSummary.xinc005[iwc]);
      xrangeGr->Draw("af");
      //xrangeGr->GetYaxis()->SetRangeUser(xrangeGr->GetYaxis()->GetXmin()*1.5,xrangeGr->GetYaxis()->GetXmax()*1.5);
      xrangeGr->GetXaxis()->SetTitle("Spill number");
      xrangeGr->GetYaxis()->SetTitle("x");
      if(doSiFiducial)
	line->DrawLine(minSpill,runSummary.x050[iwc],maxSpill,runSummary.x050[iwc]);
      else
	line->DrawLine(minSpill,runSummary.xinc050[iwc],maxSpill,runSummary.xinc050[iwc]);
      line->Draw();
      bsxEvolGr->Draw("p");
      txt.DrawLatex(0.15,0.93,CMSLABEL);
      txt.DrawLatex(0.8,0.93,Form("Run %d", runSummary.runNumber)); 
      if(doSiFiducial)
	c->SaveAs(outDir+Form("/x_wc%d.png",iwc));
      else
	c->SaveAs(outDir+Form("/xinc_wc%d.png",iwc));
    }
}

//
void doMIPfits(TTree* H4treeReco,TCanvas *c,TString outDir,int mode=0,Int_t fitMode=CORRELATEDNOISE,TString siChargeEst="charge_integ_smallw_mcp")
{
  TString wcCut("wc_xl_hits[0]!=0 && wc_xr_hits[0]!=0 && wc_xl_hits[1]!=0 && wc_xr_hits[1]!=0");

  TF1* fitFunc=new TF1("fitFunc",fitMode==INDEPENDENTNOISE ? sigFunc : sigFuncFixedNoiseResol,-100,400,7);
  //fitFunc->SetNpx(25);
  fitFunc->SetParName(0,"N_{noise}");
  fitFunc->SetParName(1,"#mu_{noise}");
  fitFunc->SetParName(2,fitMode==INDEPENDENTNOISE ? "#sigma_{noise}^{2}" : "#sigma_{noise}");
  fitFunc->SetParName(3,"#sigma_{signal}");
  fitFunc->SetParName(4,"MPV_{signal}");
  fitFunc->SetParName(5,"N_{signal}");
  fitFunc->SetParName(6,fitMode==INDEPENDENTNOISE ? "#sigma_{noise}^{2}" : "#sigma_{noise}");


  for(Int_t ich=2; ich<4; ich++)
    {      
      Float_t baseline(runSummary.pedestals.find(ich)!=runSummary.pedestals.end() ? runSummary.pedestals[ich] : 0.0);
      Float_t baselineWidth(runSummary.pedestalsSigma.find(ich)!=runSummary.pedestalsSigma.end() ? runSummary.pedestalsSigma[ich] : 0.0);

      Int_t otherSiPadCh(ich==2? 3 : 2);
      Float_t baselineOtherSiPad(runSummary.pedestals.find(otherSiPadCh)!=runSummary.pedestals.end() ? runSummary.pedestals[otherSiPadCh] : 0.0);
      Float_t baselineWidthOtherSiPad(runSummary.pedestalsSigma.find(otherSiPadCh)!=runSummary.pedestalsSigma.end() ? runSummary.pedestalsSigma[otherSiPadCh] : 0.0);

      Float_t baselineMCP(runSummary.pedestals.find(1)!=runSummary.pedestals.end() ? runSummary.pedestals[1] : 0.0);
      Float_t baselineWidthMCP(runSummary.pedestalsSigma.find(1)!=runSummary.pedestalsSigma.end() ? runSummary.pedestalsSigma[1] : 0.0);

      for(Int_t nSigma=2; nSigma<=4; nSigma+=1)
	{
	  TString trigCut(Form("wave_max[1]-%f>50*%f && ",baselineMCP,baselineWidthMCP));
	  if(mode==0)
	    {
	      Float_t xmin(runSummary.xinc005[1]),xmax(runSummary.xinc095[1]),ymin(runSummary.yinc005[1]),ymax(runSummary.yinc095[1]);
	      if(nSigma==3)
		{
		  xmin=runSummary.x005[1]; xmax=runSummary.x095[1];
		  ymin=runSummary.y005[1]; ymax=runSummary.y095[1];
		}
	      if(nSigma==4)
		{
		  xmin=runSummary.x050[1]-1.25; xmax=runSummary.x050[1]+1.25;
		  ymin=runSummary.y050[1]-1.25; ymax=runSummary.y050[1]+1.25;
		}
	      trigCut += Form("wc_recox[1]>%f && wc_recox[1]<%f && wc_recoy[1]>%f && wc_recoy[1]<%f",xmin,xmax,ymin,ymax);
	    }
	  if(mode==1)
	    {
	      Float_t xmin(runSummary.x005[1]),xmax(runSummary.x095[1]),ymin(runSummary.y005[1]),ymax(runSummary.y095[1]);
	      trigCut += Form("wc_recox[1]>%f && wc_recox[1]<%f && wc_recoy[1]>%f && wc_recoy[1]<%f &&",xmin,xmax,ymin,ymax);
	      trigCut += Form("%s[%d]-%f>%d*%f",siChargeEst.Data(),otherSiPadCh,baselineOtherSiPad,nSigma,baselineWidthOtherSiPad);
	    }
	  
	  H4treeReco->Draw(Form("%s[%d]-%f >> signal(50,-100,400)",siChargeEst.Data(),ich,baseline),wcCut+" && "+trigCut);

	  TH1F* signal=(TH1F*)gROOT->FindObject("signal");
	  signal->Sumw2();
	  
	  fitFunc->SetParLimits(0,0,signal->Integral()*10000);
	  fitFunc->FixParameter(1,0);	  
	  fitFunc->FixParameter(2,baselineWidth);
	  fitFunc->SetParLimits(3,5,400);
	  fitFunc->SetParLimits(4,0,400);
	  fitFunc->SetParLimits(5,0,signal->Integral()*10000);
	  fitFunc->FixParameter(6,baselineWidth);
	  
	  if(fitMode!=FIXEDNOISE)
	    {
	      fitFunc->SetParLimits(1,-baselineWidth,baselineWidth);
	      fitFunc->SetParLimits(2,baselineWidth*0.8,baselineWidth*1.2);
	      if(fitMode==CORRELATEDNOISE)
		fitFunc->FixParameter(6,baselineWidth);
	      else
		fitFunc->SetParLimits(6,baselineWidth*0.8,baselineWidth*1.2);	     
	    }

	  signal->Fit("fitFunc","BRMQ");


	  if(mode==1)
	    {
	      if(ich==2)
		{
		  runSummary.mip1[nSigma]=fitFunc->GetParameter(4);
		  runSummary.mip1Unc[nSigma]=fitFunc->GetParError(4);
		  runSummary.chi2ndof1[nSigma]=fitFunc->GetChisquare()/fitFunc->GetNDF();
		}
	      else
		{
		  runSummary.mip2[nSigma]=fitFunc->GetParameter(4);
		  runSummary.mip2Unc[nSigma]=fitFunc->GetParError(4);
		  runSummary.chi2ndof2[nSigma]=fitFunc->GetChisquare()/fitFunc->GetNDF();
		}
	    }


	  signal->SetMarkerStyle(20);
	  signal->SetMarkerSize(0.8);

	  c->Clear();
	  signal->Draw("PE");
	  signal->GetYaxis()->SetTitle("Events");
	  signal->GetXaxis()->SetTitle("Integrated charge [ADC]");
	  fitFunc->SetLineWidth(3);
	  fitFunc->Draw("SAME");
	  Double_t fitparams[7];
	  fitFunc->GetParameters(fitparams);

	  TF1 *noi=new TF1("noi","gaus",-50,200);
	  noi->SetParameter(0,fitparams[0]);
	  noi->SetParameter(1,fitparams[1]);
	  noi->SetParameter(2,fitparams[2]);
	  noi->SetLineColor(kMagenta);
	  noi->SetLineStyle(2);
	  //noi->SetNpx(1000);
	  noi->Draw("SAME");

	  TLatex txt;
	  txt.SetNDC();
	  txt.SetTextFont(42);
	  txt.SetTextSize(0.03);
	  txt.DrawLatex(0.15,0.93,CMSLABEL);
	  txt.DrawLatex(0.8,0.93,Form("Run %d", runSummary.runNumber)); 
	  if(mode==0)
	    {
	      if(nSigma==2) txt.DrawLatex(0.15,0.88,Form("Si%d : |#Delta x|<5-95%% inc.", ich) );
	      if(nSigma==3) txt.DrawLatex(0.15,0.88,Form("Si%d : |#Delta x|<5-95%% fid.", ich) );
	      if(nSigma==4) txt.DrawLatex(0.15,0.88,Form("Si%d : |#Delta x|<1.25", ich) );
	    }
	  else
	    {
	      txt.DrawLatex(0.15,0.88,Form("Si%d : Si%d>%d#sigma_{noise}(Si%d)",ich,otherSiPadCh,nSigma,otherSiPadCh));
	    }

	  c->SaveAs(outDir+Form("/mipsignal_%d_nsig%d_mode%d.png",ich,nSigma,mode));
	}
    }

}

//
void printSummary(TString outDir)
{
  ofstream fOut;
  fOut.open((outDir+"/summary.txt").Data());

  //print run summary
  fOut << "------------------------------------------------------------------------------------" << endl;
  fOut << "RUN " << runSummary.runNumber << endl;
  fOut << "\t Noise levels" << endl;
  for(std::map<Int_t, Float_t>::iterator it=runSummary.pedestals.begin(); it!=runSummary.pedestals.end(); it++)
    {
      fOut << "\t\t ch" << it->first << " res: " << it->second << " sigma: " << runSummary.pedestalsSigma[it->first] << endl;
    }
  fOut << "\t Inclusive beamspot from wire chambers " << endl;
  for(std::map<Int_t, Float_t>::iterator it=runSummary.xinc005.begin(); it!=runSummary.xinc005.end(); it++)
    {
      fOut << "\t\t ch" << it->first 
	   << " x: " << runSummary.xinc005[it->first] << " < " << runSummary.xinc050[it->first] << " < " << runSummary.xinc095[it->first] 
	   << " y: " << runSummary.yinc005[it->first] << " < " << runSummary.yinc050[it->first] << " < " << runSummary.yinc095[it->first] << endl;
    }
  fOut << "\t Fiducial beamspot from wire chambers (>2sigma in Si)" << endl;
  for(std::map<Int_t, Float_t>::iterator it=runSummary.x005.begin(); it!=runSummary.x005.end(); it++)
    {
      fOut << "\t\t ch" << it->first 
	   << " x: " << runSummary.x005[it->first] << " < " << runSummary.x050[it->first] << " < " << runSummary.x095[it->first] 
	   << " y: " << runSummary.y005[it->first] << " < " << runSummary.y050[it->first] << " < " << runSummary.y095[it->first] << endl;
    }
  fOut << "\t MIP peak estimates" << endl;
  for(std::map<Int_t, Float_t>::iterator it=runSummary.mip1.begin(); it!= runSummary.mip1.end(); it++)
    {
      fOut << "\t\t Si #" << it->first
	   << " Si#1:" << runSummary.mip1[it->first] << " +/- " << runSummary.mip1Unc[it->first]  
	   << " Si#2:" << runSummary.mip2[it->first] << " +/- " << runSummary.mip2Unc[it->first] << endl ;
    }
  fOut << "------------------------------------------------------------------------------------" << endl;

  fOut.close();
}

//
RunSummary_t analyzeRun(TString url="root://eoscms//eos/cms/store/cmst3/group/hgcal/TimingTB_H2_Jul2015/RECO/3eb93c8/RECO_3373.root",TString outDir="~/www/X0PbRuns",Int_t fitMode=CORRELATEDNOISE,TString siChargeEst="charge_integ_smallw_mcp")
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111111);
  gStyle->SetStatW(0.20);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.9);
  gStyle->SetStatH(0.15);
  gROOT->SetBatch(true);

  //open file
  TFile *fIn = TFile::Open(url);
  TTree* H4treeReco=(TTree*)fIn->Get("H4treeReco");
  UInt_t runNumber;
  H4treeReco->SetBranchAddress("runNumber",    &runNumber);
  H4treeReco->GetEntry(0);
  runSummary.runNumber=runNumber;

  outDir+=Form("/Run%d/",runSummary.runNumber);
  gSystem->Exec("mkdir -p "+outDir);

  TCanvas* c=new TCanvas("c","c",500,500);
  c->SetTopMargin(0.02);
  c->SetRightMargin(0.02);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);

  runSummary.pedestals.clear();
  runSummary.pedestalsSigma.clear();
  runSummary.x005.clear();
  runSummary.x050.clear();
  runSummary.x095.clear();
  runSummary.y005.clear();
  runSummary.y050.clear();
  runSummary.y095.clear();
  runSummary.xinc005.clear();
  runSummary.xinc050.clear();
  runSummary.xinc095.clear();
  runSummary.yinc005.clear();
  runSummary.yinc050.clear();
  runSummary.yinc095.clear();

  //pedestals (2nd time computes residuals)
  doPedestals(H4treeReco,c,outDir,siChargeEst);
  doPedestals(H4treeReco,c,outDir,siChargeEst);

  //beamspot
  determineFiducialBeamSpotFromWireChambers(H4treeReco,c,outDir,false,siChargeEst);
  determineFiducialBeamSpotFromWireChambers(H4treeReco,c,outDir,true,siChargeEst);

  //MIP fits
  doMIPfits(H4treeReco,c,outDir,0,fitMode,siChargeEst);
  doMIPfits(H4treeReco,c,outDir,1,fitMode,siChargeEst);

  printSummary(outDir);

  fIn->Close();

  //all done
  return runSummary;
}

//wrapper for all the MIP runs
void makeRunSummary(Int_t fitMode=FIXEDNOISE,Int_t saveMIPFor=2,TString outDir="~/www/HGCal/FastTimeTB/MIPCalibrationRuns",TString siChargeEst="charge_integ_smallw_mcp")
{
  if(fitMode==FIXEDNOISE)       outDir += "_fixednoise";
  if(fitMode==CORRELATEDNOISE)  outDir += "_correlatednoise";
  if(fitMode==INDEPENDENTNOISE) outDir += "_indnoise";
  outDir += "_mipAt"; outDir+=saveMIPFor;
  outDir += "_";      outDir+=siChargeEst;
  gSystem->Exec("mkdir -p " +outDir);

  TFile *fIn=TFile::Open(outDir+"/summary.root","RECREATE");
  TNtuple *ntuple = new TNtuple("summary",
				"summary",
				"run:xcen:xcenUnc:ycen:ycenUnc:noise_Si1:noiseUnc_Si1:noiseSigma_Si1:noiseSigmaUnc_Si1:noise_Si2:noiseUnc_Si2:noiseSigma_Si2:noiseSigmaUnc_Si2:mip_Si1:mipUnc_Si1:chi2ndof_Si1:mip_Si2:mipUnc_Si2:chi2ndof_Si2");
  ntuple->SetDirectory(0);

  Int_t runs[]    = {3370,3371,3372,3373,3374,3375,3376,3363,3369,3346};
  for(size_t i=0; i<sizeof(runs)/sizeof(Int_t); i++)
    {
      TString url("root://eoscms//eos/cms/store/cmst3/group/hgcal/TimingTB_H2_Jul2015/RECO/3eb93c8/RECO_");
      url += runs[i];
      url += ".root";
      analyzeRun(url,outDir,fitMode,siChargeEst);
      Float_t vars[]={(Float_t) runs[i],
		      runSummary.x050[1],
		      (Float_t) 0.5*fabs(runSummary.x095[1]-runSummary.x005[1]),
		      runSummary.y050[1],
		      (Float_t )0.5*fabs(runSummary.y095[1]-runSummary.y005[1]),
		      runSummary.pedestals[2],
		      runSummary.pedestalsUnc[2],
		      runSummary.pedestalsSigma[2],
		      runSummary.pedestalsSigmaUnc[2],
		      runSummary.pedestals[3],
		      runSummary.pedestalsUnc[3],
		      runSummary.pedestalsSigma[3],
		      runSummary.pedestalsSigmaUnc[3],
		      runSummary.mip1[saveMIPFor],
		      runSummary.mip1Unc[saveMIPFor],
		      runSummary.chi2ndof1[saveMIPFor],
		      runSummary.mip2[saveMIPFor],
		      runSummary.mip2Unc[saveMIPFor],
		      runSummary.chi2ndof2[saveMIPFor]
      };
      ntuple->Fill(vars);
    }

  fIn->cd();
  ntuple->Write();
  fIn->Close();
}


