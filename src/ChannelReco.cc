#include "interface/ChannelReco.h"
#include <iostream>
//
ChannelReco::ChannelReco(JSONWrapper::Object &cfg, PlotType ptype, Bool_t doPlot, Bool_t doProfile) :
  ptype_(ptype),
  doPlot_(doPlot),
  doProfile_(doProfile),
  plot_(0),
  waveform_(0)
{
  //read configuration
  name_              = cfg["name"].toString();
  group_             = cfg["group"].toInt();
  module_            = cfg["channel"].toInt();
  searchWindowUp_       = cfg["searchWindowUp"].toInt();
  searchWindowLo_       = cfg["searchWindowLo"].toInt();
  baselineSlope_       = cfg["baselineSlope"].toDouble();
  searchWindowTriggerRelative_ = cfg["searchWindowTriggerRelative"].toInt();
  searchWindowAfterUp_       = cfg["searchWindowAfterUp"].toDouble();
  searchWindowAfterLo_       = cfg["searchWindowAfterLo"].toDouble();
  chargeWindowLo_       = cfg["chargeWindowLo"].toDouble();
  chargeWindowUp_       = cfg["chargeWindowUp"].toDouble();
  absoluteTimeDelta_       = cfg["absoluteTimeDelta"].toDouble();
  mcpTimeDelta_       = cfg["mcpTimeDelta"].toDouble();
  spyWindowUp_       = cfg["spyWindowUp"].toInt();
  spyWindowLo_       = cfg["spyWindowLo"].toInt();
  cfdWindowLo_       = cfg["cfdWindowLo"].toDouble();
  thrForTiming_ = cfg["thrForTiming"].toDouble();
  pedestalWindowLo_   = cfg["pedestalWindowLo"].toInt();
  pedestalWindowUp_   = cfg["pedestalWindowUp"].toInt();
  thrForPulseInversion_ = cfg["thrForPulseInversion"].toInt();
  samplesToInterpolateAtMax_ = cfg["samplesToInterpolateAtMax"].toInt();
  samplesToInterpolateForCFD_ = cfg["samplesToInterpolateForCFD"].toInt();
  samplesToInterpolateForTD_ = cfg["samplesToInterpolateForTD"].toInt();
  smallChargeWindowSize_ = cfg["smallChargeWindowSize"].toInt();
  largeChargeWindowSize_ = cfg["largeChargeWindowSize"].toInt();
  //save configuration values in a histogram
  chRecoCfgH_ = new TH1F(name_+"_cfg",";"+name_+";Value",25,0,25);
  chRecoCfgH_->SetDirectory(0);
  chRecoCfgH_->SetBinContent(1,group_);                  chRecoCfgH_->GetXaxis()->SetBinLabel(1,"group");
  chRecoCfgH_->SetBinContent(2,module_);                 chRecoCfgH_->GetXaxis()->SetBinLabel(2,"channel");
  chRecoCfgH_->SetBinContent(3,searchWindowUp_);         chRecoCfgH_->GetXaxis()->SetBinLabel(3,"searchWindowUp");
  chRecoCfgH_->SetBinContent(4,searchWindowLo_);         chRecoCfgH_->GetXaxis()->SetBinLabel(4,"searchWindowLo");
  chRecoCfgH_->SetBinContent(5,searchWindowAfterUp_);         chRecoCfgH_->GetXaxis()->SetBinLabel(5,"searchWindowAfterUp");
  chRecoCfgH_->SetBinContent(6,searchWindowAfterLo_);         chRecoCfgH_->GetXaxis()->SetBinLabel(6,"searchWindowAfterLo");
  chRecoCfgH_->SetBinContent(7,cfdWindowLo_);            chRecoCfgH_->GetXaxis()->SetBinLabel(7,"cfdWindowLo");
  chRecoCfgH_->SetBinContent(8,thrForTiming_);           chRecoCfgH_->GetXaxis()->SetBinLabel(8,"thrForTiming");
  chRecoCfgH_->SetBinContent(9,pedestalWindowLo_);       chRecoCfgH_->GetXaxis()->SetBinLabel(9,"pedestalWindowLo");
  chRecoCfgH_->SetBinContent(10,pedestalWindowUp_);       chRecoCfgH_->GetXaxis()->SetBinLabel(10,"pedestalWindowUp");
  chRecoCfgH_->SetBinContent(11,thrForPulseInversion_);   chRecoCfgH_->GetXaxis()->SetBinLabel(11,"thrForPulseInversion");
  chRecoCfgH_->SetBinContent(12,spyWindowUp_);           chRecoCfgH_->GetXaxis()->SetBinLabel(12,"spyWindowUp");
  chRecoCfgH_->SetBinContent(13,spyWindowLo_);           chRecoCfgH_->GetXaxis()->SetBinLabel(13,"spyWindowLo");
  chRecoCfgH_->SetBinContent(14,samplesToInterpolateAtMax_);       chRecoCfgH_->GetXaxis()->SetBinLabel(14,"samplesToInterpolateAtMax");
  chRecoCfgH_->SetBinContent(15,samplesToInterpolateForCFD_);       chRecoCfgH_->GetXaxis()->SetBinLabel(15,"samplesToInterpolateForMax");
  chRecoCfgH_->SetBinContent(16,samplesToInterpolateForTD_);       chRecoCfgH_->GetXaxis()->SetBinLabel(16,"samplesToInterpolateForTD");
  chRecoCfgH_->SetBinContent(17,searchWindowTriggerRelative_);       chRecoCfgH_->GetXaxis()->SetBinLabel(17,"searchWindowTriggerRelative");
  chRecoCfgH_->SetBinContent(18,chargeWindowLo_);       chRecoCfgH_->GetXaxis()->SetBinLabel(18,"chargeWindowLo");
  chRecoCfgH_->SetBinContent(19,chargeWindowLo_);       chRecoCfgH_->GetXaxis()->SetBinLabel(19,"chargeWindowLo");
  chRecoCfgH_->SetBinContent(20,absoluteTimeDelta_);       chRecoCfgH_->GetXaxis()->SetBinLabel(20,"absoluteTimeDelta");
  chRecoCfgH_->SetBinContent(21,smallChargeWindowSize_);       chRecoCfgH_->GetXaxis()->SetBinLabel(21,"smallChargeWindowSize");
  chRecoCfgH_->SetBinContent(22,largeChargeWindowSize_);       chRecoCfgH_->GetXaxis()->SetBinLabel(22,"largeChargeWindowSize");
  chRecoCfgH_->SetBinContent(23,mcpTimeDelta_);       chRecoCfgH_->GetXaxis()->SetBinLabel(23,"mcpTimeDelta");
  chRecoCfgH_->SetBinContent(24,baselineSlope_);       chRecoCfgH_->GetXaxis()->SetBinLabel(24,"baselineSlope");

  //  std::cout << "\t " << name_ << " will be reconstructed from group=" << group_ << " channel=" << module_ << std::endl;
  Print();
}

void ChannelReco::Print()
{
  std::cout << "=====  Dump Channel Config  =====" << std::endl;
  std::cout << "name\t" << name_ << std::endl;
  std::cout << "group\t" << group_ << std::endl;
  std::cout << "channel\t" << module_ << std::endl;
  std::cout << "Pedestal WindowLo\t" << this->GetPedestalWindowLo() << std::endl;
  std::cout << "Pedestal WindowUp\t" << this->GetPedestalWindowUp() << std::endl;
  std::cout << "Inversion Threshold\t" << this->GetThrForPulseInversion() << std::endl;
  std::cout << "Baseline Slope\t" << this->GetBaselineSlope() << std::endl;
  std::cout << "Max Search WindowLo\t" << this->GetSearchWindowLo() << std::endl;
  std::cout << "Max Search WindowUp\t" << this->GetSearchWindowUp() << std::endl;
  std::cout << "Max Search TriggerRelative\t" << this->GetSearchWindowTriggerRelative() << std::endl;
  std::cout << "Max Search WindowAfterLo\t" << this->GetSearchWindowAfterLo() << std::endl;
  std::cout << "Max Search WindowAfterUp\t" << this->GetSearchWindowAfterUp() << std::endl;
  std::cout << "SaveWaveform WindowLo\t" << this->GetSpyWindowLo() << std::endl;
  std::cout << "SaveWaveform WindowUp\t" << this->GetSpyWindowUp() << std::endl;
  std::cout << "CFD WindowLo\t" << this->GetCFDWindowLo() << std::endl;
  std::cout << "TD Value\t" << this->GetThrForTiming() << std::endl;
  std::cout << "Samples to interpolate at MAX\t" << this->GetSamplesToInterpolateAtMax() << std::endl;
  std::cout << "Samples to interpolate for CFD\t" << this->GetSamplesToInterpolateForCFD() << std::endl;
  std::cout << "Samples to interpolate for TD\t" << this->GetSamplesToInterpolateForTD() << std::endl;
  std::cout << "Charge WindowLo\t" << this->GetChargeWindowLo() << std::endl;
  std::cout << "Charge WindowUp\t" << this->GetChargeWindowUp() << std::endl;
  std::cout << "Small Charge WindowSize\t" << this->GetSmallChargeWindowsSize() << std::endl;
  std::cout << "Large Charge WindowSize\t" << this->GetLargeChargeWindowsSize() << std::endl;
  std::cout << "Absolute Delta Time (wrt Trigger)\t" << this->GetAbsoluteTimeDelta() << std::endl;
  std::cout << "Absolute Delta Time (wrt MCP)\t" << this->GetMCPTimeDelta() << std::endl;
}
//
void ChannelReco::Fill(float valX, float valY, int i) 
{
  if(plot_==0 || !doPlot_)
    {
      return;
    }
  if (ptype_!=kPlot1D && ptype_!=kPlotGraph && ptype_!=kPlot1DProf) 
    {
      std::cout << "[ChannelReco::Fill] caught wrong pptype_ in " << name_.Data() << std::endl; 
      return;
    }

  if (ptype_==kPlot1D)
    {
      if (i<0) x_.push_back(valX);
      else 
	{
	  if (x_.size()<=i) x_.resize(i+1);
	  x_.at(i)=valX;
	}
      (dynamic_cast<TH1F*>(plot_))->Fill(valX);
    }
  if (ptype_==kPlot1DProf)
    {
      if (i<0) { x_.push_back(valX); y_.push_back(valY); }
      else
	{
	  if (x_.size()<=i)
	    {
	      x_.resize(i+1);
	      y_.resize(i+1);
	    }
	  x_.at(i)=valX;
	  y_.at(i)=valY;
	}
      (dynamic_cast<TProfile*>(plot_))->Fill(valX,valY);
    }
  if (ptype_==kPlotGraph)
    {
      (dynamic_cast<TGraph*>(plot_))->SetPoint(group_,module_,valX);
    }
}

//
void ChannelReco::Fill2D(float valX, float valY, float valZ, int i) 
{
  if(plot_==0 || !doPlot_)
    {
      return;
    }
  if (ptype_!=kPlot2D && ptype_!=kPlot2DProf) 
    {
      std::cout << "[ChannelReco::Fill2D] caught wrong pptype_ in " << name_.Data() << std::endl; 
      return;
    }

  if (ptype_==kPlot2D)
    {
      if (i<0) {x_.push_back(valX); y_.push_back(valY);}
      else
	{
	  if (x_.size()<=i) 
	    {
	      x_.resize(i+1);
	      y_.resize(i+1);
	    }
	  x_.at(i)=valX;
	  y_.at(i)=valY;
	}
      (dynamic_cast<TH2F*>(plot_))->Fill(valX,valY);
    }
  if (ptype_==kPlot2DProf)
    {
      if (i<0)
	{
	  x_.push_back(valX); 
	  y_.push_back(valY); 
	  z_.push_back(valZ);
	}
      else {
	if (x_.size()<=i) 
	  {
	    x_.resize(i+1);
	    y_.resize(i+1);
	    z_.resize(i+1);
	  }
	x_.at(i)=valX;
	y_.at(i)=valY;
	z_.at(i)=valZ;
      }
      (dynamic_cast<TProfile2D*>(plot_))->Fill(valX,valY,valZ);
    }
}


//
ChannelReco::~ChannelReco()
{
  if (plot_)     delete plot_;
  if (waveform_) delete waveform_;
}
