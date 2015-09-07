#ifndef _varplot_h_
#define _varplot_h_

#include "TObject.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TString.h"

#include "interface/WaveformUtils.hpp"
#include "interface/WaveformFit.hpp"
#include "interface/JSONWrapper.h"

class ChannelReco{

public:

  enum PlotType { kNull, kPlot1D, kPlot1DProf, kPlot2D, kPlot2DProf, kPlotGraph, kPlotHistory };

  ChannelReco(JSONWrapper::Object &cfg, PlotType ptype, Bool_t doPlot, Bool_t doProfile);
 
  UInt_t Size()                                               { return x_.size(); }
  std::vector<float> &GetX()                                  { return x_; }
  float GetX(UInt_t i)                                        { return x_.at(i); }
  std::vector<float> &GetY()                                  { return y_; }
  float GetY(UInt_t i)                                        { return y_.at(i); }
  std::pair<float,float> GetXY(uint i)                        { return std::pair<float,float>(x_.at(i),y_.at(i)); }
  std::pair<std::vector<float> *, std::vector<float>*> Get2D(){ return std::pair<std::vector<float>*, std::vector<float>*>(&x_,&y_); }
  void ClearVectors()                                         { x_.clear(); y_.clear(); }
  void SetPlot(TObject* plot)                                 { plot_=plot; }
  TObject* GetPlot()                                          { return plot_; }
  TString GetName()                                           { return name_; }
  //Int
  Int_t GetGroup()                                            { return group_; }
  Int_t GetModule()                                           { return module_; }
  Int_t GetThrForPulseInversion()                             { return thrForPulseInversion_; }
  Int_t GetPedestalWindowLo()                                 { return pedestalWindowLo_; }
  Int_t GetPedestalWindowUp()                                 { return pedestalWindowUp_; }
  Int_t GetSearchWindowUp()                                   { return searchWindowUp_; }
  Int_t GetSearchWindowLo()                                   { return searchWindowLo_; }  
  Int_t GetSearchWindowTriggerRelative()                      { return searchWindowTriggerRelative_; }  
  Int_t GetSpyWindowLo()                                      { return spyWindowLo_; }
  Int_t GetSpyWindowUp()                                      { return spyWindowUp_; }
  Int_t GetSamplesToInterpolateAtMax()                        { return samplesToInterpolateAtMax_; } 
  Int_t GetSamplesToInterpolateForCFD()                       { return samplesToInterpolateForCFD_; } 
  Int_t GetSamplesToInterpolateForTD()                       { return samplesToInterpolateForTD_; } 
  Int_t GetSmallChargeWindowsSize()                           { return smallChargeWindowSize_; };
  Int_t GetLargeChargeWindowsSize()                           { return largeChargeWindowSize_; };
  //Float
  Float_t GetCFDWindowLo()                                    { return cfdWindowLo_; }
  Float_t GetBaselineSlope()                                  { return baselineSlope_; }
  Float_t GetSearchWindowAfterLo()                            { return searchWindowAfterLo_; }
  Float_t GetSearchWindowAfterUp()                            { return searchWindowAfterUp_; }
  Float_t GetChargeWindowLo()                                 { return chargeWindowLo_; };                
  Float_t GetChargeWindowUp()                                 { return chargeWindowUp_; };                
  Float_t GetAbsoluteTimeDelta()                              { return absoluteTimeDelta_; };                
  Float_t GetThrForTiming()                                   { return thrForTiming_; }
  Float_t GetMCPTimeDelta()                                   { return mcpTimeDelta_; }
  Waveform *GetWaveform()                                     { return waveform_; }
  TH1F *GetConfigSummary()                                    { return chRecoCfgH_; }
  void SetWaveform(Waveform *w)                               { waveform_=w; }
  void Fill(float valX, float valY, int i=-1);
  void Fill2D(float valX, float valY, float valZ, int i);
  void Print();
  
  ~ChannelReco();

 private:
  
  Int_t group_, module_,searchWindowUp_,searchWindowLo_,spyWindowLo_,spyWindowUp_;
  Float_t cfdWindowLo_;
  Float_t baselineSlope_;
  Float_t searchWindowAfterLo_, searchWindowAfterUp_;
  Float_t chargeWindowUp_, chargeWindowLo_, absoluteTimeDelta_;
  Int_t smallChargeWindowSize_;
  Int_t largeChargeWindowSize_;
  Int_t searchWindowTriggerRelative_;
  Int_t thrForPulseInversion_,pedestalWindowLo_,pedestalWindowUp_;
  Int_t samplesToInterpolateAtMax_,samplesToInterpolateForCFD_,samplesToInterpolateForTD_;
  Float_t thrForTiming_;
  Float_t mcpTimeDelta_;
  int ptype_;
  bool doPlot_, doProfile_;
  TString name_;
  std::vector<float> x_,y_,z_;
  TH1F *chRecoCfgH_;
  TObject *plot_;
  Waveform *waveform_;
};

#endif
