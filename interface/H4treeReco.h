#ifndef H4treeReco_h
#define H4treeReco_h

#include "interface/H4tree.h"
#include "interface/ChannelPlot.h"
#include "interface/JSONWrapper.h"

#include "TFile.h"
#include "TString.h"

#include <set>

class H4treeReco : public H4tree
{

 public:
  H4treeReco(TChain *,JSONWrapper::Object *cfg,TString outUrl="H4treeRecoOut.root");
  void Loop(); 
  void FillTDC();
  void FillWaveforms();
  ~H4treeReco();
  
 private:

  void InitDigi(JSONWrapper::Object *cfg);

  typedef std::pair<UInt_t,UInt_t>         GroupChannelKey_t;
  std::map<GroupChannelKey_t,ChannelPlot*> chPlots_;

  //TDC readings
  UInt_t wcXl_, wcXr_, wcYd_, wcYu_;
  UInt_t MaxTdcChannels_,MaxTdcReadings_;
  std::vector< std::vector<Float_t> > tdc_readings_;
  float tdc_recox_, tdc_recoy_;

  UInt_t nActiveDigitizerChannels_;
  UInt_t maxch_;
  Float_t group_[100],ch_[100];
  Float_t pedestal_[100],         pedestalRMS_[100];
  Float_t wave_max_[100];
  Float_t charge_integ_max_[100], charge_integ_max30_[100], charge_integ_max50_[100], charge_integ_[100];
  Float_t t_max_[100],            t_max_frac30_[100],       t_max_frac50_[100];

  TTree *recoT_;
  TH1F *groupNamesH_, *channelNamesH_;
  TFile *fOut_;
};

#endif
