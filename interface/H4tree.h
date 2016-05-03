//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul 15 23:32:17 2015 by ROOT version 5.34/30
// from TTree H4tree/H4 testbeam tree
// found on file: root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/TimingTB_H2_Jul2015/raw/DataTree/3374/14.root
//////////////////////////////////////////////////////////

#ifndef H4tree_h
#define H4tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class H4tree {
public :
   TChain         *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          runNumber;
   UInt_t          spillNumber;
   UInt_t          evtNumber;
   UInt_t          evtTimeDist;
   UInt_t          evtTimeStart;
   UInt_t          nEvtTimes;
   ULong64_t       evtTime[1];   //[nEvtTimes]
   UInt_t          evtTimeBoard[1];   //[nEvtTimes]
   UInt_t          nAdcChannels;
   UInt_t          adcBoard[32];   //[nAdcChannels]
   UInt_t          adcChannel[32];   //[nAdcChannels]
   UInt_t          adcData[32];   //[nAdcChannels]
   UInt_t          nTdcChannels;
   UInt_t          tdcBoard[22];   //[nTdcChannels]
   UInt_t          tdcChannel[22];   //[nTdcChannels]
   UInt_t          tdcData[22];   //[nTdcChannels]
   UInt_t          nDigiSamples;
   UInt_t          digiFrequency[36865];   //[nDigiSamples]
   UInt_t          digiGroup[36865];   //[nDigiSamples]
   UInt_t          digiChannel[36865];   //[nDigiSamples]
   UInt_t          digiSampleIndex[36865];   //[nDigiSamples]
   Float_t         digiSampleValue[36865];   //[nDigiSamples]
   UInt_t          digiBoard[36865];   //[nDigiSamples]
   UInt_t          nScalerWords;
   UInt_t          scalerWord[8];   //[nScalerWords]
   UInt_t          scalerBoard[8];   //[nScalerWords]
   UInt_t          nPatterns;
   UInt_t          pattern[1];   //[nPatterns]
   UInt_t          patternBoard[1];   //[nPatterns]
   UInt_t          patternChannel[1];   //[nPatterns]
   UInt_t          nTriggerWords;
   UInt_t          triggerWords[1];   //[nTriggerWords]
   UInt_t          triggerWordsBoard[1];   //[nTriggerWords]

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_spillNumber;   //!
   TBranch        *b_evtNumber;   //!
   TBranch        *b_evtTimeDist;   //!
   TBranch        *b_evtTimeStart;   //!
   TBranch        *b_nEvtTimes;   //!
   TBranch        *b_evtTime;   //!
   TBranch        *b_evtTimeBoard;   //!
   TBranch        *b_nAdcChannels;   //!
   TBranch        *b_adcBoard;   //!
   TBranch        *b_adcChannel;   //!
   TBranch        *b_adcData;   //!
   TBranch        *b_nTdcChannels;   //!
   TBranch        *b_tdcBoard;   //!
   TBranch        *b_tdcChannel;   //!
   TBranch        *b_tdcData;   //!
   TBranch        *b_nDigiSamples;   //!
   TBranch        *b_digiFrequency;   //!
   TBranch        *b_digiGroup;   //!
   TBranch        *b_digiChannel;   //!
   TBranch        *b_digiSampleIndex;   //!
   TBranch        *b_digiSampleValue;   //!
   TBranch        *b_digiBoard;   //!
   TBranch        *b_nScalerWords;   //!
   TBranch        *b_scalerWord;   //!
   TBranch        *b_scalerBoard;   //!
   TBranch        *b_nPatterns;   //!
   TBranch        *b_pattern;   //!
   TBranch        *b_patternBoard;   //!
   TBranch        *b_patternChannel;   //!
   TBranch        *b_nTriggerWords;   //!
   TBranch        *b_triggerWords;   //!
   TBranch        *b_triggerWordsBoard;   //!

   H4tree(TChain *chain=0);
   virtual ~H4tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
