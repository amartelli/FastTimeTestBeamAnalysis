#include "interface/histoFuncT.h"
#include "TF1.h"

#include <iostream>

histoFuncT::histoFuncT(TH1F* histo, int& tStart, int& tStop):
  tStart_p(tStart), tStop_p(tStop)
{
  histo_p = histo;

  //  std::cout << " fitTemplate defined region " << histo_p->GetBinCenter(1) << "  " << histo_p->GetBinCenter(histo_p->GetNbinsX()) << std::endl;
  //  std::cout << " tStart_p = " << tStart_p << " tStop_p = " << tStop_p <<  std::endl; 
  // int& tStart, int& tStop in case of saturated pulses
}


histoFuncT::~histoFuncT(void)
{}


double histoFuncT::operator()(double* x, double* par){
  //  std::cout << " operator () " << std::endl;

  double xx = (x[0]);
  //    double xx = par[1]* (x[0] - par[2]);

  double xMin = histo_p->GetBinCenter(1);
  double xMax = histo_p->GetBinCenter(histo_p->GetNbinsX());

  //  std::cout << " xx = " << xx << " xMin = " << xMin << " xMax =  " << xMax << std::endl;

  //TF1 range outside fitTemplate
  if( (xx < xMin) || (xx >= xMax) ) return 1.e-10;

  // par[0] = y scale
  // par[1] = x scale
  // par[2] = x shift
  // par[3] = y shift

  //if saturated pulse can fit with template the pulses over threshold
  if(tStart_p != 0 && tStop_p != 0 && xx > tStart_p && (xx) < tStop_p){
    //    std::cout << " rejected xx = " << xx  << std::endl;                                                                                     
    TF1::RejectPoint();
    return 0;
  }
  else{
    //    std::cout << " xx = " << xx << std::endl;
    int bin = histo_p->FindBin(xx);
    int bin1 = 0;
    int bin2 = 0;
    
    if(xx >= histo_p->GetBinCenter(bin)){
      bin1 = bin;
      bin2 = bin+1;
    }    
    else{
      bin1 = bin-1;
      bin2 = bin;
    }

    //    std::cout << " bin1 = " << bin1 << " bin2 = " << bin2 << std::endl;
    
    //interpolated point for the fit
    // double y1 = histo_p->GetBinContent(bin1);
    // double y2 = histo_p->GetBinContent(bin2);
    // double x1T = histo_p->GetBinCenter(bin1 - par[2]);
    // double x2T = histo_p->GetBinCenter(bin2 - par[2]); 
    double x1 = histo_p->GetBinCenter(bin1);
    double x2 = histo_p->GetBinCenter(bin2);
    double y1T = histo_p->GetBinContent(par[1]*(bin1 - par[2]));
    double y2T = histo_p->GetBinContent(par[1]*(bin2 - par[2]));
    
    // std::cout << " x1 = " << x1 << " x2 = " << x2 
    // 	      << " y1T = " << y1T << " y2T = " << y2T << std::endl;

    double m = 1. * (y2T - y1T) / par[1] / (x2 - x1);
    //    double m = 1. * (y2 - y1) / (x2 - x1);
    
    /* 
    if( (y1 + m * (xx - x1)) < 1.e-10)                                                                                                              
      return 1.e-10;                                                                                                                                
    */

    return par[0] * (y1T + m * par[1] *(xx - x1));
  } 


  return 1.e-10;
}
